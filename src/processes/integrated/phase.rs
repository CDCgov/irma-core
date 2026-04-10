//! Phase clustering assignment.
//!
//! Reads in the variants table and associated distance matrix for a gene.
//! Assigns a phase number to each variant by creating a single-linkage
//! agglomerative clustering tree from the distance matrix and cutting the
//! branches at a specified height to provide phase groups. If minor variants
//! are phased, they will belong to the same group and have the same phase
//! number. If variants are not in-phase, they will have different phase
//! numbers. If no variants are phased, each minor variant will have its own
//! unique phase number. Phase numbers are added to the variants table written
//! at the end of an IRMA assembly.

use crate::io::{InputOptions, OutputOptions};
use clap::Args;
use kodama::{Method, linkage};
use std::{
    collections::{HashMap, hash_map::Entry},
    fmt::Display,
    io::{BufRead, Write},
    path::PathBuf,
    str::Split,
};
use zoe::{data::err::ResultWithErrorContext, iter_utils::ProcessResultsExt};

#[derive(Args, Debug)]
pub struct PhaseArgs {
    /// Gene variants file `<gene>-variants.txt`
    pub variants_file: PathBuf,
    /// Gene variants square matrix file `<gene>-EXPENRD.sqm`
    pub sqm_file:      PathBuf,
    /// Optional height value for phase clustering. The default height was
    /// selected empirically on known mixture data to render 95% sensitivity and
    /// specificity.
    #[arg(short, long, value_parser = validate_height, default_value_t = 0.78)]
    pub tree_height:   f64,
}

/// Location of the "Position" column in the `variants_file`, 0-indexed
const POSITION_COLUMN: usize = 1;
/// Location of the "Minority_Allele" column in the `variants_file`, 0-indexed
const MINORITY_ALLELE_COLUMN: usize = 4;
// TODO: Replace with max operator when it becomes const compatible
/// Minimum number of columns expected in the `variants_file`
const MIN_COLUMNS: usize = POSITION_COLUMN + MINORITY_ALLELE_COLUMN.saturating_sub(POSITION_COLUMN) + 1;

pub fn phase_process(args: PhaseArgs) -> std::io::Result<()> {
    let mut variants_file_lines = InputOptions::new_from_path(&args.variants_file).use_file().open()?.lines();

    let Some(header) = variants_file_lines.next().transpose()? else {
        return Err(std::io::Error::other(format!(
            "File is empty: '{}'",
            args.variants_file.display()
        )));
    };
    validate_header(&header).with_file_context("Failed to validate header from variants file", &args.variants_file)?;

    let mut variants_file_table = Vec::new();
    for (line_ind, line) in variants_file_lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let variants_file_line = VariantsFileLine::try_from(line).with_file_context(
            format!(
                "Failed to parse line number {line_num} from variants file",
                line_num = line_ind + 2
            ),
            &args.variants_file,
        )?;

        variants_file_table.push(variants_file_line);
    }

    let mut variants_file_writer = OutputOptions::new_from_path(&args.variants_file).use_file().open()?;

    writeln!(variants_file_writer, "{header}\tPhase", header = header.trim())?;

    // If there is only a single variant, skip clustering calculation. Single
    // variant is assigned phase number `1`.
    if let [single_row] = variants_file_table.as_slice() {
        writeln!(variants_file_writer, "{single_row}\t1")?;
    } else if !variants_file_table.is_empty() {
        let variants_matrix_reader = InputOptions::new_from_path(&args.sqm_file).use_file().open()?;
        // Phase clustering calculation and assignment happens here.
        let variants_matrix = variants_matrix_reader.lines().process_results(|lines| {
            VariantsMatrix::from_sqm_file(lines, args.tree_height, variants_file_table.len())
                .with_file_context("Cannot parse the .sqm file", &args.sqm_file)
        })??;

        for line in variants_file_table {
            let Some(&var_mat_ind) = variants_matrix.tag_index.get(&(line.position, line.minority_allele)) else {
                return Err(std::io::Error::other(format!(
                    "Cannot find variant position \"{position}\" with minority allele \"{min_allele}\" in file: '{sqm_file}'",
                    position = line.position,
                    min_allele = line.minority_allele as char,
                    sqm_file = args.sqm_file.display()
                )));
            };
            let phase_num = variants_matrix.phase_group[var_mat_ind];

            writeln!(variants_file_writer, "{line}\t{phase_num}")?
        }
    }

    variants_file_writer.flush()
}

/// Validates the type and range of the `tree_height` argument.
fn validate_height(tree_height: &str) -> Result<f64, String> {
    match tree_height.parse::<f64>() {
        Ok(c) if (0.0..=1.0).contains(&c) => Ok(c),
        Ok(_) => Err("Value must be between 0.0 and 1.0".to_string()),
        Err(e) => Err(format!("{e}")),
    }
}

/// Validates the header in the `<gene>-variants.txt` file.
///
/// Helps validate the same column location assumptions when parsing lines in
/// [`VariantsFileLine`] for the variants file table. Returns an error early if
/// any of the assumptions are false.
fn validate_header(input_header: &str) -> std::io::Result<()> {
    let header_vec = input_header.split('\t').collect::<Vec<_>>();

    if header_vec.len() < MIN_COLUMNS {
        Err(std::io::Error::other(format!(
            "Expected at least {MIN_COLUMNS} columns, found {cols}.",
            cols = header_vec.len()
        )))
    } else if header_vec[POSITION_COLUMN] != "Position" {
        Err(std::io::Error::other(format!(
            "Header has misplaced column labeled \"Position\". Expected column number {pos_col}.",
            pos_col = POSITION_COLUMN + 1
        )))
    } else if header_vec[MINORITY_ALLELE_COLUMN] != "Minority_Allele" {
        Err(std::io::Error::other(format!(
            "Header has misplaced column labeled \"Minority_Allele\". Expected column number {ma_col}.",
            ma_col = MINORITY_ALLELE_COLUMN + 1
        )))
    } else if header_vec.into_iter().any(|header| header == "Phase") {
        Err(std::io::Error::other("Input file already has a \"Phase\" column."))
    } else {
        Ok(())
    }
}

/// A structure that holds the information from one line of the `variants_file`.
struct VariantsFileLine {
    line:            String,
    position:        usize,
    minority_allele: u8,
}

impl TryFrom<String> for VariantsFileLine {
    type Error = std::io::Error;

    fn try_from(s: String) -> std::io::Result<Self> {
        let split_line = s.split('\t').collect::<Vec<_>>();

        if split_line.len() < MIN_COLUMNS {
            return Err(std::io::Error::other(format!(
                "Expected at least {MIN_COLUMNS} items, found {line_items}.",
                line_items = split_line.len()
            )));
        }

        let (position, minority_allele) = (split_line[POSITION_COLUMN], split_line[MINORITY_ALLELE_COLUMN].as_bytes());

        let &[minority_allele] = minority_allele else {
            return Err(std::io::Error::other(format!(
                "Minority Allele field does not contain a single character: \"{min_allele}\".",
                min_allele = split_line[MINORITY_ALLELE_COLUMN]
            )));
        };

        if !minority_allele.is_ascii_alphabetic() {
            return Err(std::io::Error::other(format!(
                "Minority Allele field is not ASCII alphabetic: \"{min_allele}\".",
                min_allele = split_line[MINORITY_ALLELE_COLUMN]
            )));
        }

        let position = position
            .parse::<usize>()
            .with_context(format!("Failed to parse variant position as integer: \"{position}\".",))?;

        Ok(VariantsFileLine {
            line: s,
            position,
            minority_allele,
        })
    }
}

impl Display for VariantsFileLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.line.trim())
    }
}

/// A structure that holds the information from the square matrix `sqm_file`.
struct VariantsMatrix {
    /// Distance matrix between variants in condensed format. Stores the values
    /// of the upper triangle of the matrix above the diagonal in a vector.
    condensed_matrix:   Vec<f64>,
    number_of_variants: usize,
    /// A map from the variant position and minority allele (tag) to the variant
    /// index.
    tag_index:          HashMap<(usize, u8), usize>,
    /// The phase group of each variant.
    phase_group:        Vec<usize>,
}

impl VariantsMatrix {
    /// Reads each line of the `sqm_file` and parses the values into the
    /// variants matrix.
    fn from_sqm_file<I>(lines: I, tree_height: f64, number_of_variants: usize) -> std::io::Result<Self>
    where
        I: Iterator<Item = String>, {
        let mut variants_matrix = VariantsMatrix {
            number_of_variants,
            condensed_matrix: Vec::with_capacity((number_of_variants) * (number_of_variants - 1) / 2),
            tag_index: HashMap::with_capacity(number_of_variants),
            phase_group: vec![0; number_of_variants],
        };

        let mut row_ind = 0;
        for row in lines {
            if row.trim().is_empty() {
                continue;
            }

            if row_ind + 1 > number_of_variants {
                return Err(std::io::Error::other(format!(
                    "The matrix is not square. There are {number_of_variants} columns, but there are more rows than that."
                )));
            }

            let mut row_iter = row.split('\t');
            let tag = row_iter.next().expect("Split is always non-empty");
            variants_matrix.assign_variant_tag(tag, row_ind).with_context(format!(
                "Failed to parse variant position with minority allele on line number {line_no}",
                line_no = row_ind + 1
            ))?;
            variants_matrix
                .parse_matrix_row(row_iter, row_ind)
                .with_context(format!("Failed to parse matrix row number {row_num}", row_num = row_ind + 1))?;

            row_ind += 1;
        }

        match row_ind {
            0 => return Err(std::io::Error::other("Found no variants in .sqm file.")),
            r if r < number_of_variants => {
                return Err(std::io::Error::other(format!(
                    "The matrix is not square. There are {number_of_variants} columns, but there are only {row_ind} rows."
                )));
            }
            _ => {}
        }

        variants_matrix.group_by_tree_height(tree_height);
        Ok(variants_matrix)
    }

    /// Parses the values of a single line in the `sqm_file`, and updates the
    /// [`VariantsMatrix`] with the values.
    ///
    /// This only parses the upper triangular entries of the matrix, ignoring the
    /// diagonal and entries below it.
    fn parse_matrix_row(&mut self, row_iter: Split<'_, char>, row_ind: usize) -> std::io::Result<()> {
        let mut row_len = 0;

        for (col_ind, val) in row_iter
            .enumerate()
            .inspect(|(col_ind, _)| row_len = *col_ind + 1)
            .skip(row_ind + 1)
        {
            let val = val.parse::<f64>().with_context(format!(
                "Failed to parse matrix value \"{val}\" at column number {col_num}.",
                col_num = col_ind + 1
            ))?;

            self.condensed_matrix.push(val);
        }

        if row_len != self.number_of_variants {
            if row_ind == 0 {
                Err(std::io::Error::other(format!(
                    "Incorrect number of variants found. This row contained {row_len}."
                )))
            } else {
                Err(std::io::Error::other(format!(
                    "Incorrect number of variants found. Row 1 contained {expected} values but this row contained {row_len}.",
                    expected = self.number_of_variants
                )))
            }
        } else {
            Ok(())
        }
    }

    /// Assigns the specified row index to the given variant (represented as a
    /// tag).
    ///
    /// The tag from the `.sqm` file should be of the form
    /// `<position><minority_allele>` where `<position>` is numeric digits and
    /// `<minority_allele>` is always a single character.
    fn assign_variant_tag(&mut self, tag: &str, row_ind: usize) -> std::io::Result<()> {
        let Some((&min_allele, position)) = tag.as_bytes().split_last() else {
            return Err(std::io::Error::other("Variant field is empty."));
        };

        if !min_allele.is_ascii_alphabetic() {
            return Err(std::io::Error::other(format!(
                "The first field must end in an alphabetic character as the minority allele. Found: \"{tag}\""
            )));
        }

        let position_str = str::from_utf8(position).with_context(format!(
            "Variant position and allele must contain ASCII data. Found: \"{tag}\""
        ))?;

        let position = position_str.parse::<usize>().with_context(format!(
            "Failed to parse variant position as integer. Found: \"{position_str}\""
        ))?;

        match self.tag_index.entry((position, min_allele)) {
            Entry::Vacant(entry) => entry.insert(row_ind),
            Entry::Occupied(_) => {
                return Err(std::io::Error::other(format!(
                    "Duplicate variant position and minority allele \"{tag}\" found."
                )));
            }
        };

        Ok(())
    }

    /// Creates a dendrogram and clusters phases according to a given
    /// dissimilarity threshold. Assigns the resulting phase groups to each
    /// variant.
    fn group_by_tree_height(&mut self, tree_height: f64) {
        let dendrogram = linkage(&mut self.condensed_matrix, self.number_of_variants, Method::Single);
        let mut tree = UnionFindTree::new(self.number_of_variants);

        for (step_number, step) in dendrogram
            .steps()
            .iter()
            .take_while(|step| step.dissimilarity <= tree_height)
            .enumerate()
        {
            tree.union(step.cluster1, step.cluster2, step.size, step_number);
        }

        self.assign_phase_groups(tree);
    }

    /// Sorts the phase clusters by size so the largest is labeled `1`, the
    /// second largest is labeled `2`, etc. Assigns phase clusters to each
    /// variant.
    fn assign_phase_groups(&mut self, tree: UnionFindTree) {
        let mut clusters_ordered = tree.cluster_sizes.iter().collect::<Vec<_>>();
        clusters_ordered.sort_by_key(|b| std::cmp::Reverse(b.1));

        let mut root_to_phase = HashMap::with_capacity(clusters_ordered.len());
        for (label, (&root, _)) in clusters_ordered.into_iter().enumerate() {
            root_to_phase.insert(root, label + 1);
        }

        // If not part of a cluster, label sequentially.
        let mut new_phase_number = root_to_phase.len() + 1;
        for i in 0..self.number_of_variants {
            match root_to_phase.get(&tree.find(i)) {
                Some(&phase) => self.phase_group[i] = phase,
                None => {
                    self.phase_group[i] = new_phase_number;
                    new_phase_number += 1;
                }
            }
        }
    }
}

/// A union-find tree structure for storing the phase clusters.
struct UnionFindTree {
    /// Stores the relationship between variants. `0..N` stores the roots of the
    /// variants, `N..2N-1` stores the internal connecting nodes.
    canonical_clusters: Vec<usize>,
    /// A map from the given variant root to the size of the cluster that root
    /// belongs to.
    cluster_sizes:      HashMap<usize, usize>,
    /// Number of outer leaves in the tree.
    number_of_variants: usize,
}

impl UnionFindTree {
    fn new(number_of_variants: usize) -> Self {
        Self {
            canonical_clusters: (0..2 * number_of_variants - 1).collect::<Vec<_>>(),
            cluster_sizes: HashMap::new(),
            number_of_variants,
        }
    }

    /// Takes the index of a variant and returns the index of the canonical
    /// cluster that variant belongs to.
    fn find(&self, mut cluster: usize) -> usize {
        while cluster != self.canonical_clusters[cluster] {
            cluster = self.canonical_clusters[cluster];
        }
        cluster
    }

    /// Takes two variant clusters and merges them into the same cluster. The
    /// initially larger cluster becomes the canonical cluster for both.
    fn union(&mut self, cluster1: usize, cluster2: usize, merged_size: usize, step_num: usize) {
        let parent_cluster1 = self.find(cluster1);
        let parent_cluster2 = self.find(cluster2);

        let parent_size1 = self.cluster_sizes.get(&parent_cluster1);
        let parent_size2 = self.cluster_sizes.get(&parent_cluster2);

        let canonical_cluster = match (parent_size1, parent_size2) {
            (Some(&size1), Some(&size2)) => {
                if size1 >= size2 {
                    self.cluster_sizes.remove(&parent_cluster2);
                    parent_cluster1
                } else {
                    self.cluster_sizes.remove(&parent_cluster1);
                    parent_cluster2
                }
            }
            (_, Some(_)) => parent_cluster2,
            _ => parent_cluster1,
        };

        let new_cluster_idx = step_num + self.number_of_variants;
        self.canonical_clusters[parent_cluster1] = canonical_cluster;
        self.canonical_clusters[parent_cluster2] = canonical_cluster;
        self.canonical_clusters[new_cluster_idx] = canonical_cluster;
        self.cluster_sizes.insert(canonical_cluster, merged_size);
    }
}
