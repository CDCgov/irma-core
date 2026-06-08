//! Phase clustering assignment.
//!
//! Reads in the variants table and associated distance matrix for a gene.
//! Assigns a phase number to each variant by clustering variants connected by
//! distances at or below a specified height. The variants table is rewritten
//! with a trailing `Phase` column. If there is a single variant, it is assigned
//! phase number `1` without reading the matrix.

use clap::Args;
use irma_records::io::{InputOptions, OutputOptions};
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
    validate_header(&header).with_path_context("Failed to validate header from variants file", &args.variants_file)?;

    let mut variants_file_table = Vec::new();
    for (line_ind, line) in variants_file_lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let variants_file_line = VariantsFileLine::try_from(line).with_path_context(
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
        let variant_phases = variants_matrix_reader.lines().process_results(|lines| {
            VariantPhases::from_sqm_file(lines, args.tree_height, variants_file_table.len())
                .with_path_context("Cannot parse the .sqm file", &args.sqm_file)
        })??;

        for line in variants_file_table {
            let Some(&phase_num) = variant_phases.phase_by_tag.get(&VariantTag {
                position:   line.position,
                min_allele: line.minority_allele,
            }) else {
                return Err(std::io::Error::other(format!(
                    "Cannot find variant position \"{position}\" with minority allele \"{min_allele}\" in file: '{sqm_file}'",
                    position = line.position,
                    min_allele = line.minority_allele as char,
                    sqm_file = args.sqm_file.display()
                )));
            };

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
/// Ensures the `Position` and `Minority_Allele` columns are present at the
/// locations expected by [`VariantsFileLine`], and rejects input that already
/// contains a `Phase` column.
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

/// A row from the `variants_file`, preserving the original line for output.
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

/// Final phase assignment indexed by variant position and minority allele.
struct VariantPhases {
    phase_by_tag: HashMap<VariantTag, usize>,
}

impl VariantPhases {
    /// Reads the `sqm_file`, parses variant tags and upper-triangle distances,
    /// validates the expected square dimensions, and assigns phase groups.
    fn from_sqm_file<I>(lines: I, tree_height: f64, number_of_variants: usize) -> std::io::Result<Self>
    where
        I: Iterator<Item = String>, {
        let mut variants_matrix = VariantsMatrix {
            number_of_variants,
            tag_index: HashMap::with_capacity(number_of_variants),
            variant_tags: Vec::with_capacity(number_of_variants),
            tree: UnionFindTree::new(number_of_variants),
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
                .parse_matrix_row(row_iter, row_ind, tree_height)
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

        let phase_group = variants_matrix.assign_phase_groups();
        Ok(variants_matrix.into_variant_phases(phase_group))
    }
}

/// Parsed square matrix data used while computing phase assignments.
struct VariantsMatrix {
    number_of_variants: usize,
    /// A map from the variant position and minority allele (tag) to the variant
    /// index.
    tag_index:          HashMap<VariantTag, usize>,
    /// The variant position and minority allele at each matrix row index.
    variant_tags:       Vec<VariantTag>,
    /// Connected components of variants whose pairwise distance is at or below
    /// the selected tree height.
    tree:               UnionFindTree,
}

impl VariantsMatrix {
    /// Parses the matrix values from one `sqm_file` row and merges variants
    /// whose pairwise distance is at or below `tree_height`.
    ///
    /// This only parses the upper triangular entries of the matrix, ignoring the
    /// diagonal and entries below it.
    fn parse_matrix_row(&mut self, row_iter: Split<'_, char>, row_ind: usize, tree_height: f64) -> std::io::Result<()> {
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

            if val <= tree_height {
                self.tree.union(row_ind, col_ind);
            }
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

    /// Records the relationship between an `.sqm` variant tag and its row
    /// index.
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

        match self.tag_index.entry(VariantTag { position, min_allele }) {
            Entry::Vacant(entry) => {
                entry.insert(row_ind);
                self.variant_tags.push(VariantTag { position, min_allele });
            }
            Entry::Occupied(_) => {
                return Err(std::io::Error::other(format!(
                    "Duplicate variant position and minority allele \"{tag}\" found."
                )));
            }
        };

        Ok(())
    }

    /// Labels phase clusters by size (descending), then by the lowest variant
    /// position, then by the lowest minority allele byte value.
    fn assign_phase_groups(&mut self) -> Vec<usize> {
        struct ClusterSummary {
            /// Size of the cluster.
            size:              usize,
            /// Smallest variant position within the cluster.
            min_position:      usize,
            /// Lowest minority allele byte value.
            lowest_min_allele: u8,
        }

        let mut clusters_by_root: HashMap<usize, ClusterSummary> = HashMap::with_capacity(self.number_of_variants);
        for (variant_index, tag) in self.variant_tags.iter().enumerate() {
            let root = self.tree.find(variant_index);

            clusters_by_root
                .entry(root)
                .and_modify(|cluster| {
                    cluster.size += 1;
                    cluster.min_position = cluster.min_position.min(tag.position);
                    cluster.lowest_min_allele = cluster.lowest_min_allele.min(tag.min_allele);
                })
                .or_insert_with(|| ClusterSummary {
                    size:              1,
                    min_position:      tag.position,
                    lowest_min_allele: tag.min_allele,
                });
        }

        let mut clusters_ordered = clusters_by_root.into_iter().collect::<Vec<_>>();
        clusters_ordered.sort_unstable_by_key(|(_, cluster)| {
            (
                std::cmp::Reverse(cluster.size),
                cluster.min_position,
                // breaks ties for same-size, same-position clusters
                cluster.lowest_min_allele,
            )
        });

        let mut phase_by_root = vec![0; self.tree.parents.len()];
        for (label, (root, _)) in clusters_ordered.into_iter().enumerate() {
            phase_by_root[root] = label + 1;
        }

        let mut phase_group = vec![0; self.number_of_variants];
        for (variant_index, phase) in phase_group.iter_mut().enumerate() {
            *phase = phase_by_root[self.tree.find(variant_index)];
        }

        phase_group
    }

    /// Converts row-indexed phase assignments into final tag-indexed lookup.
    fn into_variant_phases(self, phase_group: Vec<usize>) -> VariantPhases {
        VariantPhases {
            phase_by_tag: self.variant_tags.into_iter().zip(phase_group).collect(),
        }
    }
}

/// Position and minority allele identifier for a variant.
#[derive(Eq, PartialEq, Hash)]
struct VariantTag {
    /// Variant's position.
    position:   usize,
    /// Variant's minority allele.
    min_allele: u8,
}

/// A union-find tree structure for storing the phase clusters.
struct UnionFindTree {
    /// Parent links for each variant row index.
    parents: Vec<usize>,
    /// Component sizes for roots.
    sizes:   Vec<usize>,
}

impl UnionFindTree {
    fn new(number_of_variants: usize) -> Self {
        Self {
            parents: (0..number_of_variants).collect::<Vec<_>>(),
            sizes:   vec![1; number_of_variants],
        }
    }

    /// Returns the canonical root cluster for `cluster`, compressing the path
    /// on a second pass.
    fn find(&mut self, cluster: usize) -> usize {
        let mut root = cluster;
        while root != self.parents[root] {
            root = self.parents[root];
        }

        let mut node = cluster;
        while node != root {
            let parent = self.parents[node];
            self.parents[node] = root;
            node = parent;
        }

        root
    }

    /// Merges two variants into the same connected component.
    fn union(&mut self, cluster1: usize, cluster2: usize) {
        let mut parent_cluster1 = self.find(cluster1);
        let mut parent_cluster2 = self.find(cluster2);

        if parent_cluster1 == parent_cluster2 {
            return;
        }

        if self.sizes[parent_cluster1] < self.sizes[parent_cluster2] {
            std::mem::swap(&mut parent_cluster1, &mut parent_cluster2);
        }

        self.parents[parent_cluster2] = parent_cluster1;
        self.sizes[parent_cluster1] += self.sizes[parent_cluster2];
    }
}
