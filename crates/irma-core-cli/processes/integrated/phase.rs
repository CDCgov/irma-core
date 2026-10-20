//! Phase distance calculation and group assignment.
//!
//! The central IRMA script supplies two TSV files exported from `call.pl`. The
//! Frequencies TSV has header columns `Position`, `Allele`, and `Frequency`.
//! Positions are positive and 1-based. The Counted Patterns TSV has header
//! columns `Pattern` and `Count`. Pattern characters follow sorted unique
//! positions in the frequency TSV.
//!
//! All four [legacy distance
//! measures]https://wonder.cdc.gov/controller/pdf/Updated_Phasing_Methods.pdf)
//! are computed once per pair, retaining the lower triangle in memory. Phase
//! groups are determined by identifying connected components where the EXPENRD
//! distance is below a threshold. Only explicitly selected matrices are
//! written. After validating all inputs and writing the matrices, the variants
//! table is overwritten in place with a trailing `Phase` column.

use clap::Args;
use irma_records::io::{Finish, InputOptions, OutputOptions};
use rayon::{
    ThreadPoolBuilder,
    iter::{IndexedParallelIterator, ParallelIterator},
    slice::ParallelSliceMut,
};
use std::{
    cmp::Reverse,
    collections::{BTreeMap, HashMap},
    fmt::{self, Display},
    io::{BufRead, Write},
    num::{NonZeroU64, NonZeroUsize},
    ops::Index,
    path::{Path, PathBuf},
};
use zoe::data::{
    err::ResultWithErrorContext,
    nucleotides::{CheckNucleotides, IsValidDNA::IupacWithGapsUc},
};

#[derive(Args, Debug)]
pub struct PhaseArgs {
    /// Prefix for the variants table and selected square matrices.
    pub ctype:               PathBuf,
    /// TSV with Position, Allele, Frequency columns.
    #[arg(long)]
    pub variant_frequencies: PathBuf,
    /// TSV with Pattern, Count columns.
    #[arg(long)]
    pub counted_patterns:    PathBuf,
    /// Maximum number of local worker threads available.
    #[arg(short, long, default_value = "1")]
    pub threads:             NonZeroUsize,
    /// Maximum EXPENRD distance for connected-component phase grouping.
    #[arg(short, long, value_parser = validate_threshold, default_value_t = 0.78)]
    pub distance_threshold:  f64,
    /// Comma-separated matrices: mutuald,jaccard,expenrd,njointp
    /// (case-insensitive). At least one valid name is required.
    #[arg(short, long, required = true, value_delimiter = ',', value_parser = parse_matrix_name)]
    pub matrices:            Vec<DistanceType>,
}

/// Validates the type and range of the `distance_threshold` argument.
fn validate_threshold(value: &str) -> Result<f64, String> {
    match value.parse::<f64>() {
        Ok(distance) if (0.0..=1.0).contains(&distance) => Ok(distance),
        _ => Err(String::from("Distance threshold must be a finite number between 0 and 1")),
    }
}

/// Parses a single matrix name, ignoring case.
fn parse_matrix_name(value: &str) -> Result<DistanceType, String> {
    DistanceType::ALL_TYPES
        .into_iter()
        .find(|matrix| matrix.name().eq_ignore_ascii_case(value))
        .ok_or_else(|| format!("Unknown matrix name {value}: expected mutuald, jaccard, expenrd, or njointp"))
}

pub fn phase_process(args: PhaseArgs) -> std::io::Result<()> {
    let variants_path = suffixed_path(&args.ctype, "-variants.txt");
    let table_reader = InputOptions::new_from_path(&variants_path).use_file().open()?;
    let table = VariantsTable::read(table_reader).with_path_context("Failed to parse variants table", &variants_path)?;

    match table.rows.len() {
        0 => return table.overwrite(&variants_path, std::iter::empty()),
        1 => return table.overwrite(&variants_path, std::iter::once(1)),
        _ => {}
    }

    let variants_reader = InputOptions::new_from_path(&args.variant_frequencies).use_file().open()?;
    let variants = read_variant_frequencies(variants_reader)
        .with_path_context("Failed to parse variant frequencies input", args.variant_frequencies)?;

    let site_count = variants.last().map_or(0, |variant| variant.site_index + 1);
    let patterns_reader = InputOptions::new_from_path(&args.counted_patterns).use_file().open()?;
    let patterns = read_counted_patterns(patterns_reader, &variants, site_count)
        .with_path_context("Failed to parse counted patterns input", args.counted_patterns)?;

    let num_threads = args.threads.get();
    let distances = if num_threads == 1 || row_offset(variants.len()) < num_threads {
        calculate_distances_sequential(&variants, &patterns, args.distance_threshold)?
    } else {
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .with_context("Failed to build Rayon thread pool")?;
        pool.install(|| calculate_distances_parallel(&variants, &patterns, args.distance_threshold))?
    };

    let phases = assign_phases(&variants, &distances.connected_pairs);

    write_matrices(
        &args.ctype,
        &variants,
        &distances.values,
        &DistanceType::select(&args.matrices),
    )?;

    let row_indices = table.match_variants(&variants)?;
    table.overwrite(&variants_path, row_indices.into_iter().map(|i| phases[i]))
}

/// A variant allele and its frequency, ordered by position and allele.
struct Variant {
    /// Positive, 1-based position of the variant
    position:   usize,
    /// Nucleotide allele.
    allele:     u8,
    /// Frequency of this allele in the frequency TSV.
    frequency:  f64,
    /// 0-based index of this variant's position among sorted unique positions.
    site_index: usize,
}

/// Reads and validates variant frequencies, assigning an index to each unique
/// position.
fn read_variant_frequencies(variant_frequencies: impl BufRead) -> std::io::Result<Vec<Variant>> {
    let mut lines = variant_frequencies.lines();
    validate_tsv_header(lines.next().transpose()?, "Position\tAllele\tFrequency")?;

    let mut variants = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        // TODO: replace with a `try` block if it's ever stabilized
        let result = (|| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let [pos, base, freq] = fields.as_slice() else {
                return Err(std::io::Error::other(
                    r#"Expected exactly three fields for "Position", "Allele", and "Frequency""#,
                ));
            };

            let position = parse_position(pos)?;
            let allele = parse_allele(base)?;
            let frequency = freq.parse::<f64>().map_err(std::io::Error::other)?;

            if !(0.0..=1.0).contains(&frequency) {
                return Err(std::io::Error::other("Variant frequency must be finite and between 0 and 1"));
            }

            let variant = Variant {
                position,
                allele,
                frequency,
                // placeholder, gets updated later in a second pass
                site_index: 0,
            };
            if variants.insert((position, allele), variant).is_some() {
                return Err(std::io::Error::other(format!(
                    "Duplicate variant {position}{allele}",
                    allele = char::from(allele)
                )));
            }

            Ok(())
        })();

        result.with_context(format!("Frequency TSV line {}", index + 2))?;
    }

    // Update each `site_index` field placeholder here while collecting
    let mut previous_position = None;
    let mut site_count = 0;
    let variants = variants
        .into_values()
        .map(|mut variant| {
            if previous_position != Some(variant.position) {
                previous_position = Some(variant.position);
                site_count += 1;
            }
            variant.site_index = site_count - 1;
            variant
        })
        .collect::<Vec<_>>();

    Ok(variants)
}

/// A pattern across unique variant positions and its observation count.
struct CountedPattern {
    /// Bases in sorted unique-position order.
    pattern: Box<[u8]>,
    /// Number of observations represented by this pattern.
    count:   u64,
}

/// Reads counted patterns, combining duplicate patterns and excluding patterns
/// containing only missing (`.`) or unknown (`N`) bases.
fn read_counted_patterns(
    counted_patterns: impl BufRead, variants: &[Variant], site_count: usize,
) -> std::io::Result<Vec<CountedPattern>> {
    let mut lines = counted_patterns.lines();
    validate_tsv_header(lines.next().transpose()?, "Pattern\tCount")?;

    let mut patterns = BTreeMap::<Box<[u8]>, u64>::new();
    let mut total_count = 0_u64;
    for (index, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        // TODO: another `try` block stand-in
        let result = (|| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let &[pattern, count] = fields.as_slice() else {
                return Err(std::io::Error::other(
                    r#"Expected exactly two fields for "Pattern" and "Count""#,
                ));
            };

            if variants.is_empty() || pattern.len() != site_count {
                return Err(std::io::Error::other(
                    "Patterns require variants and exactly one character per distinct variant position",
                ));
            }

            if !pattern.as_bytes().is_valid_dna(IupacWithGapsUc) {
                return Err(std::io::Error::other(
                    r#"The "Pattern" field contains an unsupported character"#,
                ));
            }

            let count = count
                .parse::<NonZeroU64>()
                .with_context("The pattern count must be a positive u64 integer")?
                .get();

            total_count = total_count
                .checked_add(count)
                .ok_or_else(|| std::io::Error::other("Total pattern count overflow"))?;

            // Patterns consisting only of missing/N sites are discarded
            // according to `call.pl`
            if !pattern.bytes().all(|base| b".N".contains(&base)) {
                let accumulated = patterns.entry(pattern.as_bytes().into()).or_default();
                *accumulated = accumulated
                    .checked_add(count)
                    .ok_or_else(|| std::io::Error::other("Pattern count overflow"))?;
            }

            Ok(())
        })();

        result.with_context(format!("Pattern TSV line {}", index + 2))?;
    }

    let counted_patterns = patterns
        .into_iter()
        .map(|(pattern, count)| CountedPattern { pattern, count })
        .collect();

    Ok(counted_patterns)
}

/// Requires the first TSV line to exactly match `expected`.
fn validate_tsv_header(header: Option<String>, expected: &str) -> std::io::Result<()> {
    let Some(header) = header else {
        return Err(std::io::Error::other(format!("Missing TSV header; expected {expected}")));
    };
    if header != expected {
        return Err(std::io::Error::other(format!(
            "Expected TSV header {expected}, found {header}"
        )));
    }
    Ok(())
}

/// Parses a single uppercase `A`, `C`, `G`, or `T` allele.
fn parse_allele(value: &str) -> std::io::Result<u8> {
    match *value.as_bytes() {
        [allele @ (b'A' | b'C' | b'G' | b'T')] => Ok(allele),
        _ => Err(std::io::Error::other("Variant allele must be one uppercase A, C, G, or T")),
    }
}

/// Parses a positive, 1-based variant position.
fn parse_position(value: &str) -> std::io::Result<usize> {
    let position = value
        .parse::<std::num::NonZeroUsize>()
        .with_context(format!("Expected a positive, 1-based position, found {value}"))?
        .get();

    Ok(position)
}

/// Appends `suffix` to the path's final component.
fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

// These column locations match the variants table in IRMA v1.4.0.
/// Index of the "Position" column in the `-variants.txt` table.
const POSITION_COLUMN: usize = 1;
/// Index of the "Minority_Allele" column in the `-variants.txt` table.
const MINORITY_ALLELE_COLUMN: usize = 4;
// TODO: Replace with max operator when it becomes const compatible
/// Minimum number of columns expected in the `-variants.txt` table
const MIN_COLUMNS: usize = if MINORITY_ALLELE_COLUMN > POSITION_COLUMN {
    MINORITY_ALLELE_COLUMN + 1
} else {
    POSITION_COLUMN + 1
};

/// The variants table header and parsed rows.
struct VariantsTable {
    /// Original TSV header.
    header: String,
    /// Data rows in their original order.
    rows:   Vec<VariantsFileLine>,
}

impl VariantsTable {
    /// Imports a `-variants.txt` table
    fn read(reader: impl BufRead) -> std::io::Result<Self> {
        let mut lines = reader.lines();

        let header = lines
            .next()
            .transpose()?
            .ok_or_else(|| std::io::Error::other("Variants table is empty"))?;

        Self::validate_variants_header(&header)?;

        let columns = header.split('\t').count();
        let mut rows = Vec::new();
        for (index, line) in lines.enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }

            if line.split('\t').count() != columns {
                return Err(std::io::Error::other(format!(
                    "Variants table line {} has incorrect column count",
                    index + 2
                )));
            }

            rows.push(VariantsFileLine::try_from(line).with_context(format!("Variants table line {}", index + 2))?);
        }
        Ok(Self { header, rows })
    }

    /// Validates the header in the `<ctype>-variants.txt` file.
    ///
    /// Ensures the `Position` and `Minority_Allele` columns are present at the
    /// locations expected by [`VariantsFileLine`], and rejects input that
    /// already contains a `Phase` column.
    fn validate_variants_header(header: &str) -> std::io::Result<()> {
        let fields = header.split('\t').collect::<Vec<_>>();
        if fields.len() < MIN_COLUMNS
            || fields[POSITION_COLUMN] != "Position"
            || fields[MINORITY_ALLELE_COLUMN] != "Minority_Allele"
        {
            return Err(std::io::Error::other(
                "Variants table requires Position in column 2 and Minority_Allele in column 5",
            ));
        }
        if fields.contains(&"Phase") {
            return Err(std::io::Error::other("Input variants table already has a Phase column"));
        }
        Ok(())
    }

    /// Matches each variants-table row to an entry in `variants` by position
    /// and allele.
    ///
    /// Returns indices into `variants` in table-row order, so per-variant
    /// results can be written back in the same order as the table.
    ///
    /// ## Errors
    ///
    /// Returns an error if a table row has no matching entry, if a table
    /// variant appears more than once, or if any entry in `variants` is not
    /// represented in the table.
    fn match_variants(&self, variants: &[Variant]) -> std::io::Result<Vec<usize>> {
        let mut by_tag = variants
            .iter()
            .enumerate()
            .map(|(i, v)| ((v.position, v.allele), i))
            .collect::<HashMap<_, _>>();
        let mut indices = Vec::with_capacity(self.rows.len());

        for row in &self.rows {
            let index = by_tag.remove(&(row.position, row.minority_allele)).ok_or_else(|| {
                std::io::Error::other(format!(
                    "Unexpected or duplicate table variant {}{}",
                    row.position,
                    char::from(row.minority_allele)
                ))
            })?;
            indices.push(index);
        }

        if !by_tag.is_empty() {
            return Err(std::io::Error::other(
                "Frequency TSV contains variants missing from the variants table",
            ));
        }

        Ok(indices)
    }

    /// Overwrites the table with a trailing `Phase` column in row order.
    fn overwrite(&self, path: &Path, phases: impl ExactSizeIterator<Item = usize>) -> std::io::Result<()> {
        debug_assert_eq!(
            phases.len(),
            self.rows.len(),
            "Phase assignment count does not match variants table row count",
        );

        let mut writer = OutputOptions::new_from_path(path).use_file().open()?;
        writeln!(writer, "{}\tPhase", self.header)?;
        for (row, phase) in self.rows.iter().zip(phases) {
            writeln!(writer, "{}\t{phase}", row.line)?;
        }
        writer.finish()
    }
}

/// A variants-table row, preserving its original text and parsed match key.
struct VariantsFileLine {
    /// Original tab-separated row text.
    line:            String,
    /// Positive, 1-based position used to match the frequency TSV.
    position:        usize,
    /// Minority allele used to match the frequency TSV.
    minority_allele: u8,
}

impl TryFrom<String> for VariantsFileLine {
    type Error = std::io::Error;

    fn try_from(line: String) -> std::io::Result<Self> {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < MIN_COLUMNS {
            return Err(std::io::Error::other(format!("Expected at least {MIN_COLUMNS} columns")));
        }
        let position = parse_position(fields[POSITION_COLUMN])?;
        let minority_allele = parse_allele(fields[MINORITY_ALLELE_COLUMN])?;
        Ok(Self {
            line,
            position,
            minority_allele,
        })
    }
}

/// Pairwise distances in packed lower-triangle order and threshold-qualified
/// pairs for phase assignment.
struct DistanceMatrices {
    /// Four distances per variant pair, excluding the diagonal.
    values:          Vec<CalculatedDistances>,
    /// Variant-index pairs whose EXPENRD distance is at most the threshold.
    connected_pairs: Vec<(usize, usize)>,
}

/// Calculates the packed pairwise distance matrix and qualifying pairs
/// serially.
fn calculate_distances_sequential(
    variants: &[Variant], patterns: &[CountedPattern], distance_threshold: f64,
) -> std::io::Result<DistanceMatrices> {
    let mut distances = vec![CalculatedDistances::SAME; row_offset(variants.len())];
    let mut connected_pairs = Vec::new();

    for i in 1..variants.len() {
        let row_start = row_offset(i);
        for j in 0..i {
            let value = distance(&variants[i], &variants[j], patterns);
            if value.expenrd <= distance_threshold {
                connected_pairs.push((i, j));
            }
            distances[row_start + j] = value;
        }
    }

    Ok(DistanceMatrices {
        values: distances,
        connected_pairs,
    })
}

/// Computes all pairwise [`Distances`] between `variants`.
///
/// Each off-diagonal pair is computed once with [`distance`] and stored in a
/// flattened lower-triangular matrix of per-pair structs, storing row `i`
/// followed by row `i + 1`, with `i` entries for pairs `(i, j)` where `j < i`.
/// Selected square matrix outputs are derived from these shared values by
/// [`write_matrices`].
///
/// Equal-pair chunks are processed in parallel. Within each chunk, row segments
/// share a fixed first variant and iterate over the corresponding column slice.
///
/// Pairs whose EXPENRD distance is at most `distance_threshold` are collected
/// during calculation, in packed-matrix order, for phase assignment.
fn calculate_distances_parallel(
    variants: &[Variant], patterns: &[CountedPattern], distance_threshold: f64,
) -> std::io::Result<DistanceMatrices> {
    let pairs = row_offset(variants.len());
    let mut distances = vec![CalculatedDistances::SAME; pairs];

    // This number was tuned empirically by local benchmarks based on expected
    // typical-sized inputs. If large numbers of variants are expected (>2,000),
    // `32` may improve speed by ~1%.
    const CHUNKS_PER_WORKER: usize = 8;
    let target_chunks = rayon::current_num_threads().saturating_mul(CHUNKS_PER_WORKER);
    let chunk_size = distances.len().div_ceil(target_chunks).max(1);

    let pairs_by_chunk = distances
        .par_chunks_mut(chunk_size)
        .enumerate()
        .map(|(chunk_index, chunk)| {
            let start = chunk_index * chunk_size;

            // Convert the packed index to strict lower-triangle coordinates.
            let mut i = (1 + 8 * start).isqrt().div_ceil(2);
            let mut j = start - row_offset(i);

            let mut connected_pairs = Vec::new();
            let mut remaining = chunk;

            while !remaining.is_empty() {
                let segment_len = remaining.len().min(i - j);
                let (segment, rest) = remaining.split_at_mut(segment_len);
                let first = &variants[i];
                for (column, (value, second)) in segment.iter_mut().zip(&variants[j..j + segment_len]).enumerate() {
                    *value = distance(first, second, patterns);
                    if value.expenrd <= distance_threshold {
                        connected_pairs.push((i, j + column));
                    }
                }

                remaining = rest;
                i += 1;
                j = 0;
            }
            connected_pairs
        })
        .collect::<Vec<_>>();

    Ok(DistanceMatrices {
        values:          distances,
        connected_pairs: pairs_by_chunk.into_iter().flatten().collect(),
    })
}

/// Starting offset of a row in a packed strict lower triangle.
fn row_offset(row: usize) -> usize {
    let previous = row.saturating_sub(1);
    if row.is_multiple_of(2) {
        (row / 2) * previous
    } else {
        row * (previous / 2)
    }
}

/// The four distance measures calculated for a variant pair.
#[derive(Clone, Debug)]
struct CalculatedDistances {
    mutuald: f64,
    jaccard: f64,
    expenrd: f64,
    njointp: f64,
}

impl CalculatedDistances {
    const NONE: Self = Self {
        mutuald: 1.0,
        jaccard: 1.0,
        expenrd: 1.0,
        njointp: 1.0,
    };

    const SAME: Self = Self {
        mutuald: 0.0,
        jaccard: 0.0,
        expenrd: 0.0,
        njointp: 0.0,
    };
}

impl Index<DistanceType> for CalculatedDistances {
    type Output = f64;

    fn index(&self, index: DistanceType) -> &Self::Output {
        match index {
            DistanceType::Mutuald => &self.mutuald,
            DistanceType::Jaccard => &self.jaccard,
            DistanceType::Expenrd => &self.expenrd,
            DistanceType::Njointp => &self.njointp,
        }
    }
}

/// A distance measure and its matrix output name.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum DistanceType {
    Mutuald,
    Jaccard,
    Expenrd,
    Njointp,
}

impl DistanceType {
    const ALL_TYPES: [Self; 4] = [Self::Mutuald, Self::Jaccard, Self::Expenrd, Self::Njointp];

    const fn name(self) -> &'static str {
        match self {
            Self::Mutuald => "MUTUALD",
            Self::Jaccard => "JACCARD",
            Self::Expenrd => "EXPENRD",
            Self::Njointp => "NJOINTP",
        }
    }

    /// Deduplicates the `parsed_matrices`.
    fn select(parsed_matrices: &[Self]) -> Vec<Self> {
        Self::ALL_TYPES
            .into_iter()
            .filter(|matrix_type| parsed_matrices.contains(matrix_type))
            .collect()
    }
}

impl Display for DistanceType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Computes the four legacy distance measures between two variants from the
/// counted `patterns`.
///
/// The returned [`Distances`] holds, in [`Distance::VARIANTS`] order:
///
/// 1. MUTUALD
/// 2. JACCARD
/// 3. EXPENRD
/// 4. NJOINTP
///
/// A pattern covers the pair unless either site is missing (`.`). Ambiguous
/// bases, `N`, and gaps (`-`) count toward coverage but never match an allele.
///
/// The following special cases apply:
///
/// - Two variants at the same site return [`Distances::SAME`] if they share an
///   allele, and [`Distances::NONE`] otherwise.
/// - If either variant has a frequency of zero, or there is no coverage, no
///   joint observation, or no observation of either allele, this returns
///   [`Distances::NONE`].
///
/// All measures are symmetric, so the order of `first` and `second` does not
/// matter.
fn distance(first: &Variant, second: &Variant, patterns: &[CountedPattern]) -> CalculatedDistances {
    if first.site_index == second.site_index {
        return if first.allele == second.allele {
            CalculatedDistances::SAME
        } else {
            CalculatedDistances::NONE
        };
    }
    if first.frequency == 0.0 || second.frequency == 0.0 {
        return CalculatedDistances::NONE;
    }

    let (mut total, mut joint, mut first_count, mut second_count) = (0.0, 0.0, 0.0, 0.0);

    for counted_pattern in patterns {
        let count = counted_pattern.count as f64;
        let (a, b) = (
            counted_pattern.pattern[first.site_index],
            counted_pattern.pattern[second.site_index],
        );
        if a == b'.' || b == b'.' {
            continue;
        }
        total += count;
        if a == first.allele {
            first_count += count;
        }
        if b == second.allele {
            second_count += count;
        }
        if a == first.allele && b == second.allele {
            joint += count;
        }
    }
    if total == 0.0 || joint == 0.0 || first_count == 0.0 || second_count == 0.0 {
        return CalculatedDistances::NONE;
    }

    let joint = joint / total;
    let first_local = first_count / total;
    let second_local = second_count / total;
    let minimum = first_local.min(first.frequency).min(second_local.min(second.frequency));
    let first_max = first_local.max(first.frequency);
    let second_max = second_local.max(second.frequency);

    let mutuald = 1.0 - joint * joint / (first_max * second_max);
    let jaccard = 1.0 - joint / (first_max + second_max - joint);
    let expenrd = if total <= 20.0 || mutuald == 0.0 || jaccard == 0.0 {
        1.0 - (joint * minimum) / (first_max * second_max)
    } else {
        jaccard
    };
    let njointp = 1.0 - 2.0 * joint;

    CalculatedDistances {
        mutuald,
        jaccard,
        expenrd,
        njointp,
    }
}

/// Assigns each variant a phase number from the connected components of the
/// graph defined by the threshold-qualified `connected_pairs`.
///
/// Each pair contains two indices into `variants`. Linking is transitive, so
/// two variants share a phase whenever a chain of qualifying pairs connects
/// them, even if their own distance is above the threshold.
///
/// Phases are numbered from 1, and components are ordered by:
///
/// 1. Size, largest first
/// 2. Smallest position in the component
/// 3. Smallest allele in the component
///
/// The returned vector has one phase per variant, in the same order as
/// `variants`.
fn assign_phases(variants: &[Variant], connected_pairs: &[(usize, usize)]) -> Vec<usize> {
    let mut components = UnionFindTree::new(variants.len());
    for &(i, j) in connected_pairs {
        components.union(i, j);
    }

    let roots = (0..variants.len()).map(|i| components.find(i)).collect::<Vec<_>>();
    let mut summaries = HashMap::new();
    for (variant, &root) in variants.iter().zip(&roots) {
        let summary = summaries.entry(root).or_insert((0_usize, usize::MAX, u8::MAX));
        summary.0 += 1;
        summary.1 = summary.1.min(variant.position);
        summary.2 = summary.2.min(variant.allele);
    }
    let mut ordered = summaries.into_iter().collect::<Vec<_>>();
    ordered.sort_unstable_by_key(|&(root, (size, position, allele))| (Reverse(size), position, allele, root));
    let mut phase_by_root = vec![0; variants.len()];
    for (index, (root, _)) in ordered.into_iter().enumerate() {
        phase_by_root[root] = index + 1;
    }
    roots.into_iter().map(|root| phase_by_root[root]).collect()
}

/// Union-find tree structure for grouping variants connected by qualifying
/// pairs.
struct UnionFindTree {
    /// Parent variant index for each node. Roots point to themselves.
    parents: Vec<usize>,
    /// Component size stored at each root.
    sizes:   Vec<usize>,
}

impl UnionFindTree {
    fn new(n: usize) -> Self {
        Self {
            parents: (0..n).collect(),
            sizes:   vec![1; n],
        }
    }

    /// Merges the components containing `first` and `second`.
    fn union(&mut self, first: usize, second: usize) {
        let (mut a, mut b) = (self.find(first), self.find(second));
        if a != b {
            if self.sizes[a] < self.sizes[b] {
                std::mem::swap(&mut a, &mut b);
            }
            self.parents[b] = a;
            self.sizes[a] += self.sizes[b];
        }
    }

    /// Finds a component root while compressing the traversed path.
    fn find(&mut self, mut index: usize) -> usize {
        while self.parents[index] != index {
            self.parents[index] = self.parents[self.parents[index]];
            index = self.parents[index];
        }
        index
    }
}

/// Writes selected packed distance measures as symmetric square matrices.
///
/// Lower-triangle entries advance sequentially. Upper-triangle entries skip
/// between packed rows with an increasing stride.
fn write_matrices(
    prefix: &Path, variants: &[Variant], distances: &[CalculatedDistances], matrix_types: &[DistanceType],
) -> std::io::Result<()> {
    if row_offset(variants.len()) != distances.len() {
        return Err(std::io::Error::other("Matrix size must match the number of variants"));
    }

    let mut outputs = matrix_types
        .iter()
        .map(|&matrix_type| {
            let path = suffixed_path(prefix, &format!("-{matrix_type}.sqm"));
            let writer = OutputOptions::new_from_path(&path).use_file().open()?;
            Ok((matrix_type, writer))
        })
        .collect::<std::io::Result<Vec<_>>>()?;

    let mut idx = 0;
    for (row, variant) in variants.iter().enumerate() {
        // The upper triangle has the values and linear indexes of the lower
        // triangle reflected about the diagonal. In the lower triangle, each
        // row has one more entry than the previous. `idx + row` is the start of
        // packed reflected row `row + 1` which corresponds to the column of the
        // first upper triangle entry.
        let mut upper_triangular_idx = idx + row;
        for (_, writer) in &mut outputs {
            write!(
                writer,
                "{position}{allele}",
                position = variant.position,
                allele = char::from(variant.allele)
            )?;
        }

        for column in 0..variants.len() {
            let pair_distances = if column < row {
                let pair_distances = &distances[idx];
                idx += 1;
                pair_distances
            } else if column == row {
                &CalculatedDistances::SAME
            } else {
                // For each upper-triangle column, skip `column - 1` entries to
                // reach the position in the next packed reflected row. The `-1`
                // is there because the reflected row immediately before column
                // has `column - 1` entries.
                upper_triangular_idx += column - 1;
                &distances[upper_triangular_idx]
            };
            for (matrix_type, writer) in &mut outputs {
                write!(writer, "\t{}", pair_distances[*matrix_type])?;
            }
        }

        for (_, writer) in &mut outputs {
            writeln!(writer)?;
        }
    }

    for (_, output) in outputs {
        output.finish()?;
    }

    Ok(())
}
