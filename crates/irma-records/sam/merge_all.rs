use crate::{
    hashing::get_hasher,
    paired::{ReadSide, get_molecular_id_side},
    sam::{PairedMergeStats, SamMergeablePairs},
};
use foldhash::fast::SeedableRandomState;
use std::collections::{HashMap, hash_map::IntoValues};
use zoe::data::{sam::SamData, views::Len};

/// A record awaiting its mate, boxed to keep the hash table compact.
type UnpairedAln = (ReadSide, Box<SamData>);
type UnpairedMap = HashMap<String, UnpairedAln, SeedableRandomState>;

/// Adds `read` to the pool of available mates, or removes its mate and returns
/// the pair. A record without a read side fills whichever side its mate lacks.
/// If neither has a side, the lesser `qname` is R1.
///
/// ## Validity
///
/// - Each molecule has at most two records, with at most one explicit R1 and
///   one explicit R2. Otherwise, a repeated side is dropped while its first
///   record awaits a mate, and redundant records after a completed pair can
///   lead to new, duplicate pairs or singletons.
fn pair_record(available_mates: &mut UnpairedMap, read: SamData) -> Option<(SamData, SamData)> {
    let (id, side) = match get_molecular_id_side(&read.qname, '0') {
        Some((id, '1')) => (id, ReadSide::R1),
        Some((id, '2')) => (id, ReadSide::R2),
        Some((id, _)) => (id, ReadSide::None),
        None => (read.qname.as_str(), ReadSide::None),
    };

    let Some((mate_side, mate)) = available_mates.remove(id) else {
        available_mates.insert(id.to_string(), (side, Box::new(read)));
        return None;
    };

    // Should not happen, but we preserve legacy behavior.
    if mate_side == side && side != ReadSide::None {
        available_mates.insert(id.to_string(), (side, mate));
        return None;
    }

    Some(
        // If both reads in the pair do not have a side, the merged `qname` is
        // deterministically based on lexicographic order of the qnames.
        if ReadSide::R1 == side || mate_side == ReadSide::R2 || (side == mate_side && read.qname < mate.qname) {
            (read, *mate)
        } else {
            (*mate, read)
        },
    )
}

/// An iterator that pairs up records by molecule. Each pair is yielded as soon
/// as its second mate is read, and records left without a mate are yielded once
/// the input is exhausted.
struct Pairer<I> {
    iter:  I,
    cache: PairerState,
}

/// Models state machine for streaming read pairer.
enum PairerState {
    /// Reading from the iterator and making pairs.
    Reading(UnpairedMap),
    /// For draining unpaired reads.
    Draining(IntoValues<String, UnpairedAln>),
}

impl<I> Iterator for Pairer<I>
where
    I: Iterator<Item = SamData>,
{
    type Item = (SamData, Option<SamData>);

    fn next(&mut self) -> Option<Self::Item> {
        use PairerState::*;

        match &mut self.cache {
            Reading(available_mates) => {
                for read in self.iter.by_ref() {
                    if let Some((r1, r2)) = pair_record(available_mates, read) {
                        return Some((r1, Some(r2)));
                    }
                }

                self.cache = Draining(std::mem::take(available_mates).into_values());
                self.next()
            }
            Draining(remaining) => {
                let (_, unpaired) = remaining.next()?;
                Some((*unpaired, None))
            }
        }
    }
}

/// Merges `r1` and its mate `r2`, if any. Returns the merged record, or the
/// records unchanged if there is no mate or either lacks quality scores (for
/// which IRMA does not define merging).
fn merge_mates(
    r1: SamData, r2: Option<SamData>, reference: &[u8], bowtie_format: bool, stats: &mut PairedMergeStats,
) -> [Option<SamData>; 2] {
    let Some(r2) = r2 else {
        return [Some(r1), None];
    };

    if r1.qual.is_empty() || r2.qual.is_empty() {
        [Some(r1), Some(r2)]
    } else {
        let merged = r1.merge_pair_using_reference(&r2, reference, bowtie_format, stats);
        [Some(merged), None]
    }
}

/// Merges read pairs that are already matched up by molecule, using the
/// `reference` sequence and accumulating observations in `stats`.
///
/// Each item is an R1 record and its R2 mate, if any. The R1 record provides
/// the merged qname. A pair is yielded unmerged if either record lacks quality
/// scores, since IRMA does not define merging for that case.
///
/// ## Validity
///
/// - Record sequences must be uppercase.
/// - Records belong to the supplied reference.
//  TODO: possible entry point for `i-assemble`
pub fn merge_pairs(
    pairs: impl IntoIterator<Item = (SamData, Option<SamData>)>, reference: &[u8], bowtie_format: bool,
    stats: &mut PairedMergeStats,
) -> impl Iterator<Item = SamData> {
    let reference = reference.to_ascii_uppercase();
    pairs
        .into_iter()
        .flat_map(move |(r1, r2)| merge_mates(r1, r2, &reference, bowtie_format, stats))
        .flatten()
}

/// Pairs up R1 and R2 `records` by molecule and merges each pair as in
/// [`merge_pairs`], yielding it as soon as its second mate is read. Records
/// without a mate are yielded once `records` is exhausted.
///
/// The read side is taken from the qname when it is explicit. A record without
/// a read side fills whichever side its mate lacks. If neither mate has a side,
/// the one with the lesser qname is treated as R1.
///
/// Merged content is independent of input order; emission order follows pair
/// completion, followed by unmatched records in hash table order.
///
/// ## Validity
///
/// - Record sequences must be uppercase.
/// - Records belong to the supplied reference.
/// - Each molecule has at most two records, with at most one explicit R1 and
///   one explicit R2. Either or both records may lack a read side.
pub fn pair_and_merge(
    records: impl IntoIterator<Item = SamData>, reference: &[u8], bowtie_format: bool, stats: &mut PairedMergeStats,
) -> impl Iterator<Item = SamData> {
    let reference = reference.to_ascii_uppercase();
    let pairs = Pairer {
        iter:  records.into_iter(),
        cache: PairerState::Reading(HashMap::with_hasher(get_hasher())),
    };

    pairs
        .flat_map(move |(r1, r2)| merge_mates(r1, r2, &reference, bowtie_format, stats))
        .flatten()
}
