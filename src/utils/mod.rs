use foldhash::{SharedSeed, fast::SeedableRandomState};
use std::env;

const SEED_ENV_VAR: &str = "IRMA_SEED";

fn get_seed() -> Option<u64> {
    env::var(SEED_ENV_VAR).ok().map(|s| s.bytes().fold(0, |a, b| a ^ b) as u64)
}

pub fn get_hasher() -> SeedableRandomState {
    match get_seed() {
        Some(seed) => SeedableRandomState::with_seed(seed, SharedSeed::global_random()),
        None => SeedableRandomState::random(),
    }
}

/// Takes a fastq header and returns the molecular ID and side (for paired
/// reads)
pub fn get_molecular_id_side(s: &str, default_side: char) -> Option<(&str, char)> {
    let (the_id, the_side) = if s.contains(' ') {
        let mut pieces = s.split(' ');
        let id = pieces.next().unwrap_or_default();

        if !(id.starts_with("SRR") || id.starts_with("DRR") || id.starts_with("ERR")) || !id.contains('.') {
            // Illumina format
            (
                id,
                pieces
                    .next()
                    .unwrap_or_default()
                    .split(':')
                    .next()
                    .unwrap_or_default()
                    .chars()
                    .next(),
            )
        } else if let Some(index) = id.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA format, read side included
            let (new_id, side) = id.split_at(index);
            (new_id, side.chars().nth(1))
        } else {
            // SRA format, no read side
            (id, Some(default_side))
        }
    } else if let Some(index) = s.find('/') {
        // Legacy Illumina
        let (new_id, side) = s.split_at(index);
        (new_id, side.chars().nth(1))
    } else if (s.starts_with("SRR") || s.starts_with("DRR") || s.starts_with("ERR")) && s.contains('.') {
        let mut pieces = s.split('_');
        let id = pieces.next().unwrap_or_default();

        if let Some(index) = id.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA with read side
            let (new_id, side) = id.split_at(index);
            (new_id, side.chars().nth(1))
        } else {
            // SRA, no read side
            (id, Some(default_side))
        }
    } else {
        // IRMA Illumina legacy output
        let mut indices = s.match_indices(':');
        let (left, right) = (indices.nth(5), indices.next());
        if let (Some((start, _)), Some((stop, _))) = (left, right)
            && let Some(us) = s[start..stop].find('_')
        {
            let underscore_index = start + us;
            (&s[..underscore_index], s[..stop].chars().next_back())
        } else {
            return None;
        }
    };

    if let (id, Some(side @ '0'..='3')) = (the_id, the_side) {
        Some((id, side))
    } else {
        Some((the_id, default_side))
    }
}
