use foldhash::{SharedSeed, fast::SeedableRandomState};
use std::env;

pub(crate) mod paired_reads;
pub(crate) mod trimming;

const SEED_ENV_VAR: &str = "IRMA_SEED";

fn get_seed() -> Option<u64> {
    env::var(SEED_ENV_VAR).ok().map(|s| s.bytes().fold(0, |a, b| a ^ b) as u64)
}

pub fn get_hasher() -> SeedableRandomState {
    match get_seed() {
        Some(seed) => SeedableRandomState::with_seed(seed, SharedSeed::global_fixed()),
        None => SeedableRandomState::random(),
    }
}
