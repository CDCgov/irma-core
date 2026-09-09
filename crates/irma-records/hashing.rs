use foldhash::{SharedSeed, fast::SeedableRandomState};
use std::env;

const SEED_ENV_VAR: &str = "IRMA_SEED";

/// Attempts to parse the environment variable `IRMA_SEED` into a `u64` to use
/// as a seed. If a seed exists but is unable to be parsed, this falls back to a
/// hashing algorithm.
pub fn get_seed() -> Option<u64> {
    env::var(SEED_ENV_VAR)
        .ok()
        .map(|value| value.parse::<u64>().unwrap_or_else(|_| seed_from_string(&value)))
}

fn seed_from_string(value: &str) -> u64 {
    let mut result = [0u8; 8];

    let (chunks, remainder) = value.as_bytes().as_chunks::<8>();
    // takes the seed in chunks of 8 bytes, and XORs the ith byte of all chunks
    // against eachother
    for chunk in chunks {
        for (position, byte) in chunk.iter().enumerate() {
            result[position] ^= byte;
        }
    }
    // handle the bytes in the remainder
    for (position, byte) in remainder.iter().enumerate() {
        result[position] ^= byte;
    }
    // then folds the bytes back into a u64 by shifting each byte by 8*i, then
    // ORing it against the accumulator
    result
        .into_iter()
        .enumerate()
        .fold(0u64, |seed, (position, byte)| seed | (u64::from(byte) << (8 * position)))
}

/// Creates a hasher based on the seed provided in the environment variable
/// `IRMA_SEED`, or falls back to creating a random one.
pub fn get_hasher() -> SeedableRandomState {
    match get_seed() {
        Some(seed) => SeedableRandomState::with_seed(seed, SharedSeed::global_fixed()),
        None => SeedableRandomState::random(),
    }
}

#[test]
fn test_seed_from_string() {
    // i worked this one out by hand
    let string = "ACGTACGTTGCATGCA";
    assert_eq!(seed_from_string(string), 1514339863296738325);
}
