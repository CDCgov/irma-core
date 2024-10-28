use foldhash::fast::{FixedState, RandomState};
use std::{
    borrow::Borrow,
    collections::{
        hash_map::{Entry, IntoIter},
        HashMap,
    },
    env,
    hash::Hash,
};

const SEED_ENV_VAR: &str = "IRMA_SEED";

pub fn get_seed() -> Option<u64> {
    env::var(SEED_ENV_VAR).ok().map(|s| s.bytes().fold(0, |a, b| a ^ b) as u64)
}

pub enum SeedableFoldHashMap<K, V> {
    Seeded(HashMap<K, V, FixedState>),
    Random(HashMap<K, V, RandomState>),
}

impl<K, V> SeedableFoldHashMap<K, V> {
    // TODO: Inline???
    pub fn new(seed: Option<u64>) -> SeedableFoldHashMap<K, V> {
        match seed {
            Some(seed) => SeedableFoldHashMap::Seeded(HashMap::with_hasher(FixedState::with_seed(seed))),
            None => SeedableFoldHashMap::Random(HashMap::with_hasher(RandomState::default())),
        }
    }
}

impl<K, V> SeedableFoldHashMap<K, V>
where
    K: Eq + Hash,
{
    #[inline]
    pub fn get<Q>(&self, k: &Q) -> Option<&V>
    where
        K: Borrow<Q>,
        Q: Hash + Eq + ?Sized, {
        match &self {
            SeedableFoldHashMap::Seeded(map) => map.get(k),
            SeedableFoldHashMap::Random(map) => map.get(k),
        }
    }

    #[inline]
    pub fn insert(&mut self, k: K, v: V) -> Option<V> {
        match self {
            SeedableFoldHashMap::Seeded(map) => map.insert(k, v),
            SeedableFoldHashMap::Random(map) => map.insert(k, v),
        }
    }

    #[inline]
    pub fn entry(&mut self, key: K) -> Entry<'_, K, V> {
        match self {
            SeedableFoldHashMap::Seeded(map) => map.entry(key),
            SeedableFoldHashMap::Random(map) => map.entry(key),
        }
    }
}

impl<K, V> IntoIterator for SeedableFoldHashMap<K, V> {
    type Item = (K, V);
    type IntoIter = IntoIter<K, V>;

    #[inline]
    fn into_iter(self) -> IntoIter<K, V> {
        match self {
            SeedableFoldHashMap::Seeded(map) => map.into_iter(),
            SeedableFoldHashMap::Random(map) => map.into_iter(),
        }
    }
}
