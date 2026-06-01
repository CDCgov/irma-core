use std::ops::{Add, AddAssign};

/// [`PairedMergeStats`] holds statistics related to read pair merging operations.
#[derive(Copy, Clone, Debug, Default)]
pub struct PairedMergeStats {
    /// Total number of overlapping or paired bases
    pub observations:    u64,
    /// Paired bases that agree with each other but disagree with consensus
    pub true_variations: u64,
    /// Paired bases that disagree with each other
    pub variant_errors:  u64,
    /// Paired bases where one is a deletion and one is not
    pub deletion_errors: u64,
    /// Total number of paired insertions, in agreement or otherwise
    pub insert_obs:      u64,
    /// Total number of mismatching paired insertions, including disagreement in
    /// insertion presence
    pub insert_errors:   u64,
}

impl Add for PairedMergeStats {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self {
            observations:    self.observations + other.observations,
            true_variations: self.true_variations + other.true_variations,
            variant_errors:  self.variant_errors + other.variant_errors,
            deletion_errors: self.deletion_errors + other.deletion_errors,
            insert_obs:      self.insert_obs + other.insert_obs,
            insert_errors:   self.insert_errors + other.insert_errors,
        }
    }
}

impl AddAssign for PairedMergeStats {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            observations:    self.observations + other.observations,
            true_variations: self.true_variations + other.true_variations,
            variant_errors:  self.variant_errors + other.variant_errors,
            deletion_errors: self.deletion_errors + other.deletion_errors,
            insert_obs:      self.insert_obs + other.insert_obs,
            insert_errors:   self.insert_errors + other.insert_errors,
        }
    }
}
