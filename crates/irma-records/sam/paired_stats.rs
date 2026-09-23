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
