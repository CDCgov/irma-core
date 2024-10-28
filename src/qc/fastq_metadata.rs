use std::simd::prelude::*;

#[derive(Debug)]
pub(crate) struct FastQMetadata {
    pub(crate) passed_qc_count:               usize,
    pub(crate) passed_len_count:              usize,
    pub(crate) observed_q_max:                Option<f32>,
    pub(crate) observed_raw_reads:            Simd<usize, 2>,
    pub(crate) observed_max_read_len:         usize,
    pub(crate) observed_max_clipped_read_len: usize,
}

impl FastQMetadata {
    pub(crate) fn new() -> Self {
        FastQMetadata {
            passed_qc_count:               0,
            passed_len_count:              0,
            observed_q_max:                None,
            observed_raw_reads:            Simd::splat(0),
            observed_max_read_len:         0,
            observed_max_clipped_read_len: 0,
        }
    }
}
