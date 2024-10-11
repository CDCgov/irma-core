#[derive(Debug)]
pub(crate) struct FastQMetadata {
    pub(crate) passed_qc_count:              usize,
    pub(crate) passed_len_count:             usize,
    pub(crate) dataset_q_max:                Option<f32>,
    pub(crate) dataset_max_read_len:         Option<usize>,
    pub(crate) dataset_max_clipped_read_len: Option<usize>,
}


impl FastQMetadata {
    pub(crate) fn new() -> Self {
        FastQMetadata {
            passed_qc_count:              0,
            passed_len_count:             0,
            dataset_q_max:                None,
            dataset_max_read_len:         None,
            dataset_max_clipped_read_len: None,
        }
    }

    pub(crate) fn merge(mut self, other: Self) -> Self {
        self.passed_qc_count += other.passed_qc_count;
        self.passed_len_count += other.passed_len_count;
        self.dataset_q_max = match (self.dataset_q_max, other.dataset_q_max) {
            (Some(a), Some(b)) => Some(if a > b { a } else { b }),
            (Some(a), None) | (None, Some(a)) => Some(a),
            _ => None,
        };
        self.dataset_max_read_len = match (self.dataset_max_read_len, other.dataset_max_read_len) {
            (Some(a), Some(b)) => Some(a.max(b)),
            (Some(a), None) | (None, Some(a)) => Some(a),
            _ => None,
        };
        self.dataset_max_clipped_read_len = match (self.dataset_max_clipped_read_len, other.dataset_max_clipped_read_len) {
            (Some(a), Some(b)) => Some(a.max(b)),
            (Some(a), None) | (None, Some(a)) => Some(a),
            _ => None,
        };
        self
    }
}
