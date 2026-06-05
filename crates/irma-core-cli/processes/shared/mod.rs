pub mod trimming;

use zoe::search::ByteSubstringMut;

/// Replaces tabs with spaces in a String.
///
/// The legacy xflate format (XFL) uses tabs as a delimiter and these cannot be
/// present in the FASTQ header.
pub(crate) fn replace_tabs_with_spaces(header: &mut String) {
    // Safety: tab to space replacements are both ASCII and therefore do not
    // disrupt UTF-8 boundaries
    unsafe {
        header.as_mut_vec().replace_all_bytes(b'\t', b' ');
    }
}
