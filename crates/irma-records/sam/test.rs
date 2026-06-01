use crate::sam::{
    ExpandedCigar,
    merge_pairs::{SamMergeablePairs, make_merged_qname},
};
use zoe::data::{cigar::Cigar, sam::SamData};

#[test]
fn merge_no_mismatch() {
    let s1 = SamData::new(
        "s1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "8M".try_into().unwrap(),
        b"AAAAAGGC".into(),
        b"FFFFFFFF".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "s2".to_string(),
        0,
        "ref".to_string(),
        11,
        30,
        "8M".try_into().unwrap(),
        b"GCGGTTTT".into(),
        b"FFFFFFFF".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTT";

    //   12345678901234567890
    //r  TTTTAAAAAGGCGGTTTT..
    //s1 ....AAAAAGGC........
    //s2 ......... GCGGTTTT..
    //m  ....AAAAAGGCGGTTTT..

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, b"AAAAAGGCGGTTTT".into());
    assert_eq!(m.qual, b"FFFFFFFFFFFFFF".try_into().unwrap());
    assert_eq!(Cigar::try_from(b"14M").unwrap(), m.cigar);
}

#[test]
fn merge_with_mismatch() {
    let s1 = SamData::new(
        "s1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "8M".try_into().unwrap(),
        b"AAAAAGGC".into(),
        b"FFFFFFFF".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "s2".to_string(),
        0,
        "ref".to_string(),
        10,
        30,
        "8M".try_into().unwrap(),
        b"GCGGCTTT".into(),
        b"FFFFFFFA".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTT";

    //   12345678901234567890
    //r  TTTTAAAAAGGCGGTTTT..
    //s1 ....AAAAAGGC........
    //s2 .........GCGGCTTT...
    //m  ....AAAAAGGCGCTTT...

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, b"AAAAAGGCGCTTT".into());
    assert_eq!(m.qual, b"FFFFFFFFFFFFA".try_into().unwrap());
    assert_eq!(Cigar::try_from(b"13M").unwrap(), m.cigar);
}

#[test]
fn merge_no_overlap() {
    let s1 = SamData::new(
        "s1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "8M".try_into().unwrap(),
        b"AAAAAGGC".into(),
        b"HHHHHHHH".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "s2".to_string(),
        0,
        "ref".to_string(),
        15,
        30,
        "8M".try_into().unwrap(),
        b"TTTTAGGA".into(),
        b"IIIIIIII".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTTAGGA";

    //   1234567890123456789012
    //r  TTTTAAAAAGGCGGTTTTAGGA
    //s1 ....AAAAAGGC..........
    //s2 ..............TTTTAGGA
    //m  ....AAAAAGGC..TTTTAGGA

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, b"AAAAAGGCTTTTAGGA".into());
    assert_eq!(m.qual, b"HHHHHHHHIIIIIIII".try_into().unwrap());
    assert_eq!(Cigar::try_from(b"8M2N8M").unwrap(), m.cigar);
}

#[test]
fn merge_with_clipping() {
    let s1 = SamData::new(
        "SRR42.1.1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "1H10M".try_into().unwrap(),
        b"AAAAAGGCGG".into(),
        b"FFFFFEEEEE".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "SRR422.1.2".to_string(),
        0,
        "ref".to_string(),
        10,
        30,
        "5H5M5S".try_into().unwrap(),
        b"GGGGGTTTTT".into(),
        b"EEEEE!!!!!".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTT";

    //   12345678901234567890
    //r  TTTTAAAAAGGCGGTTTT..
    //s1 ...hAAAAAGGCGG......
    //s2 ....hhhhhGGGGGsssss.
    //m  ...hAAAAAGGCGGhhhhh.

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, s1.seq);
    assert_eq!(m.qual, s1.qual);
    assert_eq!(Cigar::try_from(b"1H10M5H").unwrap(), m.cigar);
}

#[test]
fn merge_with_clipping2() {
    let s1 = SamData::new(
        "SRR42.1.1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "1H10M2H".try_into().unwrap(),
        b"AAAAAGGCGG".into(),
        b"FFFFFEEEEE".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "SRR422.1.2".to_string(),
        0,
        "ref".to_string(),
        10,
        30,
        "3H2S5M5S3H".try_into().unwrap(),
        b"TTGGGGGTTTTT".into(),
        b"EEEEEEE!!!!!".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTT";

    //   1234567890123456789012
    //r  TTTTAAAAAGGCGGTTTT....
    //s1 ...hAAAAAGGCGGhh......
    //s2 ....hhhssGGGGGssssshhh
    //m  ...hAAAAAGGCGGhhhhhhhh

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, s1.seq);
    assert_eq!(m.qual, s1.qual);
    assert_eq!(Cigar::try_from(b"1H10M8H").unwrap(), m.cigar);
}

#[test]
fn merge_with_clipping_past_left() {
    let s1 = SamData::new(
        "s1".to_string(),
        0,
        "ref".to_string(),
        5,
        30,
        "8H10M2H".try_into().unwrap(),
        b"AAAAAGGCGG".into(),
        b"FFFFFEEEEE".try_into().unwrap(),
    );

    let s2 = SamData::new(
        "s2".to_string(),
        0,
        "ref".to_string(),
        10,
        30,
        "3H2S5M5S3H".try_into().unwrap(),
        b"TTGGGGGTTTTT".into(),
        b"EEEEEEE!!!!!".try_into().unwrap(),
    );

    let reference = b"TTTTAAAAAGGCGGTTTT";

    //   6789123456789012345678
    //r  ....TTTTAAAAAGGCGGTTTT
    //s1 hhhhhhhhAAAAAGGCGGhh......
    //s2 ........hhhssGGGGGssssshhh
    //m  hhhhhhhhAAAAAGGCGGhhhhhhhh

    let (m, _) = s1.merge_pair_using_reference(&s2, reference, false);
    assert_eq!(m.seq, s1.seq);
    assert_eq!(m.qual, s1.qual);
    assert_eq!(Cigar::try_from(b"8H10M8H").unwrap(), m.cigar);
}

static QNAMES: [&str; 26] = [
    "SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
    "SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    "SRR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
    "SRR26182418.1.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    "A00350:691:HCKYLDSX3:2:2119:23863:2456/2",
    "A00350:691:HCKYLDSX3:2:2119:23863:2456/1",
    "M02989:9:000000000-L4PJL:1:2112:9890:15606 1:N:0:AACGCACGAG+GCCTCGGATA",
    "M02989:9:000000000-L4PJL:1:2112:9890:15606 2:N:0:AACGCACGAG+GCCTCGGATA",
    "NS500500:69:HKJFLAFX5:1:11204:14878:14643 1:N:0:TTCTCGTGCA+CTCTGTGTAT",
    "NS500500:69:HKJFLAFX5:1:11204:14878:14643 2:N:0:TTCTCGTGCA+CTCTGTGTAT",
    "A01000:249:HJFFWDRX2:1:2107:24605:18082 1:N:0:TAGGCATG+ATAGCCTT",
    "A01000:249:HJFFWDRX2:1:2107:24605:18082 2:N:0:TAGGCATG+ATAGCCTT",
    "M02989:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
    "M02989:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
    "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
    "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
    "SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
    "SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
    "SRR26182418.1.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
    "SRR26182418.1.2_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
    "SRR26182418.1 1:N:18:NULL",
    "SRR26182418.1.1 1:N:18:NULL",
    "ERR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
    "DRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    "ERR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
    "DRR26182418.2.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
];

#[test]
fn test_make_merged_qname() {
    let merged = [
        "SRR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "SRR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "SRR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "SRR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "A00350:691:HCKYLDSX3:2:2119:23863:2456/3",
        "A00350:691:HCKYLDSX3:2:2119:23863:2456/3",
        "M02989:9:000000000-L4PJL:1:2112:9890:15606 3:N:0:AACGCACGAG+GCCTCGGATA",
        "M02989:9:000000000-L4PJL:1:2112:9890:15606 3:N:0:AACGCACGAG+GCCTCGGATA",
        "NS500500:69:HKJFLAFX5:1:11204:14878:14643 3:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "NS500500:69:HKJFLAFX5:1:11204:14878:14643 3:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "A01000:249:HJFFWDRX2:1:2107:24605:18082 3:N:0:TAGGCATG+ATAGCCTT",
        "A01000:249:HJFFWDRX2:1:2107:24605:18082 3:N:0:TAGGCATG+ATAGCCTT",
        "M02989:9:000000000-L4PJL:1:2114:17393:19614_3:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989:9:000000000-L4PJL:1:2114:17393:19614_3:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_3:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_3:N:0:CTCTGCAGCG+GATGGATGTA",
        "SRR26182418.1.3_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "SRR26182418.1.3_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "SRR26182418.1.3_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "SRR26182418.1.3_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "SRR26182418.1.3 1:N:18:NULL",
        "SRR26182418.1.3 1:N:18:NULL",
        "ERR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "DRR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "ERR26182418.1.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "DRR26182418.2.3 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    ];

    for (i, o) in QNAMES.iter().enumerate() {
        assert_eq!(make_merged_qname(o), merged[i], "'{o}'");
    }
}

#[test]
fn test_expand() {
    let cigar = Cigar::from_slice_unchecked("4S10M2I2D3M4H4P");
    let expanded: ExpandedCigar = "SSSSMMMMMMMMMMIIDDMMMHHHHPPPP".into();
    assert_eq!(ExpandedCigar::from(cigar), expanded);
}

#[test]
fn test_condense_cigar() {
    let cigars: [&str; 5] = ["4S10M2I2D3M4H4P", "", "3M2D1M", "255M", "3M1D4I8X9=4M"];

    for c in cigars {
        let cigar = Cigar::from_slice_unchecked(c);
        assert_eq!(ExpandedCigar::from(cigar.clone()).condense_to_cigar(), cigar);
    }
}
