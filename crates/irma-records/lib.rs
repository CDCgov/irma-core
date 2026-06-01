#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]
#![feature(portable_simd, try_trait_v2, int_format_into)]

pub mod fastq;
pub mod hashing;
pub mod io;
pub mod paired;
pub mod sam;
