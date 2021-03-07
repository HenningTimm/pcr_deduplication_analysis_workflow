//! Documentation for Rust Bio Tools
use clap::{load_yaml, value_t};
use log::LevelFilter;

use clap::App;
use fern;
use std::error::Error;

pub mod common;
pub mod fastq;

fn main() -> Result<(), Box<dyn Error>> {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
        .version(env!("CARGO_PKG_VERSION"))
        .get_matches();

    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    match matches.subcommand() {
        ("call-consensus-reads", Some(matches)) => match matches.subcommand() {
            ("fastq", Some(matches)) => {
                fastq::call_consensus_reads::call_consensus_reads_from_paths(
                    matches.value_of("fq1").unwrap(),
                    matches.value_of("fq2").unwrap(),
                    matches.value_of("consensus-fq1").unwrap(),
                    matches.value_of("consensus-fq2").unwrap(),
                    matches.value_of("consensus-fq3"),
                    value_t!(matches, "umi-len", usize).unwrap(),
                    value_t!(matches, "max-seq-dist", usize).unwrap(),
                    value_t!(matches, "max-umi-dist", usize).unwrap(),
                    matches.is_present("umi-on-reverse"),
                    matches.is_present("verbose-read-names"),
                    if matches.is_present("insert-size") {
                        Some(value_t!(matches, "insert-size", usize).unwrap())
                    } else {
                        None
                    },
                    if matches.is_present("std-dev") {
                        Some(value_t!(matches, "std-dev", usize).unwrap())
                    } else {
                        None
                    },
                )
            },
            _ => panic!("No subcommand")
        },
        // This cannot be reached, since the matches step of
        // clap assures that a valid subcommand is provided
        _ => unreachable!(),
    }
}
