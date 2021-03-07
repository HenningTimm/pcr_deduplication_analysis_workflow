use bio::io::fastq;
use std::fs;
use std::process::Command;

/// Compare an output file to the expected output and delete the output file.
fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
        .arg(result)
        .arg(expected)
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    fs::remove_file(result).unwrap();
}

/// Compare two fastq files, ignoring the name lines
/// Reads are sorted by their sequence, which is not 100% robust
/// if mutations/ sequencing errors are considered.
fn compare_fastq(result: &str, expected: &str) {
    let result_reader = fastq::Reader::from_file(result).unwrap();
    let mut result_recs: Vec<fastq::Record> =
        result_reader.records().filter_map(Result::ok).collect();
    result_recs.sort_by_key(|x| x.seq().to_owned());
    let expected_reader = fastq::Reader::from_file(expected).unwrap();
    let mut expected_recs: Vec<fastq::Record> =
        expected_reader.records().filter_map(Result::ok).collect();
    expected_recs.sort_by_key(|x| x.seq().to_owned());

    for (result, expected) in result_recs.iter().zip(expected_recs.iter()) {
        assert_eq!(result.seq(), expected.seq());
        assert_eq!(result.qual(), expected.qual());
    }
}

#[test]
fn test_call_consensus_reads_two_cluster() {
    assert!(
        Command::new("bash")
                .arg("-c")
                .arg("target/debug/rbt call-consensus-reads fastq --umi-len 3 -u --max-umi-dist 0 --max-seq-dist 2 tests/test-consensus.fastq tests/test-consensus.fastq /tmp/test-consensus.1.fastq /tmp/test-consensus.2.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test-consensus.1.fastq",
        "tests/expected/test-consensus.1.fastq",
    );
    compare_fastq(
        "/tmp/test-consensus.2.fastq",
        "tests/expected/test-consensus.2.fastq",
    );
}

#[test]
fn test_call_consensus_reads_single_cluster() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt call-consensus-reads fastq --umi-len 3 -u --max-umi-dist 2 --max-seq-dist 2 tests/test-consensus.fastq tests/test-consensus.fastq /tmp/test-consensus_single.1.fastq /tmp/test-consensus_single.2.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test-consensus_single.1.fastq",
        "tests/expected/test-consensus_single.1.fastq",
    );
    compare_fastq(
        "/tmp/test-consensus_single.2.fastq",
        "tests/expected/test-consensus_single.2.fastq",
    );
}

#[test]
fn test_call_overlapping_consensus_reads() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt call-consensus-reads fastq --umi-len 10 --max-umi-dist 0 --max-seq-dist 8 --insert-size 450 --std-dev 50  tests/overlapping-consensus.1.fastq tests/overlapping-consensus.2.fastq /tmp/test_overlapping-consensus.1.fastq /tmp/test_overlapping-consensus.2.fastq /tmp/test_overlapping-consensus.3.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test_overlapping-consensus.1.fastq",
        "tests/expected/test_overlapping-consensus.1.fastq",
    );
    compare_fastq(
        "/tmp/test_overlapping-consensus.2.fastq",
        "tests/expected/test_overlapping-consensus.2.fastq",
    );
    compare_fastq(
        "/tmp/test_overlapping-consensus.3.fastq",
        "tests/expected/test_overlapping-consensus.3.fastq",
    );
}
