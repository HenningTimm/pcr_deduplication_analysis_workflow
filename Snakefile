rule all:
    input:
        "plots/ddRAGEdataset_violins.pdf",
        "plots/ddRAGEdataset_scatter.png",


bin_name = "rbt_with_names"
rbt_path = "rust-bio-tools"


# Compile a local version of rust bio tools
# that is able to interpret ddRAGe annotations
# in name lines.
rule compile_custom_rbt:
    input:
        f"../{rbt_path}/src/main.rs",
    output:
        rbt_bin=bin_name,
    conda:
        "envs/rust.yaml"
    params:
        target_path=f"{rbt_path}/target/release/rbt",
        manifest_path=f"{rbt_path}/Cargo.toml",
    shell:
        "cargo build --release --manifest-path {params.manifest_path} && "
        "cp {params.target_path} {output.rbt_bin}"


# Use ddrage to simulate a test dataset with labled PCR duplicates
rule simulate_pcr_duplicates:
    input:
        "barcode.txt"
    output:
        "data/{filename}_ATCACG_1.fastq",
        "data/{filename}_ATCACG_2.fastq",
        "data/{filename}_ATCACG_gt.yaml",
    conda:
        "envs/testdata.yaml"
    params:
        loci=100000,
        individuals=1,
        output_prefix="data",
        hrl_number=0,
        prob_ID=0,
        prob_het=0,
        event_probs="0.95 0.0 0.05",  #--event-probabilities 0.9 0.00 0.05 common, dropout, mut, this prevents dropouts
    shell:
        "ddrage --name {wildcards.filename} -l {params.loci} -n {params.individuals} -b barcode.txt -o {params.output_prefix} "
        "--hrl-number {params.hrl_number} --prob-incomplete-digestion {params.prob_ID} "
        "--prob-heterozygous {params.prob_het} --event-probabilities {params.event_probs} "
        "--no-singletons"


# Deduplicate the generated dataset using modified rust bio tools
rule dedup:
    input:
        fq1="data/{filename}_ATCACG_1.fastq",
        fq2="data/{filename}_ATCACG_2.fastq",
        rbt_bin=bin_name,
    output:
        fq1="dedup/{filename}_ud:{umi_max_dist}_sd:{max_seq_dist}_1.fq.gz",
        fq2="dedup/{filename}_ud:{umi_max_dist}_sd:{max_seq_dist}_2.fq.gz",
    conda:
        "envs/consensus.yaml"
    params:
        umi_length=13,
    shell:
        "./{input.rbt_bin} call-consensus-reads fastq -l {params.umi_length} "
        "-d {wildcards.umi_max_dist} -D {wildcards.max_seq_dist} "
        "--umi-on-reverse --verbose-read-names "
        "{input.fq1} {input.fq2} {output.fq1} {output.fq2}"


# Compare original fq file, deduplicated fq file and ground truth
# to evaluate how consensus reads were assembled
rule evaluate:
    input:
        gt_file="data/{filename}_ATCACG_gt.yaml",
        fq1="data/{filename}_ATCACG_1.fastq",
        fq1_dedup="dedup/{filename}_ud:{umi_max_dist}_sd:{max_seq_dist}_1.fq.gz",
    output:
        csv="results/{filename}_ud:{umi_max_dist}_sd:{max_seq_dist}_results.csv"
    conda:
        "envs/eval.yaml"
    script:
        "scripts/evaluate_deduplication.py"


rule plot_locus_results:
    input:
        csvs=expand("results/{{filename}}_ud:{umi_max_dist}_sd:{max_seq_dist}_results.csv",
                   umi_max_dist=[1, 2],
                   max_seq_dist=[1, 2, 3, 4, 6, 8],
        )
    output:
        violins_pdf="plots/{filename}_violins.pdf",
        wide_df="plots/{filename}_wide_df.csv",
    conda:
        "envs/eval.yaml"
    script:
        "scripts/plot_results.py"


rule plot_scatter:
    input:
        wide_df="plots/{filename}_wide_df.csv",
    output:
        scatter="plots/{filename}_scatter.png",
    conda:
        "envs/eval.yaml"
    script:
        "scripts/plot_scatter.py"
