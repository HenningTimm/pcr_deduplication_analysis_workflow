from collections import namedtuple, defaultdict, Counter
import pandas as pd
import dinopy as dp
import re
import pysam


class PCRRecord:

    def __init__(self, ):
        self.pcr = 0
        self.real = 0

    def __repr__(self):
        return f"real: {self.real:>10} pcr copies: {self.pcr:>10}"






def parse_info_line(inf_line):
    """parse a ddrage name line into a dictionary
    """
    line_info = {
        "pcr_copy": False,
    }

    # pattern that matches all ddRAGE read annotations
    pattern = "(\S+):'(\S*?\s??\S*?)'"
    # assemble a dictionary with all ddRAGE read infos
    for (key, value) in re.findall(pattern, inf_line):
        if key == "type" and value == "PCR copy":
            line_info["pcr_copy"] = True
            continue
        line_info[key] = value
    return line_info


def parse_fq_file(fq_file):

    pcr_counts = defaultdict(PCRRecord)
    fqr = dp.FastqReader(fq_file)

    for read in fqr.reads():

        line_info = parse_info_line(read.name.decode())
        try:
            locus = line_info["at_locus"]
        except:
            print(read, line_info)
            raise

        # count the number of real an PCR reads for this locus
        # tested against a solution with grep + wc -l
        # grep  ^@ data/ddRAGEdataset_ATCACG_1.fastq | grep "at_locus:'1'" | grep -v PCR | wc -l
        if line_info["pcr_copy"]:
            pcr_counts[locus].pcr += 1
        else:
            pcr_counts[locus].real += 1

    return pcr_counts


def get_dedup_coverages(fq_file):

    singleton_pattern = "at_locus:'(\d+)'"
    consensus_pattern = "Locus_(\d+)"
    locus_coverage_counter = Counter()

    fqr = dp.FastqReader(fq_file)
    for r in fqr.reads():
        name = r.name.decode()
        singleton = re.search(singleton_pattern, name)
        consensus = re.findall(consensus_pattern, name)
        if singleton is not None:
            locus = singleton.groups()[0]
            locus_coverage_counter[int(locus)] += 1
        if consensus:
            merged_loci = [int(s) for s in consensus]
            # print(merged_loci)
            if len(set(merged_loci)) > 1:
                raise ValueError("Overmerge!")
            else:
                locus_coverage_counter[list(set(merged_loci))[0]] += 1
        if singleton is None and not consensus:
            raise ValueError("BAD NAME. No matches")
        if singleton is not None  and consensus:
            raise ValueError("BAD NAME. Two matches")
    # print(sorted(locus_coverage_counter.items()))
    loci, coverage = zip(*locus_coverage_counter.items())
    df = pd.DataFrame({
        "locus": loci,
        "after_dedup": coverage
    })
    return df



# def compare_locus_numbers(pcr_counts, bam_file, out_file):
#     # print(pcr_counts)
#     df = pd.DataFrame(columns=["locus", "real", "pcr_copies", "after_dedup"])
#     with pysam.AlignmentFile(bam_file) as sam_file:
#         for locus, counts in pcr_counts.items():
#             if locus.startswith("singleton"):
#                 continue
#             name = f"Locus_{locus}"
#             # print(locus, name)
#             reads_after_dedup = 0
#             for _ in sam_file.fetch(name):
#                 reads_after_dedup += 1
#             # print(counts, reads_after_dedup)
#             df = df.append(
#                 {
#                     "locus": name,
#                     "real": counts.real,
#                     "pcr_copies": counts.pcr,
#                     "after_dedup": reads_after_dedup
#                 },
#                 ignore_index=True,
#             )
#     with open(out_file, "w") as csv_file:
#         csv_file.write(df.to_csv(index=False))


        
def compare_locus_numbers(pcr_counts, dedup_covs, out_file):
    # parse the default dict assembled above
    count_df = pd.DataFrame(columns=["locus", "real", "pcr_copies"])
    for locus, counts in pcr_counts.items():
        if locus.startswith("singleton"):
            continue
        count_df = count_df.append(
            {
                "locus": int(locus),
                "real": counts.real,
                "pcr_copies": counts.pcr,
            },
            ignore_index=True,
        )
    df = pd.merge(dedup_covs, count_df, left_on="locus", right_on="locus")

    with open(out_file, "w") as csv_file:
        csv_file.write(df.to_csv(index=False))


pcr_counts = parse_fq_file(snakemake.input.fq1)
dedup_covs = get_dedup_coverages(snakemake.input.fq1_dedup)
# compare_locus_numbers(pcr_counts, snakemake.input.bam, snakemake.output.csv)
compare_locus_numbers(pcr_counts, dedup_covs, snakemake.output.csv)
