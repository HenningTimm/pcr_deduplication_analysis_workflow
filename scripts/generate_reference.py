import yaml
import dinopy as dp
from collections import namedtuple

GTRecord = namedtuple("GTRecord", ["name", "seq_p5", "seq_p7", "mutations",
                                   "id_reads", "dropout", "alleles",
                                   "allele_frequencies"])
GTStats = namedtuple("GTStats", ["nr_muts", "nr_snps", "nr_inserts",
                                 "nr_deletions", "nr_loci_with_snps",
                                 "nr_loci_with_muts", "nr_loci"])

def parse_rage_gt_file(gt_file):
    """Read in a RAGE ground truth file.

    Returns:
        list: of GTRecord named tuples each of which hash the entries
        'name', 'seq_p5', 'seq_p7', 'mutations', id_reads', 'dropout'
    """
    with open(gt_file, 'r') as stream:
        try:
            # read all documents in the data
            inds, loci, *other = list(yaml.load_all(stream))
        except yaml.YAMLError as exc:
            print(exc)
    nr_muts, nr_snps, nr_inserts, nr_deletions = 0, 0, 0, 0
    nr_loci_with_snps, nr_loci_with_muts = 0, 0

    ind_inf = inds["Individual Information"]
    
    p5_enz = list(inds["Individual Information"].values())[0]["p5 overhang"]
    loc_seqs = []
    # filter out all loci with only one allele, i.e. all unmutated loci
    loci_with_snps = ((n, l) for (n, l) in loci.items()
                      if len(l["allele coverages"]) > 1)

    # print("inds", inds)
    spacer_lengths = [len(i["p5 spacer"]) for i
                      in inds["Individual Information"].values()]
    # spacer_variance = max(spacer_lengths) - min(spacer_lengths)
    overhang_lengths = [len(i["p5 overhang"]) for i
                        in inds["Individual Information"].values()]
    # overhang_variance = max(overhang_lengths) - min(overhang_lengths)
    offset = None

    

    for name, locus in loci.items():
        dropout = []
        mutations = set()
        gt_alleles = {}
        for n, ind in locus["individuals"].items():
            if ind:
                # dropout events get an empty dict,
                # hence everything that does not evaluate to False
                # is a valid entry with one or two alleles
                for nr_allele, allele in ind.items():
                    normalized_mutations = set()
                    mutations |= normalized_mutations  # extend set
                dropout.append(False)
            else:
                dropout.append(True)

        if any((mut_type == "SNP" for mut_type, _ in mutations)):
            nr_loci_with_snps += 1
            nr_loci_with_muts += 1
        elif mutations:
            nr_loci_with_muts += 1  # locus with indels only

        # compile and append a record for this locus
        id_reads = locus["id reads"]
        # gt_record = GTRecord(name, seq, mutations, id_reads, dropout)
        gt_record = GTRecord(name, p5_enz + locus["p5 seq"], locus["p7 seq"],
                             mutations, id_reads, dropout,
                             gt_alleles, locus["allele frequencies"])
        nr_muts += len(mutations)
        nr_snps += len([mut for mut in mutations if ">" in mut])
        nr_inserts += len([mut for mut in mutations if "+" in mut])
        nr_deletions += len([mut for mut in mutations if "-" in mut])
        loc_seqs.append(gt_record)


    gt_stats = GTStats(nr_muts, nr_snps, nr_inserts, nr_deletions,
                       nr_loci_with_snps, nr_loci_with_muts, len(loci))
    return loc_seqs, gt_stats, list(ind_inf.values())[0]


def generate_fasta_file(loci, stats, ind_inf, fa_file):
    print(ind_inf)
    with dp.FastaWriter(fa_file) as faw:
    
        for locus in loci:
            p5_read = f"{ind_inf['p5 bc']}{ind_inf['p5 spacer']}{ind_inf['p5 overhang']}{locus.seq_p5}"
            name = bytes(locus.name.replace(" ", "_"), 'ascii')
            faw.write_entry((p5_read, name), dtype=str)
        

loci, stats, ind_inf = parse_rage_gt_file(snakemake.input.gt_file)
generate_fasta_file(loci, stats, ind_inf, snakemake.output.fa)
