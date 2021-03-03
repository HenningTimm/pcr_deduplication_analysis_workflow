import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_scatter():
    wide_df = pd.read_csv(snakemake.input.wide_df)

    wide_df = wide_df.rename(columns={
        "Max. sequence dist.": "Max. sequence distance",
        "Max. DBR dist.": "Max. DBR distance",
    })


    g = sns.FacetGrid(wide_df, row="Max. sequence distance", col="Max. DBR distance",
                      height=2.5, aspect=2.2)
    g = g.map(plt.scatter, "real", "after_dedup", alpha=0.25)
    g.set_axis_labels("Simulated Locus Coverage", "Locus Coverage\nAfter Deduplication")
    print(g, dir(g))
    for row in g.axes:
        for x in row:
            x.plot([0, 40], [0, 40], alpha=0.75, color="r")
    
    sns.set(font_scale=0.5)
    sns.despine()
    plt.savefig(snakemake.output.scatter, dpi=200)

plot_scatter()
