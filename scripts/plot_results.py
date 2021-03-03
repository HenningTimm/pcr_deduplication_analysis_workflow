import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re


def tuftefy(ax):
    """Remove spines and tick position markers to reduce ink."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color('grey')

    ax.grid(color="#a0a0a0", alpha=0.5)
    ax.get_yaxis().grid(True)
    ax.get_xaxis().grid(False)


def plot_bars(data):
    # barchart
    fig, ax = plt.subplots()
    
    width = 0.35
    
    x_vals, labels = zip(*enumerate(data["locus"]))
    
    x_vals_stacked = [x - width / 2 for x in x_vals]
    x_vals_dedup = [x + width / 2 for x in x_vals]
    rects1 = ax.bar(x_vals_stacked, data["real"] + data["pcr_copies"], width=width, alpha=0.8, label="ref")
    rects2 = ax.bar(x_vals_stacked, data["real"], width=width, alpha=0.8, label="ref")
    rects3 = ax.bar(x_vals_dedup, data["after_dedup"], width=width, alpha=0.8, label="dedup")
    tuftefy(ax)
    plt.savefig(snakemake.output.loci_pdf)


def plot_violins(tidy_df):
    sns.set(font_scale=1.5, style="whitegrid")

    tidy_df = tidy_df.rename(columns={
        "Max. sequence dist.": "Max. sequence distance",
        "Max. DBR dist.": "Max. DBR distance",
    })

    cat = sns.catplot(x="type", y="count", row="Max. sequence distance",
                     col="Max. DBR distance", kind="violin", aspect=2.2,
                     data=tidy_df[tidy_df["type"] != "pcr_copies"])
    cat.set_axis_labels("Read Types", "Locus Coverage")
    
    sns.despine()
    plt.savefig(snakemake.output.violins_pdf)


def get_params(filename):
    filename = os.path.basename(filename)
    translation = {
        "ud": "max UMI dist",
        "sd": "max sequence dist",
        }
    params = dict()
    for key, value in re.findall("_(\S{2}):(\S)", filename):
        params[translation[key]] = value
    return params


data = None
tidy_data = None

for csv_file in snakemake.input.csvs:
    
    parameter_set_data = pd.read_csv(csv_file)
    parameter_set_data["real+pcr"] = parameter_set_data["real"] + parameter_set_data["pcr_copies"]

    tidy_df = pd.melt(parameter_set_data, id_vars=["locus"], value_vars=["real", "pcr_copies", "real+pcr", "after_dedup"], var_name="type", value_name="count")

    params = get_params(csv_file)
    tidy_df["Max. DBR dist."] = params["max UMI dist"]
    tidy_df["Max. sequence dist."] = params["max sequence dist"]
    parameter_set_data["Max. DBR dist."] = params["max UMI dist"]
    parameter_set_data["Max. sequence dist."] = params["max sequence dist"]

    if data is None:
        data = parameter_set_data
    else:
        data = data.append(parameter_set_data, ignore_index=True)

    if tidy_data is None:
        tidy_data = tidy_df
    else:
        tidy_data = tidy_data.append(tidy_df, ignore_index=True)

with open(snakemake.output.wide_df, "w") as wide_df:
    wide_df.write(data.to_csv(index=False))

plot_violins(tidy_data)

