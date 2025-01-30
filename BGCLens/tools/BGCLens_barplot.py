#!/usr/bin/env python

"""
This code has been adopted from humann DOI: https://doi.org/10.1038/s41592-018-0176-y
To Run:
$ python ./BGCLens_barplot.py -i1 <input_dir1> -i2 <input_dir2> -title <"name_of_barplot"> -o <barplot>
$ python ./BGCLens_barplot.py -i1 /Users/xinyangzhang/Downloads/bacterial/all5ResultsNum120.tsv -i2 /Users/xinyangzhang/Downloads/bacterial/FiveTargetReads120.tsv -t "life time omics"  -o ./test_barplot.pdf

"""
import pandas as pd
import os
import sys
import numpy as np
import argparse
import matplotlib.gridspec as gridspec
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE_TABLE_DELIMITER = "\t"

def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="Create a barplot for the sequencing \n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i1", "--input1",
        help="the directory of table for the sample's reads\n",
        required=True)
    parser.add_argument(
        "-i2", "--input2",
        help="the directory of table for the sample's percentage in each dataset\n",
        required=True)
    parser.add_argument(
        "-t", "--title",
        default="Distribution",
        help="the barplot's title\n",
        required=False)
    parser.add_argument(
        "-o", "--output",
        help="the barplot to save\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="output name of the barplot")

    return parser.parse_args()


def main():
    args = parse_arguments(sys.argv)

    # check for format of the tables
    input_dir1 = os.path.abspath(args.input1)
    # all5ResultsNum120.tsv
    print(input_dir1)

    # check for format of the tables
    input_dir2 = os.path.abspath(args.input2)
    # FiveTargetReads120
    print(input_dir2)

    # title_name = os.path.abspath(args.title)
    # print(title_name)

    # check for format of the output table
    output_dir = os.path.abspath(args.output)

    print(output_dir)

    df = pd.read_csv(input_dir1, sep="\t", header=0, index_col=0)
    print(df)

    df = df.sort_values(by=["Unidentified", "Viral", "Bacterial", "Fungi", "Archaea"], ascending=False)
    #TODO df = df.sort_values(by=["Human", "Viral", "Bacterial", "Fungi", "Archaea"], ascending=False)


    df1 = df.div(df.sum(axis=1), axis=0)

    lst_Viral = []
    for i in df1["Viral"]:
        lst_Viral.append(i)
    # print(lst_Viral)

    lst_Bacterial = []
    for i in df1["Bacterial"]:
        lst_Bacterial.append(i)
    # print(lst_Bacterial)

    lst_Human = []
    # TODO -  for i in df1["Human"]:
    for i in df1["Unidentified"]:
        lst_Human.append(i)
    # print(lst_Human)

    lst_Fungi = []
    for i in df1["Fungi"]:
        lst_Fungi.append(i)
    # print(lst_Fungi)

    lst_Archaea = []
    for i in df1["Archaea"]:
        lst_Archaea.append(i)
    # print(lst_Archaea)

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    # positions of the left bar-boundaries
    bar_l = [i + 1 for i in range(len(df.index))]
    plt.subplots_adjust(wspace=0, hspace=0)

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    reads = pd.read_csv(input_dir2, sep="\t", header=0, index_col=0)

    sort_list = list(df.index)

    sort_df_human = reads.loc[sort_list]
    x2 = len(list(sort_df_human.index))
    y2 = list(sort_df_human['#Total Reads'])

    # Create bars
    ax1.bar(bar_l, y2, width=1, alpha=0.99, linewidth=0, label='Reads', color="#AFAFAF")
    ax1.set_xticks([])
    ax1.set_yticks(np.linspace(0, 6559849, 4))
    ax1.set_ylabel("Reads")

    ax1.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)

    # Set the bar width
    bar_width = 1

    # Create a bar plot, in position bar_1
    ax2.bar(bar_l,
            lst_Human,
            # set the width
            width=bar_width,
            # with pre_score on the bottom
            # with the label mid score
            label='Unidentified',  # hypothetical unannotated',
            # with alpha 0.5
            alpha=0.99,
            # with color
            color='#0081a7',  # 'silver',#'#F1911E',
            linewidth=0)
    ax2.bar(bar_l,
            lst_Viral,
            # set the width
            width=bar_width,
            # with pre_score on the bottom
            # bottom=df['Unannotated'],
            # with the label mid score
            label='Viral',  # hypothetical unannotated',
            bottom=lst_Human,
            # with alpha 0.5
            alpha=0.99,
            # with color
            color='#00afb9',  # 'silver',#'#F1911E',
            linewidth=0)

    # Create a bar plot, in position bar_1
    ax2.bar(bar_l,
            lst_Bacterial,
            # set the width
            width=bar_width,
            # with the label pre score
            label='Bacterial',
            bottom=[a + b for a, b in zip(lst_Human, lst_Viral)],
            # with alpha 0.5
            alpha=0.99,
            # with color
            color='#fdfcdc',  # '#29CF66',
            linewidth=0)

    # Create a bar plot, in position bar_1
    ax2.bar(bar_l,
            lst_Fungi,
            # set the width
            width=bar_width,
            # with the label pre score
            label='Fungi',
            bottom=[a + b + c for a, b, c in zip(lst_Human, lst_Viral, lst_Bacterial)],
            # with alpha 0.5
            alpha=0.99,
            # with color
            color='#fed9b7',  # '#29CF66',
            linewidth=0)

    # Create a bar plot, in position bar_1
    ax2.bar(bar_l,
            lst_Archaea,
            # set the width
            width=bar_width,
            # with the label pre score
            label='Archaea',
            bottom=[a + b + c + d for a, b, c, d in zip(lst_Human, lst_Viral, lst_Bacterial, lst_Fungi)],
            # with alpha 0.5
            alpha=0.99,
            # with color
            color='#f07167',  # '#29CF66',
            linewidth=0)

    ax2.set_yticklabels(["0.000", "20%", "40%", "60%", "80%", "100%"])
    ax2.invert_yaxis()
    ax2.set_ylabel("Fraction")
    ax2.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    plt.suptitle(args.title)
    plt.xlabel("Sample N=" + str(x2))
    plt.savefig(output_dir, dpi=100, bbox_inches="tight")


if __name__ == "__main__":
    main()
