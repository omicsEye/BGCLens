#!/usr/bin/env python

"""
This code has been adopted from humann DOI: https://doi.org/10.1038/s41592-018-0176-y
Join a set of gene, pathway, or taxonomy tables into a single table
This module will join gene and pathway tables output by HUMAnN.
Dependencies: Biom (only required if running with .biom files)
To Run:
$ ./join_tables.py -i <input_dir> -o <gene_table.{tsv,biom}>

"""

import argparse
import sys
import os
import pandas as pd

GENE_TABLE_DELIMITER = "\t"


def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="Join BGCLens Mapping tables\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the directory of tables\n",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="the table to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join tables with this string included in the file name")
    parser.add_argument(
        "--colname",
        default="Final Best Hit",
        help="Name of the column from the input tables that you want to use in join_tables function\n",
        required=False)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args = parse_arguments(sys.argv)

    # check for format of the gene tables
    input_dir = os.path.abspath(args.input)

    print(input_dir)

    # check for format of the output table
    output_dir = os.path.abspath(args.output)

    print(output_dir)

    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." +
                 "  Please enter a new directory.")

    path_list2 = os.listdir(input_dir)
    # print(path_list2)

    # Set the first one as base
    samp2 = pd.read_csv(input_dir + "/" + path_list2[0], delimiter='\t', header=0, skiprows=1)
    # Final best Hit change to colname
    s2 = samp2[["Genome", args.colname]]
    lst = ["Genome_ID", path_list2[0].split(".tsv")[0]]
    s2.columns = lst

    # Apply for the rest others
    for i in path_list2[1:]:
        if i != ".DS_Store":
            print(i)
            samp = pd.read_csv(input_dir + "/" + i, delimiter='\t', header=0, skiprows=1)
            temp = samp[["Genome", args.colname]]
            renamelst = ["Genome_ID", i.split(".tsv")[0]]
            temp.columns = renamelst
            s2 = pd.merge(s2, temp, how='outer', on=['Genome_ID'])
    # pd.set_option('display.max_rows',None)
    # print(s2)
    s2.set_index('Genome_ID', inplace=True)
    s2['count'] = s2.apply(lambda x: x.count(), axis=1)
    s2.loc["count1"] = s2.count()

    # s2.to_csv(output_dir + "SampeToFullGenome0610.csv")  # 可以指定文件目录路径
    s2.to_csv(output_dir)


if __name__ == "__main__":
    main()
