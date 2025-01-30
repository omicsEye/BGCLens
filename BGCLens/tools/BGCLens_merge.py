import pandas as pd
import argparse
import csv
import sys
import os


def merge(file_path, outputDirFile):
    # Get all file names from the result folder
    path_list = os.listdir(file_path)
    path_name = []
    for i in path_list:
        path_name.append(i.split("-sam")[0])
    print(sorted(path_name))

    # Set the first one as base
    samp2 = pd.read_csv(file_path+sorted(path_name)[0]+"-sam-report.tsv", delimiter='\t', header=0, skiprows=1)
    s2 = samp2[["Genome", "Final Guess"]]
    firstFile = str(sorted(path_name)[0])
    s2.columns = ["Genome_ID",firstFile]
    print(s2)

    # Apply for the rest others
    for i in sorted(path_name)[1:]:
        if i != ".DS_Store":
            print(i)
            samp = pd.read_csv(file_path + i + "-sam-report.tsv", delimiter='\t',
                               header=0, skiprows=1)
            temp = samp[["Genome", "Final Guess"]]
            renamelst = ["Genome_ID", i]
            temp.columns = renamelst
            print(temp)
            s2 = pd.merge(s2, temp, how='outer', on=['Genome_ID'])
    # pd.set_option('display.max_rows',None)
    print(s2)
    s2.set_index('Genome_ID', inplace=True)
    # s2['count'] = s2.apply(lambda x: x.count(), axis=1)
    # s2.loc["count1"] = s2.count()

    # put the path you wish to save the output
    # s2.to_csv("./BGCLens_BGC_All_Results/BGC_Frankel_Merged.tsv")
    s2.to_csv(outputDirFile)

    return


argp = argparse.ArgumentParser(prog="merge_BGCLens_tables.py",
                               description="Performs a table join on one or more BGCLens output files.")
argp.add_argument("--file_path", help="Name of file containing the paths to the files to combine")
argp.add_argument('--outputDirFile', metavar="output.txt", help="Name of output file containing the path")

argp.usage = (argp.format_usage() + "\nPlease make sure to supply file paths to the files to combine.\n\n" +
              "A wildcard to indicate all .txt files that start with Table can be used as follows:\n" +
              "    ./merge_BGCLens_tables.py Table*.txt > output.txt")


# def main():
#     args = argp.parse_args()
#
#     if args.l:
#         args.aistms = [x.strip().split()[0] for x in open(args.l)]
#
#     if not args.aistms:
#         print('merge_BGCLens_tables: no inputs to merge!')
#         return
#
#     if args.o and os.path.exists(args.o) and not args.overwrite:
#         print('merge_BGCLens_tables: output file "{}" exists, specify the --overwrite param to ovrewrite it!'.format(
#             args.o))
#         return
#
#     merge(args.aistms, open(args.o, 'w') if args.o else sys.stdout, args.gtdb_profiles)


if __name__ == '__main__':
    file_path = "/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/BGCLens_BGC_All_Results/BGC_Frankel/"
    outputDirFile = "/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/BGC_Frankel_MergedTEST.tsv"
    merge(file_path, outputDirFile)



