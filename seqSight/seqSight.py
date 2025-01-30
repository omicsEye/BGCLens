import os
import sys
import re
import logging
import shutil
import argparse
import warnings
import os, sys

# seqSightdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# print(seqSightdir)
# sys.path.insert(0, seqSightdir)

from seqSight.tools import seqSight_map, seqSight_id


class StoreAsListAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))



# class testseqSightMapOptions:
#     MAX_REF_FILE_SIZE = 4.3e9
#     verbose = False
#     outDir = ""
#     indexDir = ""
#     numThreads = 8
#     outAlignFile = "Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/outAlign.sam"
#     inReadFile = "/Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/ex1.fastq"
#     inReadFilePair1 = ""
#     inReadFilePair2 = ""
#     targetRefFiles = []
#     filterRefFiles = []
#     targetIndexPrefixes = ["Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/genomes"]
#     filterIndexPrefixes = []
#     targetAlignFiles = []
#     filterAlignFiles = []
#     targetAlignParameters = None
#     filterAlignParameters = None
#     btHome = None
#     exp_tag = ""


def parse_arguments():
    parser = argparse.ArgumentParser()

    # MAP
    parser.add_argument('-U', default='', action='store', dest='map_inputread', required=False,
                        help='Input Read Fastq File (Unpaired/Single-end)')
    parser.add_argument('-1', default='', action='store', dest='map_inputread1', required=False,
                        help='Input Read Fastq File (Pair 1)')
    parser.add_argument('-2', default='', action='store', dest='map_inputread2', required=False,
                        help='Input Read Fastq File (Pair 2)')
    parser.add_argument('-targetRefFiles', default='', action='store',
                        dest='map_targetref', required=False,
                        help='Target Reference Genome Fasta Files Full Path (Comma Separated)')
    parser.add_argument('-filterRefFiles', default='', action='store',
                        dest='map_filterref', required=False,
                        help='Filter Reference Genome Fasta Files Full Path (Comma Separated)')
    parser.add_argument('-targetAlignParams', action='store',
                        dest='map_targetalignparams', default=None, required=False,
                        help='Target Mapping Bowtie2 Parameters (Default: seqSight chosen best parameters)')
    parser.add_argument('-filterAlignParams', action='store',
                        dest='map_filteralignparams', default=None, required=False,
                        help='Filter Mapping Bowtie2 Parameters (Default: Use the same Target Mapping Bowtie2 parameters)')
    parser.add_argument('-outDir', action='store', default='.',
                        dest='map_outdir', required=False,
                        help='Output Directory (Default=. (current directory))')
    parser.add_argument('-outAlign', action='store', default='outalign.sam',
                        dest='map_outalign', required=False,
                        help='Output Alignment File Name (Default=outalign.sam)')
    parser.add_argument('-indexDir', action='store', default='.',
                        dest='map_indexdir', required=False,
                        help='Index Directory (Default=. (current directory))')
    parser.add_argument('-targetIndexPrefixes', default='', action='store',
                        dest='map_targetindex', required=False,
                        help='Target Index Prefixes (Comma Separated)')
    parser.add_argument('-filterIndexPrefixes', default='', action='store',
                        dest='map_filterindex', required=False,
                        help='Filter Index Prefixes (Comma Separated)')
    parser.add_argument('-targetAlignFiles', default='', action='store',
                        dest='map_targetalign', required=False,
                        help='Target Alignment Files Full Path (Comma Separated)')
    parser.add_argument('-filterAlignFiles', default='', action='store',
                        dest='map_filteralign', required=False,
                        help='Filter Alignment Files Full Path (Comma Separated)')
    parser.add_argument('-btHome', default=None, action='store',
                        dest='map_bthome', required=False,
                        help='Full Path to Bowtie2 binary directory (Default: Uses bowtie2 in system path)')
    parser.add_argument('-numThreads', action='store', dest='map_numthreads', required=False,
                        default=8, type=int,
                        help='Number of threads to use by aligner (bowtie2) if different from default (8)')
    parser.add_argument('-expTag', action='store', default='', dest='map_exp_tag',
                        help='Experiment Tag added to files generated for identification (Default: seqSightMap)')

    # parser.add_argument('--input', '-i', help="files contains the sequences", type=str, required=True)
    # parser.add_argument('--output', '-o', help="output directory to write teh results", type=str, required=True)

    # ID
    parser.add_argument('--outMatrix', action='store_true', default=False, dest='id_out_matrix',
                          help='Output alignment matrix')
    parser.add_argument('--noUpdatedAlignFile', action='store_true', default=False,
                          dest='id_noalign', help='Do not generate an updated alignment file')
    parser.add_argument('--noDisplayCutoff', action='store_true', default=False,
                          dest='id_nocutoff', help='Do not cutoff display of genomes, even if it is insignificant')
    parser.add_argument('-scoreCutoff', action='store', default=0.01, type=float,
                          dest='id_score_cutoff', help='Score Cutoff')
    parser.add_argument('-emEpsilon', action='store', default=1e-7, type=float,
                          dest='id_emEpsilon', help='EM Algorithm Epsilon cutoff')
    parser.add_argument('-maxIter', action='store', default=50, type=int,
                          dest='id_maxIter', help='EM Algorithm maximum iterations')
    parser.add_argument('-piPrior', action='store', default=0, type=int, dest='id_piPrior',
                          help='EM Algorithm Pi Prior equivalent to adding n unique reads (Default: n=0)')
    parser.add_argument('-thetaPrior', action='store', default=0, type=int, dest='id_thetaPrior',
                          help='EM Algorithm Theta Prior equivalent to adding n non-unique reads (Default: n=0)')
    parser.add_argument('-fileType', action='store', default='sam', dest='id_ali_format',
                          help='Alignment Format: sam/bl8/gnu-sam (Default: sam)')
    # parser.add_argument('-alignFile', action='store', dest='id_ali_file', required=True,
    #                       help='Alignment file path')
    parser.add_argument('-reportDir', action='store', default='.',
                        dest='id_reportDir', required=True,
                        help='Output dir for the relative abundance report(Default=. (current directory))')

    parser.add_argument('-filename', action='store', default='Report',
                        dest='id_filename', required=False,
                        help='Relative abundance table name')

    return parser.parse_args()


def main():
    # Parse arguments from command line
    warnings.filterwarnings("ignore")
    warnings.simplefilter('ignore')

    args = parse_arguments()
    print(args)

    class seqSightMapOptions:
        MAX_REF_FILE_SIZE = 4.3e9
        verbose = True
        outDir = args.map_outdir
        indexDir = args.map_indexdir
        numThreads = args.map_numthreads
        outAlignFile = args.map_outalign
        inReadFile = args.map_inputread
        inReadFilePair1 = args.map_inputread1
        inReadFilePair2 = args.map_inputread2

        if len(args.map_targetref) > 0:
            targetRefFiles = args.map_targetref.split(',')
        else:
            targetRefFiles = []

        if len(args.map_filterref) > 0:
            filterRefFiles = args.map_filterref.split(',')
        else:
            filterRefFiles = []

        if len(args.map_targetindex) > 0:
            targetIndexPrefixes = args.map_targetindex.split(',')
        else:
            targetIndexPrefixes = []

        if len(args.map_filterindex) > 0:
            filterIndexPrefixes = args.map_filterindex.split(',')
        else:
            filterIndexPrefixes = []

        if len(args.map_targetalign) > 0:
            targetAlignFiles = args.map_targetalign.split(',')
        else:
            targetAlignFiles = []

        if len(args.map_filteralign) > 0:
            filterAlignFiles = args.map_filteralign.split(',')
        else:
            filterAlignFiles = []

        targetAlignParameters = args.map_targetalignparams
        filterAlignParameters = args.map_filteralignparams
        btHome = args.map_bthome
        exp_tag = args.map_exp_tag

    #     filterRefFiles = []
    #     targetIndexPrefixes = ["Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/genomes"]
    #     filterIndexPrefixes = []
    #     targetAlignFiles = []
    #     filterAlignFiles = []
    """
    step1 - MAP function;
    Input -  ex1.fastq;
    Output - ex1.sam;
    Tool - bowtie2;
    """
    print(vars(seqSightMapOptions))
    outAlignFile_test = seqSight_map.processseqSightMap(seqSightMapOptions)
    print('MAP done')
    """
     step2 - ID function;
     Input -  sam file from the MAP function;
     Output - report.tsv;
     Tool - reassign+em;
     """
    # seqSight_id.seqSight_reassign(True, 0.01, "testset", "sam", outAlignFile_test,
    #                               "/Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData",
    #                               10, not (False), 0, 0, False, False, emEpsilon=1e-7)


    seqSight_id.seqSight_reassign(out_matrix=args.id_out_matrix,
                                  scoreCutoff=args.id_score_cutoff,
                                  expTag=args.map_exp_tag,
                                  ali_format=args.id_ali_format,
                                  ali_file=outAlignFile_test,
                                  maxIter=args.id_maxIter,
                                  upalign=args.id_noalign,
                                  piPrior=args.id_piPrior,
                                  thetaPrior=args.id_thetaPrior,
                                  noCutOff=args.id_nocutoff,
                                  verbose=False,
                                  emEpsilon=args.id_emEpsilon,
                                  output=args.id_reportDir,
                                  filename=args.id_filename)

                                  # output=args.map_outdir,

    # seqSight_id.seqSight_reassign(out_matrix, scoreCutoff, expTag, ali_format, ali_file, output, maxIter,
    # upalign, piPrior, thetaPrior, noCutOff, verbose, emEpsilon=0.01)

    return print('Reassign done!')


# main()
if __name__ == "__main__":
    main()
