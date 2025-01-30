#!/usr/bin/python

import os, sys


from BGCLens.tools.utils import samUtils


# ===========================================================
class Bowtie2Options:
    DEFAULT_OPTION = "--very-sensitive-local -k 100 --score-min L,20,1.0"  # Another option --score-min L,205,0.0
    btHome = None  # Uses bowtie2 in the system path by default# None
    verbose = True
    btBin = "bowtie2"
    btIndexer = "bowtie2-build"
    pairedReadFlag = False
    singleReadFlag = False
    bothReadFlag = False
    readFile = ""
    readFilePair1 = ""
    readFilePair2 = ""
    outAlignFile = ""
    outDir = ""  # previous:"."
    indexDir = "."
    numThreads = 8
    additionalOptions = "--very-sensitive-local -k 100 --score-min L,20,1.0"#DEFAULT_OPTION
    refFile = ""
    btIndexPrefix = ""

    def __init__(self): {}


# Runs bowtie2 alignment with the given bowtie2 options
def run_bowtie2(bowtie2Options):
    align_exist = 1
    if os.path.exists(bowtie2Options.outAlignFile):
        if bowtie2Options.verbose:
            print("Bowtie2 alignment file already exist: " + bowtie2Options.outAlignFile)
        return align_exist

    btBinPath = bowtie2Options.btBin
    # TODO BOWTIE2 BEHOME
    # if bowtie2Options.btHome is not None:
    #     print(bowtie2Options.btHome)
    #     btBinPath = bowtie2Options.btHome + os.sep + bowtie2Options.btBin
    #     print("2btBinPath", btBinPath)
    print("bowtie2Options.outDir", bowtie2Options.outDir)
    print("bowtie2Options.outAlignFile",bowtie2Options.outAlignFile)
    # latest
    outAlignFilePath = bowtie2Options.outAlignFile
    # previous
    # outAlignFilePath = bowtie2Options.outDir + os.sep + bowtie2Options.outAlignFile
    print("outAlignFilePath", outAlignFilePath)
    print("bowtie2Options.DEFAULT_OPTION", bowtie2Options.DEFAULT_OPTION)
    print("btIndexPrefix",bowtie2Options.btIndexPrefix)
    if bowtie2Options.bothReadFlag:  # newly added
        cmd = "%s -x %s -1 %s -2 %s -U %s -p %s %s -S %s" % (btBinPath,
                                                             bowtie2Options.btIndexPrefix,
                                                             bowtie2Options.readFilePair1, bowtie2Options.readFilePair2, bowtie2Options.readFile,
                                                             bowtie2Options.numThreads, bowtie2Options.DEFAULT_OPTION,
                                                             outAlignFilePath)
        print("1", cmd)
    if bowtie2Options.pairedReadFlag:
        cmd = "%s -x %s -1 %s -2 %s -p %s %s -S %s" % (btBinPath,
                                                       bowtie2Options.btIndexPrefix,
                                                       bowtie2Options.readFilePair1, bowtie2Options.readFilePair2,
                                                       bowtie2Options.numThreads, bowtie2Options.DEFAULT_OPTION,
                                                       outAlignFilePath)
        print("2", cmd)
    elif bowtie2Options.singleReadFlag:  # newly added
        cmd = "%s -x %s -U %s -p %s %s -S %s" % (btBinPath,
                                                 bowtie2Options.btIndexPrefix,
                                                 bowtie2Options.readFile,
                                                 bowtie2Options.numThreads, bowtie2Options.DEFAULT_OPTION,
                                                 outAlignFilePath)
        print("3", cmd)
    if bowtie2Options.verbose:
        print(cmd)
    return os.system(cmd)


# Create bowtie2 Index for the given reference file
def create_bowtie2_index(bowtie2Options):
    index_exist = 1
    btIndexPrefixPath = bowtie2Options.indexDir + os.sep + bowtie2Options.btIndexPrefix
    btIndexPath = btIndexPrefixPath + ".1.bt2"
    btIndexPath_large = btIndexPrefixPath + ".1.bt2l"
    if os.path.exists(btIndexPath) or os.path.exists(btIndexPath_large):
        if bowtie2Options.verbose:
            print("Bowtie2 index already exist for: " + btIndexPath)
        return index_exist
    btIndexerPath = bowtie2Options.btIndexer
    # TODO BOWTIE2 BEHOME
    # if bowtie2Options.btHome is not None:
    #     btIndexerPath = bowtie2Options.btHome + os.sep + bowtie2Options.btIndexer
    cmd = "%s %s %s" % (btIndexerPath, bowtie2Options.refFile, btIndexPrefixPath)
    if bowtie2Options.verbose:
        print(cmd)
    return os.system(cmd)


# Returns an alignment file with the reads that align in targetAlignFile and
# that do not align in filterAlignFile
def filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile):
    print("targetAlignFile", targetAlignFile)
    print("filterAlignFiles", filterAlignFiles)
    print("outAlignFile", outAlignFile)
    if filterAlignFiles[0] is not None:
        h_readId = find_readsAlignScore(filterAlignFiles)
    else:
        h_readId = {}
    with open(targetAlignFile, 'r') as in1:
        with open(outAlignFile, 'w') as out1:
            for ln in in1:
                if (ln[0] == '@' or ln[0] == '#'):
                    out1.write(ln)
                    continue
                l = ln.split('\t')
                readId = l[0]
                aScore = samUtils.findSamAlignHiScore(l)
                if h_readId is not {} :
                    if (aScore is not None):
                        score = h_readId.get(readId, None)
                        if (score == None) or (score < aScore):
                            # This read is (not/having low scores) in the filterAlignFiles
                            out1.write(ln)
    return outAlignFile


# Returns a hash table of all read ids to alignment score
# from all the given alignment files
# Used by filter_sam() function
def find_readsAlignScore(filterAlignFiles):
    h_readId = {}
    for filterAlignFile in filterAlignFiles:
        with open(filterAlignFile, 'r') as in1:
            for ln in in1:
                if (ln[0] == '@' or ln[0] == '#'):
                    continue
                l = ln.split('\t')
                readId = l[0]
                aScore = samUtils.findSamAlignHiScore(l)
                if aScore is not None:
                    score = h_readId.get(readId, None)
                    if (score == None) or (score < aScore):
                        h_readId[readId] = aScore
    return h_readId


def extractRead(appendAlignFile, readFile):
    h_readId = {}
    with open(appendAlignFile, 'r') as in1:
        with open(readFile, 'w') as out1:
            for ln in in1:
                if (ln[0] == '@' or ln[0] == '#'):
                    continue
                l = ln.split('\t')
                aScore = samUtils.findSamAlignScore(l)
                if (aScore is None):
                    continue
                readId = l[0]
                score = h_readId.get(readId, None)
                if (score is None):
                    h_readId[readId] = aScore
                    sequence = l[9]
                    quality = l[10]
                    out1.write('@%s\n%s\n+\n%s\n' % (readId, sequence, quality))
