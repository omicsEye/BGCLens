#!/usr/bin/python
import os, sys, math, shutil

BGCLensdir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, BGCLensdir)

from BGCLens.tools.utils import seqParse
from BGCLens.tools.bowtie2Wrap import bowtie2Wrap


# Main entry function to BGCLensMap that does all the processing
def processBGCLensMap(BGCLensMapOptions):
    procBGCLensMapOptions = copyBGCLensMapOptions(BGCLensMapOptions)
    # Create the output folder if it not exists
    if os.path.isdir(BGCLensMapOptions.outDir):
        if BGCLensMapOptions.verbose:
            print("The output folder" + BGCLensMapOptions.outDir + "already exists.")
    else:
        os.mkdir(BGCLensMapOptions.outDir)
        if BGCLensMapOptions.verbose:
            print("Creating the output folder" + BGCLensMapOptions.outDir)

    # Splitting reference files if bigger than MAX_REF_FILE_SIZE
    ptargetRefFiles = []
    for filePath in BGCLensMapOptions.targetRefFiles:
        if BGCLensMapOptions.verbose:
            print("Checking whether the file: " + filePath + " needs to be split")
        files = splitCheck(filePath, BGCLensMapOptions.MAX_REF_FILE_SIZE);
        for f in files:
            ptargetRefFiles.append(f)
    procBGCLensMapOptions.targetRefFiles = ptargetRefFiles
    pfilterRefFiles = []
    for filePath in BGCLensMapOptions.filterRefFiles:
        if BGCLensMapOptions.verbose:
            print("Checking whether the file: " + filePath + " needs to be split")
        files = splitCheck(filePath, BGCLensMapOptions.MAX_REF_FILE_SIZE);
        for f in files:
            pfilterRefFiles.append(f)
    procBGCLensMapOptions.filterRefFiles = pfilterRefFiles

    # Creating Index if it does not exist
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procBGCLensMapOptions.verbose
    bowtie2Options.btHome = procBGCLensMapOptions.btHome
    bowtie2Options.indexDir = procBGCLensMapOptions.indexDir
    for filePath in ptargetRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if BGCLensMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        procBGCLensMapOptions.targetIndexPrefixes.append(base)
    for filePath in pfilterRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if BGCLensMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        procBGCLensMapOptions.filterIndexPrefixes.append(base)

    # Creating the Alignment file
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procBGCLensMapOptions.verbose
    bowtie2Options.btHome = procBGCLensMapOptions.btHome
    bowtie2Options.numThreads = procBGCLensMapOptions.numThreads
    bowtie2Options.outDir = procBGCLensMapOptions.outDir
    bowtie2Options.indexDir = procBGCLensMapOptions.indexDir
    bowtie2Options.readFile = procBGCLensMapOptions.inReadFile
    bowtie2Options.readFilePair1 = procBGCLensMapOptions.inReadFilePair1
    bowtie2Options.readFilePair2 = procBGCLensMapOptions.inReadFilePair2
    if (len(procBGCLensMapOptions.inReadFilePair1) > 0 and
            len(procBGCLensMapOptions.inReadFilePair2) > 0 and
            len(procBGCLensMapOptions.inReadFile) > 0):
        bowtie2Options.bothReadFlag = True  # newly added
    elif (len(procBGCLensMapOptions.inReadFilePair1) > 0 and
          len(procBGCLensMapOptions.inReadFilePair2) > 0):
        bowtie2Options.pairedReadFlag = True
    elif (len(procBGCLensMapOptions.inReadFile) > 0):
        bowtie2Options.singleReadFlag = True  # newly added
    if procBGCLensMapOptions.targetAlignParameters is not None:
        bowtie2Options.additionalOptions = procBGCLensMapOptions.targetAlignParameters
    for indexPrefix in procBGCLensMapOptions.targetIndexPrefixes:

        # previous
        bowtie2Options.btIndexPrefix = procBGCLensMapOptions.indexDir + os.sep + indexPrefix
        # make the indexprefix files path without affected by indexDir
        # latest
        # bowtie2Options.btIndexPrefix = indexPrefix
        print("bowtie2Options.btIndexPrefix",bowtie2Options.btIndexPrefix)
        bowtie2Options.outAlignFile = procBGCLensMapOptions.exp_tag + indexPrefix + ".sam"
        if BGCLensMapOptions.verbose:
            print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
        bowtie2Wrap.run_bowtie2(bowtie2Options)
        # procBGCLensMapOptions.targetAlignFiles.append(bowtie2Options.outAlignFile)
        # previous
        procBGCLensMapOptions.targetAlignFiles.append(procBGCLensMapOptions.outDir + os.sep +
                                                        bowtie2Options.outAlignFile)
        # latest
        # procBGCLensMapOptions.targetAlignFiles.append(bowtie2Options.outAlignFile)
        # print("procBGCLensMapOptions.outDir + os.sep + bowtie2Options.outAlignFile",procBGCLensMapOptions.outDir + os.sep + bowtie2Options.outAlignFile)

    # Appending the Alignment files and Filtering
    if len(procBGCLensMapOptions.targetAlignFiles) > 1:
        appendAlignFile = procBGCLensMapOptions.outDir + os.sep + procBGCLensMapOptions.exp_tag + "appendAlign.sam"
        if BGCLensMapOptions.verbose:
            print("Appending alignment files to: " + appendAlignFile)
        append_sam_file(appendAlignFile, procBGCLensMapOptions.targetAlignFiles)
    else:
        appendAlignFile = procBGCLensMapOptions.targetAlignFiles[0]

    if len(procBGCLensMapOptions.filterIndexPrefixes) > 0:
        bowtie2Options.readFile = procBGCLensMapOptions.outDir + os.sep + procBGCLensMapOptions.exp_tag + "appendAlign.fq"
        bowtie2Options.readFilePair1 = ""
        bowtie2Options.readFilePair2 = ""
        bowtie2Options.bothReadFlag = False
        bowtie2Options.pairedReadFlag = False
        bowtie2Options.singleReadFlag = True
        if procBGCLensMapOptions.filterAlignParameters is not None:
            bowtie2Options.additionalOptions = procBGCLensMapOptions.filterAlignParameters
        bowtie2Wrap.extractRead(appendAlignFile, bowtie2Options.readFile)
        for indexPrefix in procBGCLensMapOptions.filterIndexPrefixes:
            bowtie2Options.btIndexPrefix = procBGCLensMapOptions.indexDir + os.sep + indexPrefix
            bowtie2Options.outAlignFile = procBGCLensMapOptions.exp_tag + indexPrefix + ".sam"
            if BGCLensMapOptions.verbose:
                print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
            bowtie2Wrap.run_bowtie2(bowtie2Options)
            procBGCLensMapOptions.filterAlignFiles.append(procBGCLensMapOptions.outDir + os.sep +
                                                           bowtie2Options.outAlignFile)
    # Filtering the Alignment file
    outAlignFile = procBGCLensMapOptions.outDir + os.sep + procBGCLensMapOptions.outAlignFile
    if BGCLensMapOptions.verbose:
        print("Filtering and creating the alignment: " + outAlignFile)
    if len(procBGCLensMapOptions.filterAlignFiles) > 0:
        filter_alignment(appendAlignFile, procBGCLensMapOptions.filterAlignFiles, outAlignFile)
    elif ((len(procBGCLensMapOptions.targetAlignFiles) > 1) or \
          (len(procBGCLensMapOptions.targetIndexPrefixes) > 0)):
        print("appendAlignFile:",appendAlignFile,"outAlignFile",outAlignFile)
        # os.rename(appendAlignFile, outAlignFile)
        os.rename("/Users/xinyang/Library/CloudStorage/Box-Box/BGC_Learning/bgc_genes_nt.sam", outAlignFile)
    else:  # Input appendAlignFile provided by user, hence make a copy for outAlignFile
        shutil.copy(appendAlignFile, outAlignFile)
    return outAlignFile


# Make a copy of core BGCLensMapOptions
def copyBGCLensMapOptions(BGCLensMapOptions):
    procBGCLensMapOptions = BGCLensMapOptions()
    procBGCLensMapOptions.verbose = BGCLensMapOptions.verbose
    procBGCLensMapOptions.outDir = BGCLensMapOptions.outDir
    procBGCLensMapOptions.indexDir = BGCLensMapOptions.indexDir
    procBGCLensMapOptions.numThreads = BGCLensMapOptions.numThreads
    procBGCLensMapOptions.outAlignFile = BGCLensMapOptions.outAlignFile
    procBGCLensMapOptions.inReadFile = BGCLensMapOptions.inReadFile
    procBGCLensMapOptions.inReadFilePair1 = BGCLensMapOptions.inReadFilePair1
    procBGCLensMapOptions.inReadFilePair2 = BGCLensMapOptions.inReadFilePair2
    procBGCLensMapOptions.targetRefFiles = BGCLensMapOptions.targetRefFiles
    procBGCLensMapOptions.filterRefFiles = BGCLensMapOptions.filterRefFiles
    procBGCLensMapOptions.targetIndexPrefixes = BGCLensMapOptions.targetIndexPrefixes
    procBGCLensMapOptions.filterIndexPrefixes = BGCLensMapOptions.filterIndexPrefixes
    procBGCLensMapOptions.targetAlignFiles = BGCLensMapOptions.targetAlignFiles
    procBGCLensMapOptions.filterAlignFiles = BGCLensMapOptions.filterAlignFiles
    procBGCLensMapOptions.targetAlignParameters = BGCLensMapOptions.targetAlignParameters
    procBGCLensMapOptions.filterAlignParameters = BGCLensMapOptions.filterAlignParameters
    procBGCLensMapOptions.btHome = BGCLensMapOptions.btHome
    procBGCLensMapOptions.exp_tag = BGCLensMapOptions.exp_tag
    return procBGCLensMapOptions


# If the given file size is greater than maxSize, then it splits
# and returns a list of split file paths where each file is less than maxSize
def splitCheck(filePath, maxSize):
    files = []
    print("FILE PATH:", filePath)
    fileSize = os.stat(filePath).st_size
    nSplit = 1
    if (fileSize > maxSize):
        nSplit = int(math.ceil(1.0 * fileSize / float(maxSize)))
    if nSplit == 1:
        files.append(filePath)
        return files
    (base, ext) = os.path.splitext(filePath)
    # check if we have already done this splitting
    for i in range(nSplit):
        fiPath = base + '_' + str(i) + ext
        splitReq = False
        if not os.path.exists(fiPath):
            splitReq = True
            break
    fps = []
    for i in range(nSplit):
        fiPath = base + '_' + str(i) + ext
        files.append(fiPath)
        if splitReq:
            fps.append(open(fiPath, 'w'))
    if splitReq:
        with open(filePath, 'r') as fp:
            j = 0
            if ext == '.fq':
                for r in seqParse.parse(fp, 'fastq'):
                    fps[j % nSplit].write('>%s %s\n%s\n%s\n' % (r.id, r.description, r.seq, r.qual))
                    j += 1
            else:
                for r in seqParse.parse(fp, 'fasta'):
                    fps[j % nSplit].write('>%s %s\n%s\n' % (r.id, r.description, r.seq))
                    j += 1
        for i in range(nSplit):
            fps[i].close()
    return files


def filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile):
    return bowtie2Wrap.filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile)


# Appends all the appendFiles to outfile
def append_file(outfile, appendFiles):
    with open(outfile, 'w') as out1:
        for file1 in appendFiles:
            if (file1 is not None):
                with open(file1, 'r') as in2:
                    for ln in in2:
                        out1.write(ln)


# Appends all the sam appendFiles to outfile
def append_sam_file(outfile, appendFiles):
    with open(outfile, 'w') as out1:
        # First, writing the header by merging headers from all files
        for file1 in appendFiles:
            if (file1 is not None):
                with open(file1, 'r') as in2:
                    for ln in in2:
                        if ln[0] == '@':
                            out1.write(ln)
        # Writing the body by merging body from all files
        for file1 in appendFiles:
            if (file1 is not None):
                with open(file1, 'r') as in2:
                    for ln in in2:
                        if ln[0] != '@':
                            out1.write(ln)
