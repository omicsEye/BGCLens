#!/usr/bin/python
import os, sys, math, shutil

seqSightdir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, seqSightdir)

from seqSight.tools.utils import seqParse
from seqSight.tools.bowtie2Wrap import bowtie2Wrap


# Main entry function to seqSightMap that does all the processing
def processseqSightMap(seqSightMapOptions):
    procseqSightMapOptions = copyseqSightMapOptions(seqSightMapOptions)
    # Create the output folder if it not exists
    if os.path.isdir(seqSightMapOptions.outDir):
        if seqSightMapOptions.verbose:
            print("The output folder" + seqSightMapOptions.outDir + "already exists.")
    else:
        os.mkdir(seqSightMapOptions.outDir)
        if seqSightMapOptions.verbose:
            print("Creating the output folder" + seqSightMapOptions.outDir)

    # Splitting reference files if bigger than MAX_REF_FILE_SIZE
    ptargetRefFiles = []
    for filePath in seqSightMapOptions.targetRefFiles:
        if seqSightMapOptions.verbose:
            print("Checking whether the file: " + filePath + " needs to be split")
        files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE);
        for f in files:
            ptargetRefFiles.append(f)
    procseqSightMapOptions.targetRefFiles = ptargetRefFiles
    pfilterRefFiles = []
    for filePath in seqSightMapOptions.filterRefFiles:
        if seqSightMapOptions.verbose:
            print("Checking whether the file: " + filePath + " needs to be split")
        files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE);
        for f in files:
            pfilterRefFiles.append(f)
    procseqSightMapOptions.filterRefFiles = pfilterRefFiles

    # Creating Index if it does not exist
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procseqSightMapOptions.verbose
    bowtie2Options.btHome = procseqSightMapOptions.btHome
    bowtie2Options.indexDir = procseqSightMapOptions.indexDir
    for filePath in ptargetRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if seqSightMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        procseqSightMapOptions.targetIndexPrefixes.append(base)
    for filePath in pfilterRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if seqSightMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        procseqSightMapOptions.filterIndexPrefixes.append(base)

    # Creating the Alignment file
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procseqSightMapOptions.verbose
    bowtie2Options.btHome = procseqSightMapOptions.btHome
    bowtie2Options.numThreads = procseqSightMapOptions.numThreads
    bowtie2Options.outDir = procseqSightMapOptions.outDir
    bowtie2Options.indexDir = procseqSightMapOptions.indexDir
    bowtie2Options.readFile = procseqSightMapOptions.inReadFile
    bowtie2Options.readFilePair1 = procseqSightMapOptions.inReadFilePair1
    bowtie2Options.readFilePair2 = procseqSightMapOptions.inReadFilePair2
    if (len(procseqSightMapOptions.inReadFilePair1) > 0 and
            len(procseqSightMapOptions.inReadFilePair2) > 0 and
            len(procseqSightMapOptions.inReadFile) > 0):
        bowtie2Options.bothReadFlag = True  # newly added
    elif (len(procseqSightMapOptions.inReadFilePair1) > 0 and
          len(procseqSightMapOptions.inReadFilePair2) > 0):
        bowtie2Options.pairedReadFlag = True
    elif (len(procseqSightMapOptions.inReadFile) > 0):
        bowtie2Options.singleReadFlag = True  # newly added
    if procseqSightMapOptions.targetAlignParameters is not None:
        bowtie2Options.additionalOptions = procseqSightMapOptions.targetAlignParameters
    for indexPrefix in procseqSightMapOptions.targetIndexPrefixes:
        bowtie2Options.btIndexPrefix = procseqSightMapOptions.indexDir + os.sep + indexPrefix
        bowtie2Options.outAlignFile = procseqSightMapOptions.exp_tag + indexPrefix + ".sam"
        if seqSightMapOptions.verbose:
            print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
        bowtie2Wrap.run_bowtie2(bowtie2Options)
        # procseqSightMapOptions.targetAlignFiles.append(bowtie2Options.outAlignFile)
        # previous
        # procseqSightMapOptions.targetAlignFiles.append(procseqSightMapOptions.outDir + os.sep +
        #                                                bowtie2Options.outAlignFile)
        # latest
        procseqSightMapOptions.targetAlignFiles.append(bowtie2Options.outAlignFile)
        # print("procseqSightMapOptions.outDir + os.sep + bowtie2Options.outAlignFile",procseqSightMapOptions.outDir + os.sep + bowtie2Options.outAlignFile)

    # Appending the Alignment files and Filtering
    if len(procseqSightMapOptions.targetAlignFiles) > 1:
        appendAlignFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.exp_tag + "appendAlign.sam"
        if seqSightMapOptions.verbose:
            print("Appending alignment files to: " + appendAlignFile)
        append_sam_file(appendAlignFile, procseqSightMapOptions.targetAlignFiles)
    else:
        appendAlignFile = procseqSightMapOptions.targetAlignFiles[0]

    if len(procseqSightMapOptions.filterIndexPrefixes) > 0:
        bowtie2Options.readFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.exp_tag + "appendAlign.fq"
        bowtie2Options.readFilePair1 = ""
        bowtie2Options.readFilePair2 = ""
        bowtie2Options.bothReadFlag = False
        bowtie2Options.pairedReadFlag = False
        bowtie2Options.singleReadFlag = True
        if procseqSightMapOptions.filterAlignParameters is not None:
            bowtie2Options.additionalOptions = procseqSightMapOptions.filterAlignParameters
        bowtie2Wrap.extractRead(appendAlignFile, bowtie2Options.readFile)
        for indexPrefix in procseqSightMapOptions.filterIndexPrefixes:
            bowtie2Options.btIndexPrefix = procseqSightMapOptions.indexDir + os.sep + indexPrefix
            bowtie2Options.outAlignFile = procseqSightMapOptions.exp_tag + indexPrefix + ".sam"
            if seqSightMapOptions.verbose:
                print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
            bowtie2Wrap.run_bowtie2(bowtie2Options)
            procseqSightMapOptions.filterAlignFiles.append(procseqSightMapOptions.outDir + os.sep +
                                                           bowtie2Options.outAlignFile)
    # Filtering the Alignment file
    outAlignFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.outAlignFile
    if seqSightMapOptions.verbose:
        print("Filtering and creating the alignment: " + outAlignFile)
    if len(procseqSightMapOptions.filterAlignFiles) > 0:
        filter_alignment(appendAlignFile, procseqSightMapOptions.filterAlignFiles, outAlignFile)
    elif ((len(procseqSightMapOptions.targetAlignFiles) > 1) or \
          (len(procseqSightMapOptions.targetIndexPrefixes) > 0)):
        os.rename(appendAlignFile, outAlignFile)
    else:  # Input appendAlignFile provided by user, hence make a copy for outAlignFile
        shutil.copy(appendAlignFile, outAlignFile)
    return outAlignFile


# Make a copy of core seqSightMapOptions
def copyseqSightMapOptions(seqSightMapOptions):
    procseqSightMapOptions = seqSightMapOptions()
    procseqSightMapOptions.verbose = seqSightMapOptions.verbose
    procseqSightMapOptions.outDir = seqSightMapOptions.outDir
    procseqSightMapOptions.indexDir = seqSightMapOptions.indexDir
    procseqSightMapOptions.numThreads = seqSightMapOptions.numThreads
    procseqSightMapOptions.outAlignFile = seqSightMapOptions.outAlignFile
    procseqSightMapOptions.inReadFile = seqSightMapOptions.inReadFile
    procseqSightMapOptions.inReadFilePair1 = seqSightMapOptions.inReadFilePair1
    procseqSightMapOptions.inReadFilePair2 = seqSightMapOptions.inReadFilePair2
    procseqSightMapOptions.targetRefFiles = seqSightMapOptions.targetRefFiles
    procseqSightMapOptions.filterRefFiles = seqSightMapOptions.filterRefFiles
    procseqSightMapOptions.targetIndexPrefixes = seqSightMapOptions.targetIndexPrefixes
    procseqSightMapOptions.filterIndexPrefixes = seqSightMapOptions.filterIndexPrefixes
    procseqSightMapOptions.targetAlignFiles = seqSightMapOptions.targetAlignFiles
    procseqSightMapOptions.filterAlignFiles = seqSightMapOptions.filterAlignFiles
    procseqSightMapOptions.targetAlignParameters = seqSightMapOptions.targetAlignParameters
    procseqSightMapOptions.filterAlignParameters = seqSightMapOptions.filterAlignParameters
    procseqSightMapOptions.btHome = seqSightMapOptions.btHome
    procseqSightMapOptions.exp_tag = seqSightMapOptions.exp_tag
    return procseqSightMapOptions


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
