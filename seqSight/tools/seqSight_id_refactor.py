#!/usr/bin/env python
import csv
import math
import os
import re
import sys
import argparse

from seqSight.tools.seqSightReport import seqSightReport
from seqSight.tools.utils import samUtils, seqSightUtils


def conv_align2GRmat(aliDfile, pScoreCutoff, aliFormat):
    in1 = open(aliDfile, 'r')
    U = {}
    NU = {}
    h_readId = {}
    h_refId = {}
    genomes = []
    read = []
    gCnt = 0
    rCnt = 0

    maxScore = None
    minScore = None
    for ln in in1:
        if (ln[0] == '@' or ln[0] == '#'):
            continue

        l = ln.split('\t')

        readId = l[0]
        if (aliFormat == 0 or aliFormat == 1):  # gnu-sam or sam
            if int(l[1]) & 0x4 == 4:  # bitwise FLAG - 0x4 : segment unmapped
                continue
            refId = l[2]
        elif (aliFormat == 2):  # bl8
            refId = l[1]

        if refId == '*':
            continue

        # refId=refId.split("ti:")[-1]
        mObj = re.search(r'ti\|(\d+)\|org\|([^|]+)\|', refId)
        if mObj:
            refId = "ti|" + mObj.group(1) + "|org|" + mObj.group(2)
        else:
            mObj = re.search(r'ti\|(\d+)\|', refId)
            if mObj and mObj.group(1) != "-1":
                refId = "ti|" + mObj.group(1)

        (pScore, skipFlag) = find_entry_score(ln, l, aliFormat, pScoreCutoff)

        if skipFlag:
            continue
        if ((maxScore == None) or (pScore > maxScore)):
            maxScore = pScore
        if ((minScore == None) or (pScore < minScore)):
            minScore = pScore

        gIdx = h_refId.get(refId, -1)
        if gIdx == -1:
            gIdx = gCnt
            h_refId[refId] = gIdx
            genomes.append(refId)
            gCnt += 1

        rIdx = h_readId.get(readId, -1)
        if rIdx == -1:
            # hold on this new read
            # first, wrap previous read profile and see if any previous read has a same profile with that!
            rIdx = rCnt
            h_readId[readId] = rIdx
            read.append(readId)
            rCnt += 1
            U[rIdx] = [[gIdx], [pScore], [float(pScore)], pScore]
        else:
            if (rIdx in U):
                if gIdx in U[rIdx][0]:
                    continue
                NU[rIdx] = U[rIdx]
                del U[rIdx]
            if gIdx in NU[rIdx][0]:
                continue
            NU[rIdx][0].append(gIdx)
            NU[rIdx][1].append(pScore)
            if pScore > NU[rIdx][3]:
                NU[rIdx][3] = pScore
    #			length = len(NU[rIdx][1])
    #			NU[rIdx][2] = [1.0/length]*length

    in1.close()

    if (aliFormat == 1):  # sam
        (U, NU) = samUtils.rescale_samscore(U, NU, maxScore, minScore)
    del h_refId, h_readId
    for rIdx in U:
        U[rIdx] = [U[rIdx][0][0], U[rIdx][1][0]]  # keep gIdx and score only
    for rIdx in NU:
        pScoreSum = sum(NU[rIdx][1])
        NU[rIdx][2] = [k / pScoreSum for k in NU[rIdx][1]]  # Normalizing pScore

    return U, NU, genomes, read


def initialize_probabilities(genomes):
    """
    Initialize the probabilities for the EM algorithm.

    Parameters:
    genomes (list): List of genomes.

    Returns:
    tuple: Tuple containing initialized probabilities for pi and theta.
    """
    G = len(genomes)
    initial_value = 1. / G
    pi = [initial_value for _ in genomes]
    theta = [initial_value for _ in genomes]
    return pi, theta


def update_pi(pisum, priorWeight, total_weight, piPrior, num_genomes):
    """
    Update the pi values in the M-step of the EM algorithm.

    Parameters:
    pisum (list): Summation of theta values for each genome.
    total_weight (float): Total weight of reads.
    piPrior (float): Prior for the pi parameter.
    num_genomes (int): Number of genomes.

    Returns:
    list: Updated pi values.
    """
    pip = piPrior * priorWeight
    # pi = [(k + pip) / (total_weight + pip * num_genomes) for k in pisum]
    pi = [(k + pip) / (total_weight + pip * num_genomes) for k in pisum]  ## update pi
    return pi


def update_theta(thetasum, total_weight, thetaPrior, num_genomes):
    """
    Update the theta values in the M-step of the EM algorithm.

    Parameters:
    thetasum (list): Summation of pi values for each genome.
    total_weight (float): Total weight of reads.
    thetaPrior (float): Prior for the theta parameter.
    num_genomes (int): Number of genomes.

    Returns:
    list: Updated theta values.
    """
    thetap = thetaPrior * total_weight
    theta = [(k + thetap) / (total_weight + thetap * num_genomes) for k in thetasum]
    return theta


def has_converged(current_pi, previous_pi, emEpsilon, verbose):
    """
    Check if the EM algorithm has converged based on the change in pi values.

    Parameters:
    current_pi (list): Current pi values after the latest iteration.
    previous_pi (list): Pi values from the previous iteration.
    emEpsilon (float): Convergence threshold.
    verbose (bool): Flag to enable verbose logging.

    Returns:
    bool: True if the algorithm has converged, False otherwise.
    """
    print(current_pi)
    print(previous_pi)
    # Calculate the total change in pi values between iterations
    total_change = sum(abs(current - previous) for current, previous in zip(current_pi, previous_pi))

    # Log the change if verbose mode is enabled
    if verbose:
        print(f"[Iteration]: Total change in pi values: {total_change}")

    # Check if the change is below the threshold
    return total_change <= emEpsilon


def seqSight_em(U, NU, genomes, maxIter, emEpsilon, verbose, piPrior, thetaPrior):
    """
    Perform the Expectation-Maximization algorithm on the given data.

    Parameters:
    U (dict): Dictionary containing unique reads data.
    NU (dict): Dictionary containing non-unique reads data.
    genomes (list): List of genomes.
    maxIter (int): Maximum number of iterations.
    emEpsilon (float): Convergence threshold.
    verbose (bool): Verbose output flag.
    piPrior (float): Prior for pi parameter.
    thetaPrior (float): Prior for theta parameter.

    Returns:
    tuple: Tuple containing initial and final pi and theta values, and updated NU.
    """
    # Initialize probabilities
    pi, theta = initialize_probabilities(genomes)
    initPi = pi

    # Get weights for unique and non-unique reads
    Uweights, Utotal = get_weights(U, weight_index=1)
    NUweights, NUtotal = get_weights(NU, weight_index=3)

    priorWeight = max(max(Uweights), max(NUweights))
    pisum0 = calculate_pisum0(U, genomes)
    print("pisum0", pisum0)  # [5.376234283632271e+43, 0, 0]
    lenNU = len(NU)
    if lenNU == 0:
        lenNU = 1

    for iteration in range(maxIter):
        pi_old = pi
        thetasum = [0 for _ in genomes]

        # E-Step: Update NU and calculate thetasum
        thetasum = perform_e_step(NU, pi, theta, thetasum)

        # M-Step: Update pi and theta
        pisum = [thetasum[k] + pisum0[k] for k in range(len(thetasum))]  ### calculate tally for pi
        pi = update_pi(pisum, priorWeight, Utotal + NUtotal, piPrior, len(genomes))
        theta = update_theta(thetasum, NUtotal, thetaPrior, len(genomes))

        # Check for convergence
        cutoff = 0.0
        for k in range(len(pi)):
            cutoff += abs(pi_old[k] - pi[k])
        if (cutoff <= emEpsilon or lenNU == 1):
            break
        # if has_converged(pi, iteration, emEpsilon, verbose):
        #     break

    return initPi, pi, theta, NU


def get_weights(reads_dict, weight_index):
    """
    Get weights from a dictionary of reads.

    Parameters:
    reads_dict (dict): Dictionary of reads.
    weight_index (int): Index for the weight in the reads data.

    Returns:
    tuple: Tuple containing the list of weights and the total weight.
    """
    weights = [read[weight_index] for read in reads_dict.values() if read]
    total_weight = sum(weights)
    return weights, total_weight


def calculate_pisum0(U, genomes):
    """
    Calculate the initial summation of pi values for unique reads.

    Parameters:
    U (dict): Dictionary containing unique reads data.

    Returns:
    list: List of summed pi values for each genome.
    """
    pisum0 = [0 for _ in genomes]
    for read in U.values():
        pisum0[read[0]] += read[1]
    return pisum0


def perform_e_step(NU, pi, theta, thetasum):
    """
    Perform the E-step of the EM algorithm.

    Parameters:
    NU (dict): Dictionary containing non-unique reads data.
    pi (list): List of current pi values.
    theta (list): List of current theta values.
    thetasum (list): List to accumulate theta
    Returns:
    list: Updated thetasum after the E-step.
    """
    print("pi", pi)
    for j in NU:  # for each non-unique read, j
        z = NU[j]
        genome_indices = z[0]  # a set of any genome mapping with j
        print("genome_indices", genome_indices)
        pitmp = [pi[k] for k in genome_indices]  # relevant pis for the read
        print("pitmp", pitmp)
        thetatmp = [theta[k] for k in genome_indices]  # relevant thetas for the read
        xtmp = [pitmp[k] * thetatmp[k] * z[1][k] for k in range(len(genome_indices))]
        xsum = sum(xtmp)
        xnorm = [k / xsum if xsum != 0 else 0.0 for k in xtmp]  # normalized new xs

        NU[j][2] = xnorm  # update x in NU

        for idx, genome_idx in enumerate(genome_indices):
            thetasum[genome_idx] += xnorm[idx] * NU[j][3]  # weighted tally for theta
    print("thetasum",thetasum) #thetasum [8.064351425448407e+43, 9.40311461505841e+19, 9.40311461505841e+19]


    return thetasum


def seqSight_reassign(out_matrix, scoreCutoff, expTag, ali_format, ali_file, output, maxIter,
                      upalign, piPrior, thetaPrior, noCutOff, verbose, emEpsilon=0.01):
    """
    Reassigns sequences using the EM algorithm and generates a report.

    Parameters:
    out_matrix (bool): Whether to output the matrix.
    scoreCutoff (float): Score cutoff for alignment.
    expTag (str): Tag for the experiment.
    ali_format (str): Format of the alignment file.
    ali_file (str): Path to the alignment file.
    output (str): Output directory path.
    maxIter (int): Maximum iterations for the EM algorithm.
    upalign (bool): Whether to update the alignment.
    piPrior (float): Prior value for pi in the EM algorithm.
    thetaPrior (float): Prior value for theta in the EM algorithm.
    noCutOff (bool): If True, no cutoff will be applied in reporting.
    verbose (bool): Verbose output flag.
    emEpsilon (float): Convergence epsilon for the EM algorithm.

    Returns:
    tuple: A tuple containing various outputs including the final report.
    """
    validate_alignment_file(ali_file)
    aliFormat = determine_alignment_format(ali_format, verbose)
    U, NU, genomes, reads = conv_align2GRmat(ali_file, scoreCutoff, aliFormat)
    print("U", U)
    print("NU", NU)
    print("genomes", genomes)
    print("reads", reads)
    nG = len(genomes)
    nR = len(reads)

    if verbose:
        print("Starting EM iteration...")
        print(f"(Genomes, Reads) = {len(genomes)}x{len(reads)}")
        print("Delta Change:")

    if out_matrix:
        initial_align_output(genomes, reads, U, NU, expTag, ali_file, output)

    (bestHitInitialReads, bestHitInitial, level1Initial, level2Initial) = \
        seqSightReport.computeBestHit(U, NU, genomes, reads)

    initPi, pi, _, NU = seqSight_em(U, NU, genomes, maxIter, emEpsilon, verbose, piPrior, thetaPrior)

    (bestHitFinalReads, bestHitFinal, level1Final, level2Final) = \
        seqSightReport.computeBestHit(U, NU, genomes, reads)
    # finalReport = output + os.sep + expTag + '-' + ali_format + '-report.tsv'
    finalReport = output + os.sep + "seqSight/Test/TestData/TestRefactor" + '-' + ali_format + '-report.tsv'

    print("finalReport", finalReport)
    header = ['Genome', 'Final Guess', 'Final Best Hit', 'Final Best Hit Read Numbers', \
              'Final High Confidence Hits', 'Final Low Confidence Hits', 'Initial Guess', \
              'Initial Best Hit', 'Initial Best Hit Read Numbers', \
              'Initial High Confidence Hits', 'Initial Low Confidence Hits']
    (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = seqSightReport.write_tsv_report(
        finalReport, nR, nG, pi, genomes, initPi, bestHitInitial, bestHitInitialReads,
        bestHitFinal, bestHitFinalReads, level1Initial, level2Initial, level1Final,
        level2Final, header, noCutOff)

    print("initPi", initPi)
    print("pi", pi)
    print("_", _)
    print("NU", NU)

    reAlignfile = ali_file
    if upalign:
        reAlignfile = rewrite_align(U, NU, ali_file, scoreCutoff, aliFormat, output)

    print("x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11", x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)

    #return report_file, reAlignfile
    return (finalReport, x2, x3, x4, x5, x1, x6, x7, x8, x9, x10, x11, reAlignfile)


def validate_alignment_file(ali_file):
    if os.stat(ali_file).st_size < 1.0:
        raise ValueError(f'The alignment file {ali_file} is empty.')


def determine_alignment_format(ali_format, verbose):
    format_map = {'gnu-sam': 0, 'sam': 1, 'bl8': 2}
    if ali_format not in format_map:
        raise ValueError(f"Unknown alignment format: {ali_format}")

    if verbose:
        print(f"Parsing {ali_format} file/likelihood score/reads and mapped genomes...")
    return format_map[ali_format]


def initial_align_output(ref, read, U, NU, expTag, ali_file, output):
    print("ref", ref)
    genomeId = output + os.sep + expTag + '-genomeId.txt'
    oFp = open(genomeId, 'w')
    csv_writer = csv.writer(oFp, delimiter='\n')
    csv_writer.writerows([ref])
    oFp.close()

    readId = output + os.sep + expTag + '-readId.txt'
    oFp = open(readId, 'w')
    print("read", read)
    csv_writer = csv.writer(oFp, delimiter='\n')
    csv_writer.writerows([read])
    oFp.close()


def rewrite_align(U, NU, aliDfile, pScoreCutoff, aliFormat, output):
    seqSightUtils.ensure_dir(output)
    f = os.path.basename(aliDfile)
    reAlignfile = output + os.sep + 'updated_' + f

    with open(reAlignfile, 'w') as of:
        with open(aliDfile, 'r') as in1:
            h_readId = {}
            h_refId = {}
            genomes = []
            read = []
            gCnt = 0
            rCnt = 0

            mxBitSc = 700
            sigma2 = 3
            for ln in in1:
                if (ln[0] == '@' or ln[0] == '#'):
                    of.write(ln)
                    continue

                l = ln.split('\t')

                readId = l[0]
                if (aliFormat == 0 or aliFormat == 1):  # gnu-sam or sam
                    # refId=l[2].split("ti:")[-1]
                    refId = l[2]
                    if int(l[1]) & 0x4 == 4:  # bitwise FLAG - 0x4 : segment unmapped
                        continue
                elif (aliFormat == 2):  # bl8
                    refId = l[1]

                if refId == '*':
                    continue

                mObj = re.search(r'ti\|(\d+)\|org\|([^|]+)\|', refId)
                if mObj:
                    refId = "ti|" + mObj.group(1) + "|org|" + mObj.group(2)
                else:
                    mObj = re.search(r'ti\|(\d+)\|', refId)
                    if mObj and mObj.group(1) != "-1":
                        refId = "ti|" + mObj.group(1)

                (_, skipFlag) = find_entry_score(ln, l, aliFormat, pScoreCutoff)
                if skipFlag:
                    continue

                gIdx = h_refId.get(refId, -1)
                if gIdx == -1:
                    gIdx = gCnt
                    h_refId[refId] = gIdx
                    genomes.append(refId)
                    gCnt += 1

                rIdx = h_readId.get(readId, -1)
                if rIdx == -1:
                    # hold on this new read
                    # first, wrap previous read profile and see if any previous read has a same profile with that!
                    rIdx = rCnt
                    h_readId[readId] = rIdx
                    read.append(readId)
                    rCnt += 1
                    if rIdx in U:
                        of.write(ln)
                        continue

                if rIdx in NU:
                    if (aliFormat == 0):  # gnu-sam
                        scoreComponents = l[12].split(':')
                        (upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
                        scoreComponents[2] = str(upPscore * pscoreSum)
                        if (scoreComponents[2] < pScoreCutoff):
                            continue
                        l[12] = ':'.join(scoreComponents)
                        ln = '\t'.join(l)
                        of.write(ln)
                    elif (aliFormat == 1):  # sam
                        (upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
                        if (upPscore < pScoreCutoff):
                            continue
                        if (upPscore >= 1.0):
                            upPscore = 0.999999
                        mapq2 = math.log10(1 - upPscore)
                        l[4] = str(int(round(-10.0 * mapq2)))
                        ln = '\t'.join(l)
                        of.write(ln)
                    elif (aliFormat == 2):  # bl8
                        (upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
                        score = upPscore * pscoreSum
                        if score <= 0.0:
                            continue
                        bitSc = math.log(score)
                        if bitSc > mxBitSc:
                            bitSc = mxBitSc
                        l[10] = str(bitSc * sigma2)
                        ln = '\t'.join(l)
                        of.write(ln)

    return reAlignfile


def find_updated_score(NU, rIdx, gIdx):
    try:
        index = NU[rIdx][0].index(gIdx);
    except ValueError:
        print('Value Error: %s' % gIdx)
        return (0., 0.)
    pscoreSum = 0.0
    for pscore in NU[rIdx][1]:
        pscoreSum += pscore
    pscoreSum /= 100
    upPscore = NU[rIdx][2][index]
    return (upPscore, pscoreSum)


# ===========================================================
# Internal function to calculate the score from the alignment file entries
def find_entry_score(ln, l, aliFormat, pScoreCutoff):
    mxBitSc = 700
    sigma2 = 3
    skipFlag = False
    if (aliFormat == 0):  # gnu-sam
        pScore = float(l[12].split(':')[2])
        if (pScore < pScoreCutoff):
            skipFlag = True
    elif (aliFormat == 1):  # sam
        pScore = samUtils.findSamAlignScore(l)
        if pScore is None:
            mapq = float(l[4])
            mapq2 = mapq / (-10.0)
            pScore = 1.0 - pow(10, mapq2)
            if (pScore < pScoreCutoff):
                skipFlag = True
    elif (aliFormat == 2):  # bl8
        eVal = float(l[10])
        if (eVal > pScoreCutoff):
            skipFlag = True
        bitSc = float(l[11]) / sigma2
        if bitSc > mxBitSc:
            bitSc = mxBitSc
        pScore = math.exp(bitSc)
    # pScore = int(round(pScore*100)) # Converting to integer to conserve memory space
    # if pScore < 1:
    #	skipFlag = true
    return (pScore, skipFlag)
