#!/usr/bin/env python
# -*- coding: utf-8 -*

###
# read a file of deep-sequencing mirna reads an map them to stem-loop sequences
# In: a fasta file containing the stem loops,
#        a fasta file containing the reads and their frequency
# Out: a fasta file containing the stem loop with counts for each base
#
# USAGE:
#   ./reads-to-stemloop.py -h
#
# OUTPUT
# # lines with metadata
# >hsa-mir-xxx
#   A   C   G   T   ...
#   xx  xx  xx  (counts of start sites at this position)
#   xx  xx  xx  (counts of end sites at this position)
#
# @author: Gregor Sturm (gregor.sturm@cs.tum.edu)
#
###

import sys
import warnings
import argparse
from Bio import SeqIO
import gzip
from nt_index import QuickNTIndex
from multiprocessing import Pool
from itertools import repeat

def normalizeFile(filename):
    sequences = []
    counts = []
    rpm = []
    try:
        gzReadfile = gzip.open(filename, 'r')
    except IOError, msg:
        warnings.warn(str(msg))

    for record in SeqIO.parse(gzReadfile, "fasta"):
        nReads = int(record.id.split("-")[1])             #id is of format >id-nReads
        sequences.append(str(record.seq))
        counts.append(nReads)

    baseSize = sum(counts)

    #update nReads to RPM
    for count in counts:
        rpm.append(count/float(baseSize) * 1e6)

    gzReadfile.close()
    return zip(sequences, rpm)

def getTotalSeqCount(index, bucket):
    """index is an arbitary position lying within the sequence"""
    end = sStem.find("$", index)
    start = sStem.rfind("$", 0, index)
    if end == -1:
        end = len(sStem)
    if start == -1:
        start = 0
    return sum(bucket[start:end])


def processReadfile((stemLen, filename)):
    multiple = 0        #reads that match to more than one sequence
    threshold = 0       #reads falling below the nReads threshold
    none = 0            #reads that don't match at all
    short = 0           #reads shorter than the index (not considered)
    total = 0           #reads processed
    matched = 0         #reads matching uniquely to a miRNA
    matches = 0         #reads matched, multiple matches counted multiple times
    startBucket = [0]*stemLen
    endBucket = [0]*stemLen

    multimatches = []

    sequences = normalizeFile(filename)
    for record in sequences:
        seq = record[0]
        rpm = record[1]
	if len(seq) < args.minlength:
	    short += 1
            continue
        if rpm < args.minreads:
            threshold += 1
            continue

        indexes = ntindex.all_indexes(seq)
        matches += len(indexes)
        if len(indexes) > 1:
            multimatches.append(record)
            continue
        if len(indexes) <= 0:
            none += 1
            continue
        matched += 1
	i = indexes[0]
        startBucket[i] += rpm
        endBucket[i + len(seq) -1] += rpm

    for record in multimatches:
        seq = record[0]
        rpm = record[1]
        indexes = ntindex.all_indexes(seq)
        seqCounts = {}
        for i in indexes:
            seqCounts[i] = getTotalSeqCount(i, startBucket)
        totalCounts = sum(seqCounts.values())
	if totalCounts == 0:
            #discard matches, we have no information about the relative expression
	    continue
        for i in indexes:
            ratio = seqCounts[i]/totalCounts
            startBucket[i] += ratio * rpm
            endBucket[i + len(seq) -1] += ratio * rpm

    total = len(sequences)
    multiple = len(multimatches)

    return short, none, threshold, multiple, total, matched, matches, startBucket, endBucket

if __name__ == "__main__":

## arguments
    argp = argparse.ArgumentParser(description="reads a file of deep-sequencing mirna reeds and maps them to the stem-loop sequences")
    argp.add_argument('-s', '--stemfile', type=file, action="store", dest="stemfile", required=True, help="file with mirna precursor sequences from mirbase (+30nts at each end)")
    argp.add_argument('-r', '--readfile', type=file, action="store", dest="readfile", required=True, help="file with list of files with sequencing experiments")
    argp.add_argument('-o', '--outfile', type=argparse.FileType('w'), action="store", dest="outfile", required=True)
    argp.add_argument('-t', '--minreads', type=int, action="store", dest="minreads", required=False, default=20)
    argp.add_argument('-l', '--minlength', type=int, action="store", dest="minlength", required=False, default=11)

    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))


## read stemloop file
    sStem =  ""
    for record in SeqIO.parse(args.stemfile, "fasta"):
        sStem += str(record.seq) + "$"

    ntindex = QuickNTIndex(sStem, args.minlength)

## process reads file
    p = Pool(16)
    files = [x.strip() for x in args.readfile.readlines()]
    results = p.map(processReadfile, zip(repeat(len(sStem)), files))
    shorts, nones, thresholds, multiples, totals, matcheds, matchess, startBuckets, endBuckets = zip(*results)

    short = sum(shorts)
    multiple = sum(multiples)
    threshold = sum(thresholds)
    none = sum(nones)
    total = sum(totals)
    matched = sum(matcheds)
    matches = sum(matchess)
    startBucket = [sum(i) for i in zip(*startBuckets)]
    endBucket = [sum(i) for i in zip(*endBuckets)]

## write output file
    args.outfile.write(" ".join(['#./reads-to-stemloop.py', "-r", args.readfile.name, '-s', args.stemfile.name, '-o', args.outfile.name, '-t', str(args.minreads), '-l', str(args.minlength)]) + "\n")
    args.outfile.write("#short\n")
    args.outfile.write("#" + str(short) + "\n")
    args.outfile.write("#threshold\n")
    args.outfile.write("#" + str(threshold) + "\n")
    args.outfile.write("#multiple\n")
    args.outfile.write("#" + str(multiple) + "\n")
    args.outfile.write("#none\n")
    args.outfile.write("#" + str(none) + "\n")
    args.outfile.write("#total\n")
    args.outfile.write("#" + str(total) + "\n")
    args.outfile.write("#matched\n")
    args.outfile.write("#" + str(matched) + "\n")
    args.outfile.write("#matches\n")
    args.outfile.write("#" + str(matches) + "\n")

    i = 0
    args.stemfile.seek(0)
    for record in SeqIO.parse(args.stemfile, "fasta"):
        args.outfile.write(">" + str(record.id) + "\n")
        seqlen = len(record.seq)
        args.outfile.write("\t".join(sStem[i:i+seqlen]) + "\n")
        args.outfile.write("\t".join(str(x) for x in startBucket[i:i+seqlen]) + "\n")
        args.outfile.write("\t".join(str(x) for x in endBucket[i:i+seqlen]) + "\n")
        i += seqlen + 1                                 #+1 for $

    args.stemfile.close()
    args.readfile.close()
    args.outfile.close()

