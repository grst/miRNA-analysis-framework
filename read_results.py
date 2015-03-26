##################################################
# Functions for reading the outputs of the
# ba-thesis scripts
#
# author: Gregor Sturm (gregor.sturm@cs.tum.edu)
##################################################

from Bio import SeqIO
import collections

def read_precursor_results(inputfile):
    """read the outputs of reads-to-stemloop.py, returns a list with records"""

    opened = False
    if not hasattr(inputfile, 'read'):
        opened = True
        inputfile = open(inputfile, 'r')

    results = []
    i = 0
    for line in inputfile.readlines():
        line = line.strip()
        if line[0] == "#":
            continue
        if i%4 == 0:
            tmpName = line[1:].lower()
        if i%4 == 1:
            tmpSeq = line.split("\t")
        if i%4 == 2:
            tmpStart = line.split("\t")
        if i%4 == 3:
            tmpEnd = line.split("\t")
            results.append({"seq": tmpSeq, "start": tmpStart, "end": tmpEnd, "name": tmpName})
        i += 1

    if opened:
        inputfile.close()

    return results

def read_to_dict(inputfile):
    """reads the file line by line into a dictionary such that
    it can be checked in constant time, whether a specific
    line is in the file.

    dict.get(line) == None

    Args:
        inputfile: can either be a file or a pathname."""

    opened = False
    if not hasattr(inputfile, 'read'):
        opened = True
        inputfile = open(inputfile, 'r')

    dict = {}
    for line in inputfile:
        line = line.strip().lower()
        dict[line] = True

    if opened:
        inputfile.close()

    return dict

def read_fasta_to_dict(inputfile):
    """reads a fasta file into a dictionary, such that
    dict[header] = sequence

    Args:
        inputfile: can either be a file or a pathname"""

    dict = collections.OrderedDict()
    for record in SeqIO.parse(inputfile, "fasta"):
        dict[str(record.id).lower()] = str(record.seq)

    return dict



