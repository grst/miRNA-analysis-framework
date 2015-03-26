#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################
# parser for mirbase stemloop annotations.
#
# USAGE:
#   parse_mirbase_stemloops2.py miRNAs.str
#
# was given to me by kristin (kristin.wahl@stud.ntnu.no),
# but it's not written by her. She said that it didn't
# work properly, but I just had to adjust some input/
# output issues for my purposes.
#
# original author: ?
# adapted on 2015-01-15: Gregor Sturm (gregor.sturm@cs.tum.edu)
#
#################################################################

import sys


def main(secondary_file):
    """
    Read in secondary info from local fasta file,
    converting it to 'parenthesis style'.
    Return a map of: 'miRNA_id -> secondary structure'.
    The structure is given in 'parenthesis style', see below
    """

    def rnafoldStyle(data):
        """
        Convert from mirbase secondary structure format:
                 c   -  u            c  c    gagc
        ccuauguag ggc ca caaaguggaggc cu ucuu    c
        ||||||||| ||| || |||||||||||| || ||||
        ggguacguc ccg gu guuucaccuucg ga agag    u
                 a   u  -            u  a    uaag
        to:
        (((((((((.(((((.((((((((((((.((.((((..........)))).)).)))))))))))))).))).)))))))))
        """
        res = ""
        nucs = "acgu"
        for top, bot in zip(data[0], data[1][:-1]):
            if top in nucs:
                res += "."
            elif bot in nucs:
                res += "("
        if data[1][-1] in nucs:
            res += "."
        if data[2][-1] in nucs:
            res += "."
        try:
            if data[3][-1] in nucs:
                res += "."
        except IndexError:
            print data
        for bot, top in reversed(zip(data[4][:-1], data[3])):
            if bot in nucs:
                res += "."
            elif top in nucs:
                res += ")"
        return res

    f = open(secondary_file)
    lines = [l[:-1].lower() for l in f.readlines()]
    structMap = {}
    for i in range(0, len(lines), 8):
        mid = lines[i].split(" ")[0]
        structMap[mid] = rnafoldStyle(lines[i + 2: i + 7])
    f.close()
    return structMap


if __name__ == '__main__':
    structMap = main(sys.argv[1])
    for key, dotbracket in sorted(structMap.items()):
       print key
       print dotbracket

