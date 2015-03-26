#!/usr/bin/env python
# -*- coding: utf-8 -*-

__descr__ = """ script retrieves a name, seq and cleavage site of a hairpin and
computes the corresponding CandidateHairpin Object """

import argparse
import pickle
from timer import Timer
from hairpin import CandidateHairpin

def main(name, seq, cs, constraint, outfile):
    hc = None
    with Timer("mk_candidate"):
        hc = CandidateHairpin(name, seq, cs, constraint, n_subopt=500)
        print hc.hairpin_tip_distance
        print hc.hairpin_tip_distribution[:5]
        print hc.loop_dist_distribution[:5]
        print hc.junction_distribution
        print hc.junction_distance
        print hc.loop_coords

    with Timer("pickleing"):
        pickle.dump(hc, outfile)

if __name__ == "__main__":
    argp = argparse.ArgumentParser(__descr__)
    argp.add_argument('-n', '--name', type=str,
            action="store", dest="name", required=True,
            help="Unique name of the Job")
    argp.add_argument('-s', '--sequence', type=str,
            action="store", dest="seq", required=True,
            help="RNA sequence, centered on the hairpin loop")
    argp.add_argument('-f', '--cs5p', type=int,
            action="store", dest="cs5p", required=True)
    argp.add_argument('-t', "--cs3p", type=int,
            action="store", dest="cs3p", required=True)
    argp.add_argument('-c', "--constraint", type=str,
            action="store", dest="constraint", required=False,
            default=None, help="add additional constraint for RNAfold")
    argp.add_argument('-o', '--output', type=argparse.FileType('wb'),
            action="store", dest="outfile", required=True,
            help="Output file. Will be a pickled CandidateHairpin")


    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))

    main(args.name, args.seq, (args.cs5p, args.cs3p), args.constraint, args.outfile)


