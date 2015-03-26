#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script retrieves a RNA sequence of a miRNA Hairpin and
generates probable candidates for Microprocessor cleavage
sites """

import argparse
from timer import Timer
import pickle

import multiprocessing as mp
from hairpin import Hairpin, CandidateHairpin, HairpinException

class NoCandidateException(Exception):
    pass

def get_junction_distance(obj):
    return obj.junction_distance

class MirnaController(object):

    ARGS = {
        "hairpin_tip_distance" : 25,
        "min_junction_distance" : 13,
        "candidate_window" : 4,
        #minimal distance of the hairpin from the seqence start
        "hairpin_offset" : 40,
        #nucleotides beyond the cleavage site in the CandidateHairpins
        "additional_nts" : 20,
        "n_subopt": 10,
        "threads" : 3
    }

    def __init__(self, name, seq, **kwargs):
        self.args = {key: val for key, val in self.ARGS.iteritems()}
        for key in kwargs.iterkeys():
            assert self.args.get(key) is not None, "unknown keyword {0}".format(key)
        self.args.update(kwargs)

        self.name = name
        self.seq = seq
        self.candidates = []

    def mk_candidates(self):
        """ find the hairpin loop and generate processing site candidates
        """
        print "processing candidates for {0}".format(self.name)

        input_hairpin = Hairpin(self.name, self.seq)
        boundaries = [self.args["hairpin_offset"],
                len(self.seq) - self.args["hairpin_offset"]]
        hp5p, hp3p = input_hairpin.loop_coords(boundaries)

        print "hairpin detected at {0},{1}".format(hp5p, hp3p)

        for i in range(self.args["additional_nts"], hp5p):
            cs5p = i
            cs3p = input_hairpin.get_3p_cs(cs5p)
            #check if additional nts are also available on 3p site!
            if len(self.seq) - cs3p < 20:
                continue
            candidate = CandidateHairpin("{0}_{1}".format(self.name ,  i),
                    self.seq, (cs5p, cs3p),
                    additional_nts=self.args["additional_nts"],
                    n_subopt=self.args["n_subopt"])
            self.candidates.append(candidate)

        if len(self.candidates) == 0:
            raise NoCandidateException()

        print "{0} candidates found".format(len(self.candidates))

    def filtering_hairpin_distance(self):
        """ filter candidates according to the distance
        from the processing site to the hairpin tip
        """
        print "started filtering for hairpin tip distance"

        boundaries = (self.args["hairpin_tip_distance"]
                - self.args["candidate_window"], self.args["hairpin_tip_distance"]
                + self.args["candidate_window"])

        filtered_candidates = []
        for candidate in self.candidates:
            print "candidate {0}".format(candidate.name)
            try:
                if boundaries[0] <= candidate.hairpin_tip_distance <= boundaries[1]:
                   filtered_candidates.append(candidate)
                   print "valid candidate: {0}".format(candidate.name)
            except HairpinException, msg:
                pass

        self.candidates = filtered_candidates

        if len(self.candidates) == 0:
            raise NoCandidateException("after filtering 1")

        print "{0} candidates found after filter1".format(len(self.candidates))

    def mk_subopt_features(self):
        """ calculate all features that are based on suboptimal structures.

        This process relies on RNAsubopt and yield many
        different secondary structures. This is the most
        computationally intensive step and therefore
        has to be ran in parralel.
        """

#        pool = mp.Pool(self.args["threads"])
#        pool.map(get_junction_distance,
#                self.candidates)

        #for cand in self.candidates:
        #    f = methodcaller('junction_distance')
            #print f(cand)

        print "started filtering for stem-junction-distance"
        filtered_candidates = []
        for cand in self.candidates:
            print cand.hairpin_tip_distance
            print cand.hairpin_tip_distribution[:5]
            print cand.loop_dist_distribution[:5]
            print cand.junction_distribution
            print cand.junction_distance
            print cand.loop_coords


        self.candidates = filtered_candidates




def main(name, seq, outfile):
    mc = MirnaController(name, seq, threads=16, n_subopt=500, candidate_window=7, min_junction_distance=0)
    with Timer("mk_candidates"):
        mc.mk_candidates()

    with Timer("filter hairpin distance"):
        mc.filtering_hairpin_distance()

    with Timer("filter stem junction"):
        mc.mk_subopt_features()

    pickle.dump(mc, outfile)

if __name__ == "__main__":
    argp = argparse.ArgumentParser(description="generate probable candidates for Microprocessor Cleavage")
    argp.add_argument('-n', '--name', type=str,
            action="store", dest="name", required=True,
            help="Unique name of the Job")
    argp.add_argument('-s', '--sequence', type=str,
            action="store", dest="seq", required=True,
            help="RNA sequence, centered on the hairpin loop")
    argp.add_argument('-o', '--output', type=argparse.FileType('wb'),
            action="store", dest="outfile", required=True,
            help="Output file. Will be a pickled MirnaController")

    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))

    main(args.name, args.seq, args.outfile)
