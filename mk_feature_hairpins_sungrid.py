#!/usr/bin/env python
# -*- coding: utf-8 -*-

__descr__ = """ reads Annotated/Expressed Hairpin objects from a (pickled) dataset
and submits one sungrid job each to create candidate hairpins for feature analysis
"""

import argparse
from config import TMP
from subprocess import call
import pickle
import os

def main(dataset, subfolder):
    dir = TMP+subfolder+"/"
    call(["mkdir", "-p", dir])
    all_mirs = pickle.load(dataset).itervalues()
    for mirna in all_mirs:
        name  = mirna.name
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        command = ["qsub",  "-b", "y",
            "-N", mirna.name,
            "-P", "short_proj",
            "-cwd",
            "-e", "{0}{1}.err".format(dir, name),
            "-o", "{0}{1}.out".format(dir, name),
            "/home/st/sturmg/bin/anaconda/bin/python {0}/mk_feature_hairpins.py --name {1} --sequence {2} --cs5p {3} --cs3p {4} --constraint \"{6}\" --output {5}{1}.pickle".format(scriptdir,
                name, mirna.seq, mirna.cleavage_site[0], mirna.cleavage_site[1], dir, mirna._constraint)]
        print " ".join(command)
        call(command)



if __name__ == "__main__":
    argp = argparse.ArgumentParser(description=__descr__)
    argp.add_argument("-d", "--pickled_dataset", type=file,
            action="store", dest="dataset", required=True,
            help="Path to pickled objects")
    argp.add_argument("-s", "--subfolder", type=str,
            action="store", dest="subfolder", required=True,
            help="subfolder in tmp")
    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))

    main(args.dataset, args.subfolder)
