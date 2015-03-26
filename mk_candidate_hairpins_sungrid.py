#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from read_results import *
import os
from subprocess import call
from config import TMP

def main(fasta, subfolder):
    dir = TMP + subfolder + "/"
    call(["mkdir", "-p", dir])
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    seqs = read_fasta_to_dict(fasta)
    for name, seq in seqs.iteritems():
        command = ["qsub",  "-b", "y",
            "-N", "{1}",
            "-P", "short_proj",
            "-cwd",
            "-e", "{0}{1}.err".format(dir, name),
            "-o", "{0}{1}.out".format(dir, name),
            "/home/st/sturmg/bin/anaconda/bin/python {0}/mk_candidates.py --name {1} --sequence {2} --output {3}{4}.pickle".format(scriptdir, name, seq, dir, name)]
        print " ".join(command)
        call(command)


if __name__ == "__main__":
    argp = argparse.ArgumentParser(description="read mirna sequences from a fasta file and submit them on the grid")
    argp.add_argument("-s", "--sequences", type=file,
            action="store", dest="fasta", required=True,
            help="Path to input-fasta file")
    argp.add_argument("-d", "--subfolder", type=str,
            action="store", dest="subfolder", required=True,
            help="subfolder in tmp")

    try:
        args = argp.parse_args()
    except IOError, msg:
        argp.error(str(msg))

    main(args.fasta, args.subfolder)
