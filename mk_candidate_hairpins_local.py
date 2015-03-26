#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from read_results import *
import os
from subprocess import call
from multiprocessing import Pool
from config import TMP, PYTHON

def run(command):
    call(command)

def main(fasta, subfolder):
    pool = Pool(4)
    dir = TMP + subfolder + "/"
    call(["mkdir", "-p", dir])
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    seqs = read_fasta_to_dict(fasta)
    commands = []
    for name, seq in seqs.iteritems():
        command = [PYTHON,
                "{0}/mk_candidate_hairpins.py".format(scriptdir),
                "--name", name,
                "--sequence", seq,
                "--output", "{0}{1}.pickle".format(dir, name)]
        print " ".join(command)
        commands.append(command)

    pool.map(run, commands)

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
