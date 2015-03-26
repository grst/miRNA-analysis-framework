import subprocess
import hashlib
import os
import numpy as np
from lazy import lazy
from config import *

class RnaFold(object):

    CORR = 1 #Reminder that Vienna starts indexing at 1

    def __init__(self, sequence, constraint=None, gamma=1.):
        self.initialized = False
        self.seq = sequence
        self.constraint = constraint
        self.gamma = gamma
        m = hashlib.md5()
        m.update(sequence + str(constraint) + str(gamma))
        self.name = m.hexdigest()

    def fold(self):
        """ run RNAfold and generate the dotplot files """
        command = [ RNAFOLD_BIN + "RNAfold",
                    "-p",
                    "--MEA={0}".format(self.gamma),
                    "--bppmThreshold=1e-6"]
        inputfile = [">"+self.name, self.seq]
        if self.constraint is not None:
            command.append("-C")
            inputfile.append(self.constraint)
        inputfile = "\n".join(inputfile)
        echo = subprocess.Popen(["echo", inputfile], stdout=subprocess.PIPE)
        echo.wait()
        proc = subprocess.Popen(command, stdin=echo.stdout, stdout=subprocess.PIPE)
        proc.wait()

        stdout, stderr = proc.communicate()
        assert stderr is None
        self._dotbracket = stdout.split("\n")[5].split()[0]

        self._mk_bp_propensities()
        self._mk_shannon_entropy()
        os.remove("{0}_ss.ps".format(self.name))
        os.remove("{0}_dp.ps".format(self.name))
        self.initialized = True

    def _mk_bp_propensities(self):
        """ extract the bp probatilities from the RNAfold dotplot
        postscript file"""
        bp_props = {}
        with open("{0}_dp.ps".format(self.name)) as tmp_plot:
            for line in tmp_plot.readlines():
                line.strip()
                line = line.split()
                if(len(line) == 4 and line[3] == "ubox"):
                    id1 = int(line[0])
                    id2 = int(line[1])
                    prop = float(line[2]) ** 2  #square to proability
                    assert 0 <= prop <= 1
                    bp_list1 = bp_props.get(id1)
                    bp_list2 = bp_props.get(id2)
                    if bp_list1 == None: bp_list1 = {}
                    if bp_list2 == None: bp_list2 = {}
                    bp_list1[id2] = prop
                    bp_list2[id1] = prop
                    assert 0 <= sum(bp_list1.values()) <= 1
                    assert 0 <= sum(bp_list2.values()) <= 1
                    bp_props[id1] = bp_list1
                    bp_props[id2] = bp_list2

        self._bp_props = bp_props

    def _mk_shannon_entropy(self):
        """ calculate the positional shannon entropy as described in
        Huynen, Gutell and Konings (1997)"""
        entropy = {}
        for nt_id, bp_list in self._bp_props.iteritems():
            entropy[nt_id] = sum([x * np.log(x) for x in bp_list.values() if x > 0])
            self_non_pairing = 1 - sum(self._bp_props[nt_id])
            if self_non_pairing > 0:
                entropy[nt_id] += self_non_pairing * np.log(self_non_pairing)
            entropy[nt_id] *= -1

        self._shannon_entropy = entropy


    @property
    def dotbracket(self):
        assert self.initialized
        return self._dotbracket

    @lazy
    def bp_props(self):
        """ list with dict of the bp props with other bases"""
        assert self.initialized
        # turn dict into ordered list
        bp_props = []
        for i in range(len(self.seq)):
            i += self.CORR
            bps = self._bp_props.get(i)
            if bps is None:
                bp_props.append({})
            else:
                bp_props.append(bps)

        return bp_props

    @lazy
    def shannon_entropy(self):
        assert self.initialized
        # turn dict into list
        shannon_entropy = []
        for i in range(len(self.seq)):
            i += self.CORR
            entropy = self._shannon_entropy.get(i)
            shannon_entropy.append(entropy if entropy is not None else 0)
        return shannon_entropy


    @staticmethod
    def rnasubopt(seq, n):
        """ Run rnasubopt on seq.

        Run RNAsubopt with the -p option and draw *n* sequences.
        """
        command = [RNAFOLD_BIN + "RNAsubopt", "-p{0}".format(n), "--noClosingGU"]
        echo = subprocess.Popen(["echo", seq], stdout=subprocess.PIPE)
        echo.wait()
        proc = subprocess.Popen(command, stdin=echo.stdout, stdout=subprocess.PIPE)
        proc.wait()

        stdout, stderr = proc.communicate()
        assert stderr is None
        lines = stdout.strip().split("\n")
        """ second column contains an energy """
        structures = [line.split()[0] for line in lines[1:]]
        return structures

    @staticmethod
    def quickfold(seq):
        """ just fold and return dotbracket """
        command = [RNAFOLD_BIN + "RNAfold"]
        seqs = "\n".join([">{0}".format("rnafold"), seq])
        echo = subprocess.Popen(["echo", seqs], stdout=subprocess.PIPE)
        echo.wait()
        proc = subprocess.Popen(command, stdin=echo.stdout, stdout=subprocess.PIPE)
        proc.wait()

        stdout, stderr = proc.communicate()
        assert stderr is None
        lines = stdout.strip().split("\n")
        assert len(lines) == 3

        """dotbracket is the third line of the output and contains
        an energy in the second column"""
        dotbracket = lines[2].split()[0]
        assert len(dotbracket) == len(seq)
        return dotbracket


