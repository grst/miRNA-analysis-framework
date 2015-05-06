import forgi.graph.bulge_graph as fgb
import numpy as np


class ForgiTools(object):

    CORR = 1        #constant to remember that forgi starts counting at 1

    def __init__(self):
        self._forgi = None

    def from_dotbracket(self, dotbracket):
        self._forgi = fgb.BulgeGraph()
        self._forgi.from_dotbracket(dotbracket)

    def from_forgi(self, forgi):
        self._forgi = forgi

    @property
    def forgi(self):
        return self._forgi

    def get_graph_annotation(self, key = None):
        """parse forgi graph annotation according to
        http://www.tbi.univie.ac.at/~pkerp/forgi/graph_tutorial.html#loading-a-structure-from-a-bpseq-formatted-file

        Args:
            key: the fist letter after "define". e.g. h = hairpin

        Returns:
            list of all lines with the keyword "define" and the key "key".
            Every line is split into a list containing the int-coordinates.
        """
        assert self._forgi is not None

        bg_result = self._forgi.to_bg_string().split("\n")
        bg_result = [line.split(" ") for line in bg_result]
        bg_result = filter(lambda line: line[0] == "define", bg_result)
        if(key != None):
            bg_result = filter(lambda line: line[1][0] == key, bg_result)
        coordinates = [[int(x) for x in line[2:]] for line in bg_result]
        return coordinates

    def modified_distance(self, start, end, custom_factors = None):
        """ heuristic for a molecular distance between two nucleotides
        in a hairpin. Takes bulges and loops into account.
        The given factors are fitted to minimize the relative standard
        deviation (RSD).
        """
        factors = {
                "f": lambda x: x,
                "t": lambda x: x,
                "s": lambda x: x,
                "i": lambda x: x * (4/11.),
                "h": lambda x: x/2.,
                "m": lambda x: x
        }
        if custom_factors is not None:
            factors.update(custom_factors)

        elems = [['s', 0]]      #init list with first element (type, length)

        step = -1 if end < start else 1
        for i in range(start, end, step):
            struct = self._forgi.get_node_from_residue_num(i + self.CORR)
            struct = struct[0]  #first character e.g. h0 -> h
            if elems[-1][0] == struct:
                elems[-1][1] += 1
            else:
                elems.append([struct, 1])

        dist = 0
        for struct, num in elems:
            dist += factors[struct](num)

        return dist

    def get_5p_overhang(self, pos3p, overhang):
        """ Get the corresponding molecule on the 5p strand, given a nt on the
        3p strand.

        Args:
            pos3p: the position on the 3p strand
            overhang: e.g. Drosha cs has a +2nt overhang.

        """

        i = 0
        d5p = 0
        while True:
            d5p = self._forgi.pairing_partner(pos3p + self.CORR - i - overhang)
            print d5p
            if d5p != None:
                break
            i += 1
        #if i:
        #    self.add_annotation("bulge_3p")
        d5p -= i
        d5p -= self.CORR
        return d5p

    def get_3p_overhang(self, pos5p, overhang):
        """ Get the corresponding molecule on the 3p strand, given a nt on the
        5p strand.

        See get_5p_overhang.
        """

        i = 0
        d3p = 0
        while True:
            d3p = self._forgi.pairing_partner(pos5p + self.CORR + i)
            if d3p != None:
                break
            i += 1
        #if i:
        #    self.add_annotation("bulge_5p")
        d3p += i + overhang
        d3p -= self.CORR
        return d3p

    def get_junction_distance(self, cs):
        """ find the stem junction in dotbracket backwards from cs.

        Returns the distance to from the cs to the first nt
        of the dangling ends.
        """
        loop = False

        c = cs
        d5p = -1
        for i in range(cs, 0, -1):
            struct = self._forgi.get_node_from_residue_num(i + self.CORR)
            if struct[0] == 'f':
                d5p = c
                break
            elif struct[0] == 'm':
                if loop:
                    d5p = c
                    break
                else:
                    loop = True
            else:
                loop = False
                c-=1
        if d5p < 0: d5p = 0

        return cs - d5p
