from lazy import lazy
import re
import numpy as np
from forgi_tools import ForgiTools
from rna_fold import RnaFold
from mk_mirna_objects import MirnaFactory

class HairpinException(Exception):
    pass

class Hairpin(object):
    """ Class that provides features for miRNA Hairpins.

    All Features are implemented in a lazy fashion, so that they
    are computed the first time they are requested. Then they
    are stored for later usage. All setter methods are required
    to reset all depending attributes to None.

    """

    ARGS = {}

    def __init__(self, name, seq, constraint=None,  **kwargs):
        """ init hairpin with name and sequence.

        Args:
            name: Name of the hairpin
            seq: nucleotide sequence of the hairpin
            constraint: Specify a dotbracket constraint for RNAfold (Optional)
        """
        self._name = name
        self._seq = seq.upper().replace("T", "U").strip()
        self._constraint = constraint

        self.args = {key: val for key, val in self.ARGS.iteritems()}
        for key in kwargs.iterkeys():
            assert self.args.get(key) is not None, "unknown keyword {0}".format(key)
        self.args.update(kwargs)

        assert re.search('[^ACGU]', self._seq) is None, \
                "input sequence has to be a valid DNA or RNA sequence"

    def __str__(self):
        return "\n".join([self._name, self.seq, self.dotbracket])

    @property
    def seq(self):
        return self._seq

    @property
    def length(self):
        return len(self._seq)

    @property
    def dna_seq(self):
        return self._seq.replace("U", "T")

    @property
    def name(self):
        return self._name

    @property
    def constraint(self):
        return self._constraint

    @lazy
    def dotbracket(self):
        """ secondary structure represented as ViennaRNA dotbracket sequence """
        dotbracket = self.rna_fold_tools.dotbracket
        return dotbracket

    @lazy
    def rna_fold_tools(self):
        """ return an instance of RnaFold and and run the folding """
        rnafold = RnaFold(self._seq, constraint=self._constraint)
        rnafold.fold()
        return rnafold

    @lazy
    def forgi_tools(self):
        """ return an instance of ForgiTools initialized with the dotbracket
        seq.
        """
        ft = ForgiTools()
        ft.from_dotbracket(self.dotbracket)
        return ft

    @lazy
    def bp_prop(self):
        """ list of overall basepair probabilities from RNAfold """
        bp_props = self.rna_fold_tools.bp_props
        total_bp_props = [sum(x.values()) for x in bp_props]
        return total_bp_props

    @lazy
    def shannon_entropy(self):
        """ get shannon entropy for the sequence """
        return self.rna_fold_tools.shannon_entropy

    @lazy
    def pairing_sequence(self):
        """ return a sequence representing the basepairing
        within a hairpin, such that:
            A: AT/TA
            C: CG/GC
            G: GT/TG
            U: unpaired
        """
        pmap = {"AU": "A", "CG": "C", "GU": "G", "unpaired": "U"}
        pairing_seq = []
        for i in range(len(self._seq)):
            partner = self.forgi_tools.forgi.pairing_partner(i
                    + self.forgi_tools.CORR)
            if partner is None:
               pairing_seq.append(pmap["unpaired"])
            else:
                nt1 = self._seq[i]
                nt2 = self._seq[partner - self.forgi_tools.CORR]
                letter = pmap.get(nt1 + nt2)
                if letter is None:
                    letter = pmap.get(nt2 + nt1)
                assert letter is not None
                pairing_seq.append(letter)
        return "".join(pairing_seq)

    @lazy
    def pairing_info(self):
        """ return a binary vector that indicates for every
        nucleotide, whether it is predicted to be paired (1)
        or unpaired (0)
        """
        return [0 if x == "U" else 1 for x in self.pairing_sequence]

    def _find_loop(self, forgi_tools, boundaries):
        candidates = forgi_tools.get_graph_annotation("h")
        if boundaries:
            candidates = filter(lambda c:
                    c[0] >= boundaries[0] and
                    c[1] <= boundaries[1], candidates)

        if len(candidates) == 0:
            raise HairpinException("no hairpin found")
        if len(candidates) > 1:
            raise HairpinException("multiple hairpins found")

        #the second position refers to the first position of the 3p stem
        return (candidates[0][0] - ForgiTools.CORR, candidates[0][1])

    def loop_coords(self, boundaries = None):
        """ position of the hairpin loop.

        Args:
            boundaries: Hairpin has to start and end within these
                two positions.

        Returns:
            Tuple (5p, 3p), where 5p referres to the first nt of the loop
            and 3p referres to the first nt of the 3p stem.
        """
        return self._find_loop(self.forgi_tools, boundaries)

    def get_3p_cs(self, cs5p):
        """ get the 3p processing site given a 5p processing site, based
        on the 2nt overhang """
        return self.forgi_tools.get_3p_overhang(cs5p, 2+1)

    def get_5p_cs(self, cs3p):
        """ the the 5p processing site given a 3p processing site, based
        ont the 2nt overhang """
        print self.dotbracket
        #+1 to refer to the first nt of the lower stem
        return self.forgi_tools.get_5p_overhang(cs3p, 2+1)

class AnnotatedHairpin(Hairpin):
    """ Extension of the Hairpin-Object, that can
    hold an annotation of the mature sequences and
    a list of keyword annotations.
    """

    def __init__(self, name, seq, cs5p=None, cs3p=None, constraint=None, **kwargs):
        """ Create an AnnotatedHairpin with mature sequences

        Args:
            @see Hairpin
            cs5p: 5' drosha cleavage site.
            cs3p: 3' drosha cleavage site.

            you have to provide at least one cleavage site.

        """
        super(AnnotatedHairpin, self).__init__(name, seq, constraint, **kwargs)

        self._annotation = []
        self._update_drosha_cs(cs5p, cs3p)


    def __str__(self):
        return "\n".join([
                    self.name,
                    self.seq,
                    self.dotbracket,
                    "CS: " + str(self.cleavage_site),
                    "annotation: " + ",".join(self.annotation)
            ])

    def add_annotation(self, keyword):
        """ add an annotation keyword """
        if keyword not in self._annotation:
            self._annotation.append(keyword)

    def remove_annotation(self, keyword):
        """ remove an annotation keyword. If keyword
        is not in annotation list, do nothing.
        """
        try:
            self._annotation.remove(keyword)
        except ValueError:
            pass

    @property
    def cleavage_site(self):
        """ get the cleavage site based on annotations and
        "guessing" based on the 2nt overhang
        """
        return self._cleavage_site

    @property
    def annotation(self):
        return self._annotation


    def check_integrity(self):
        """ check some basic properties that could indicate
        that something weng completely wrong
        """

        precursor_length = self.cleavage_site[1] - self.cleavage_site[0]
        assert 20 < precursor_length < 170,(
                "unrealistic precursor length %r" % precursor_length)

        assert 0 <= self.cleavage_site[0] <= self.length
        assert 0 <= self.cleavage_site[1] <= self.length
        assert self.cleavage_site[0] < self.cleavage_site[1]

    def _update_drosha_cs(self, cs5p, cs3p):
        """gets the annotated Drosha Cleavage site of a mirna-name.
        if only one strand is annotated, the other is determined
        by the 2nt 3p overhang.
        The position of the cleavage site refers to the first nucleotide of the
        hairpin and the first nucleotide of the lower 3' stem respectively
        (therefore the +1 on d3p)

        Returns:
            (5p, 3p), a tuple of the 5p and 3p Drosha cleavage site
        """

        assert cs5p is not None or cs3p is not None, "you have to provide at least one CS"

        if cs5p is not None and cs3p is not None:
            self.add_annotation("annotated_both")
        elif cs5p is not None:
            self.add_annotation("annotated_5p")
            cs3p = self.get_3p_cs(cs5p)
        elif cs3p is not None:
            self.add_annotation("annotated_3p")
            cs5p = self.get_5p_cs(cs3p)

        self._cleavage_site = (cs5p, cs3p)

class ExpressedHairpin(AnnotatedHairpin):
    """ determines the 'real' cleavage site based on expression data
    based on either mature sequences or offset reads
    """
    METHODS = ["annotated", "mor", "mature"]

    DEFAULT_ARGS = {
            "method" : "annotated",
            "prec_window" : 7,
            "strict" : True,
            "cs_offset_window" : 10
    }

    def __init__(self, annotated_hairpin, start_bucket, end_bucket, **kwargs):
        """ find the cleavage sites.

        Args:
            annotated_hairpin: preprocessed AnnotatedHairpin to which
                we want to add the expression data.
            start_bucket: list of the same length of seq, containing relative
                counts of reads starting at that position.
            end_bucket: ... counts of reads ending at that position
            method: how the cleavage site will be determined (based on
                annotation, mature sequences or offset reads (moRs)
            prec_window: size of the window that will be considered for
                precise/imprecise classification
            strict: boolean, default True. If True, the lack of expression
                for one of the methods will raise an Exception
            cs_offset_window: maximum shift of the cleavage site from
                the annotated canonical cleavage site.
            start_bucket_ht: Provide an alternative start_bucket with a
                higher threshold for mature sequences.

        Raises:
            HairpinException: when strict is True and the miRNA lacks
                in expression for mors or mature sequences.

        """

        self._name = annotated_hairpin.name
        self._seq = annotated_hairpin.seq
        self._constraint = annotated_hairpin._constraint
        self._annotation = annotated_hairpin._annotation

        #annotated (/guessed) CS
        self._cleavage_sites = {}
        self._cleavage_sites["annotated"] = annotated_hairpin.cleavage_site

        #new arguments
        self._start_bucket = map(float, start_bucket)
        self._end_bucket = map(float, end_bucket)
        self._start_bucket_ht = self._start_bucket

        if kwargs is not None and kwargs.get("start_bucket_ht") is not None:
            self._start_bucket_ht = map(float, kwargs["start_bucket_ht"])

        assert len(self._start_bucket) == len(self.seq)
        assert len(self._end_bucket) == len(self.seq)
        assert len(self._start_bucket_ht) == len(self.seq)

        self.args = self.DEFAULT_ARGS

        if kwargs is not None:
            for argname in self.args.keys():
                if kwargs.get(argname) is not None:
                    self.args[argname] = kwargs.get(argname)

        #determine CS with different methods
        self._find_cleavage_sites()

        self.method = self.args["method"]


    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method):
        """ set the 'method' to either annotated/mor/mature.

        The precision annotation is always based on the 5'
        strand. If 'annotated' is chosen, the precision
        will be based on the mature sequences
        """
        assert method in self.METHODS
        self._method = method
        self._cleavage_site = self._cleavage_sites[method]

        #add precision annotation
        cs5p = self.cleavage_site[0]
        prec_maker = HairpinPrecisionMeasures(
                self.get_cs_expression((cs5p-self.args["prec_window"],
                    cs5p+self.args["prec_window"]+1), "5p", method))

        #remove old annotation and add new one
        for prec_class in HairpinPrecisionMeasures.PREC_CLASSES:
            self.remove_annotation(prec_class)

        self.add_annotation(prec_maker.binary_classifier())

    @property
    def cleavage_sites(self):
        return self._cleavage_sites

    def _find_cleavage_sites(self):
        """ find the 'real' cleavage site in limited distance around
        the annotated one
        """
        annotated = self._cleavage_sites["annotated"]

        cs_mature5 = annotated[0] + self._get_cs_offset("5p", "mature")
        """ mature 3p based on guessing """
        cs_mature3 = self.get_3p_cs(cs_mature5)

        cs_mor5 = annotated[0] + self._get_cs_offset("5p", "mor")
        cs_mor3 = annotated[1] + self._get_cs_offset("3p", "mor")

        self._cleavage_sites["mor"] = (cs_mor5, cs_mor3)
        self._cleavage_sites["mature"] = (cs_mature5, cs_mature3)


    def _get_cs_offset(self, strand, method):
        """ get the offset to the annotated cleavage site
        based on different cleavage site detection methods
        """
        assert method in self.METHODS
        assert strand in ["3p", "5p"]
        if strand == "5p":
            cs = self._cleavage_sites["annotated"][0]
        elif strand == "3p":
            cs = self._cleavage_sites["annotated"][1]
        win = (cs - self.args["cs_offset_window"], cs + self.args["cs_offset_window"] + 1)
        exp = self.get_cs_expression(win, strand, method)
        assert len(exp) == win[1] - win[0]
        if sum(exp) == 0:
            if self.args["strict"]:
                raise Warning("%s %s cleavage site not expressed" % (method, strand))
            return 0
        else:
            self.add_annotation("exp_%s_%s" % (method, strand))
            return np.argmax(exp) - self.args["cs_offset_window"]


    def get_cs_expression(self, window, strand,  method=None):
        """ returns the corresponding expression values
        (either start or end bucket) based on the strand
        (5p, 3p) and the method (moR, mature).

        If method is set to annotated, the method will
        return the expression based on the mature
        sequence.

        window is a tuple with start and end position
        """

        if method is None: method=self.method
        assert method in self.METHODS
        assert strand in ["3p", "5p"]
        assert window[0] > 0 and window[1] <= len(self.seq)
        if method == "annotated": method = "mature" #get mature exp here

        s, e = window
        if method == "mature":
            if strand == "5p":
                return self._start_bucket_ht[s:e]
            elif strand == "3p":
                raise HairpinException(""" getting the 3p
                    strand based on  mature sequences
                    is inaccurate """)
        if method == "mor":
            if strand == "5p":
                return self._end_bucket[s-1 : e-1]
            if strand == "3p":
                return self._start_bucket[s:e]

    def check_consistency(self):
        """ add consistency annotations.

        Check if cleavage site and precision are consistent for
        the mor/mature/annotated method and add corresponding annotations.
        """

        #at the moment implemented in an ugly function in ipynb.
        assert False, "not implemented"

class HairpinPrecisionMeasures(object):
    """calculate different precision measures for windows
    around the cleavage site. The highest preak is per definition
    the cleavage site, therefore we don't have to consider
    the cleavage site position within the window
    """

    PREC_CLASSES = ["precise", "isomir", "unclassified"]

    def __init__(self, window):
        self.window = map(float, window)

    def precision1(self):
        """percentage of mirnas cleaved at *cleavage_site* in
        relation to the window of size PREC_WINDOW"""
        max_peak = max(self.window)
        window_sum = sum(self.window)
        if window_sum == 0:
            return -1
        else:
            return float(max_peak)/window_sum

    def binary_classifier(self):
        """two-class-classifier: every miRNA with one single peak
        higher than 1e2 than all others is declared as precise.
        Everything within 1e1 as unprecise. The values between
        remain unclassified to have a higher difference
        """
        prec_class = "unclassified"
        tmp_window = self.window[:]
        peak1 = self._pop_peak(tmp_window)
        peak2 = self._pop_peak(tmp_window)
        if(peak1 > 1e2 * peak2):
            prec_class = "precise"
        elif(peak1 < 1e1 * peak2):
            prec_class = "isomir"
        return prec_class

    def _pop_peak(self, window):
        """returns the current highest peak  in window
        and removes it"""
        i = np.argmax(window)
        peak = window[i]
        window[i] = 0
        return peak

    def discrete_classifier(self):
        prec_class = "unclassified"
        tmp_window = self.window[:]
        if sum(tmp_window) < 1e4:
            return prec_class
        peak1 = max(tmp_window)

        #search for "big" peaks
        big = 0
        for v in tmp_window:
            if v >= 1e-1 * peak1:
                big += 1

        #search for "small" peaks:
        small = 0
        for v in tmp_window:
            if v < 1e-1 * peak1 and v >= 1e1 and v >= 1e-5 * peak1:
                small += 1

        appendix = "_small" if small else ""
        prec_class = "".join(["peak", str(big), appendix])
        return prec_class

class CandidateHairpin(Hairpin):
    """ Extention of the Hairpin class to incorporate Features
    that depend on the cleavage site.
    """

    ARGS = {
        "n_subopt"       : 500,  #number of suboptimal structures to draw
        "additional_nts" : 20,   # #nucleotides beyond the CS
        "hairpin_offset" : 10           # distance from CS + this value
    }

    def __init__(self, name, seq, candidate_sites, constraint=None, **kwargs):
        """

        Args:
            seq: the precursor sequence
            candidate_sites: tuple with position (5p CS, 3p CS)
        """
        super(CandidateHairpin, self).__init__(name, seq, constraint)

        # this is for testing only
        self.original_candidate_site = candidate_sites[0]

        self.args = {key: val for key, val in self.ARGS.iteritems()}
        for key in kwargs.iterkeys():
            assert self.args.get(key) is not None, "unknown keyword {0}".format(key)
        self.args.update(kwargs)

        #require enough nts beyond processing site
        assert (candidate_sites[0] >= self.args["additional_nts"]
                and candidate_sites[1] <= len(self._seq) - self.args["additional_nts"]),\
                "sequence too short"

        #trim sequence and save candidate site
        add_nts = self.args["additional_nts"]
        self._seq = self._seq[ candidate_sites[0] - add_nts :
                        candidate_sites[1] + add_nts + 1]
        if constraint is not None:
            #balance constraint when shortened.
            constraint = MirnaFactory._remove_constraints(self._constraint, candidate_sites)
            self._constraint = self._constraint[ candidate_sites[0] - add_nts :
                            candidate_sites[1] + add_nts + 1]

        self.candidate_sites = (add_nts, len(self._seq) - add_nts)

        self.loop_boundaries = (self.candidate_sites[0] + self.args["hairpin_offset"],
                self.candidate_sites[1] - self.args["hairpin_offset"])

    @lazy
    def loop_coords(self):
        """@override"""
        return super(CandidateHairpin, self).loop_coords(self.loop_boundaries)

    @lazy
    def loop_distance(self):
        """ distance from the 5' cleavage site to the 5' start
        of the hairpin loop
        """
        return self.loop_coords[0] - self.candidate_sites[0]

    def _analyse_subopt_structures(self):
        if hasattr(self, "_junction_distrib"): return
        print " analysing subopt structures for {0}".format(self.name)
        dotbrackets = RnaFold.rnasubopt(self._seq, self.args["n_subopt"])
        junction_distrib = [0] * (self.candidate_sites[0] + 1) #distribution with bins
        hairpin_dists = [] #just make list here
        loop_dists = []
        ft = ForgiTools()
        for db in dotbrackets:
            ft.from_dotbracket(db)
            try:
                loop_coords = self._find_loop(ft, self.loop_boundaries)
                htd = ft.modified_distance(self.candidate_sites[0], loop_coords[1])
                hairpin_dists.append(htd)
                loop_dist = loop_coords[0] - self.candidate_sites[0]
                loop_dists.append(loop_dist)
            except HairpinException:
                pass
            junct_dist = ft.get_junction_distance(self.candidate_sites[0])
            junction_distrib[junct_dist] += 1

        self._junction_distrib = junction_distrib
        self._loop_dists = loop_dists
        self._hairpin_dists = hairpin_dists


    @lazy
    def hairpin_tip_distribution(self):
        self._analyse_subopt_structures()
        return self._hairpin_dists

    @lazy
    def loop_dist_distribution(self):
        """ distribution of the distances from the cand. site to the loop-junction
        """
        self._analyse_subopt_structures()
        return self._loop_dists

    @lazy
    def hairpin_tip_distance(self):
        """ distance that estimates the molecular distance between the cleavage site
        and the hairpin tip for the most likely secondary structure"""
        ft = self.forgi_tools
        htd = ft.modified_distance(self.candidate_sites[0], self.loop_coords[1])
        return htd

    @lazy
    def junction_distribution(self):
        """ get predicted abundances for the 5' stem-junction.

        Run RNAsubopt and draw sequences according to their energy-level.
        Run Forgi to find the stem-junction of each sample and make a distribution
        """
        self._analyse_subopt_structures()
        return self._junction_distrib

    @lazy
    def junction_distance(self):
        """ get the distance from the cleavage site to the 5'
        stem junction.

        compute the junction_distribution and take the most abundant
        junction while ignoring the "cs" position (all predictions that
        end up in not having a junction).
        """
        print "evaluating jd {0}".format(self._name)

        distrib = self.junction_distribution
        dist = np.argmax(distrib[:-1])
        return dist

    @lazy
    def precursor_length(self):
        """ length of the pre-miRNA == distance between
        the two cleavage sites
        """
        return self.candidate_sites[1] - self.candidate_sites[0]

