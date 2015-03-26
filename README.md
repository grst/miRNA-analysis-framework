# miRNA Analysis framework
This module provides a toolbox for generatinga an analysing Features
of pri-miRNA hairpins. 

# Install
## Compatibility
The package was developed and tested on Linux Mint 17.1 with Anaconda Python 2.7.9

##Prerequisites
This package relies on a series of open source libraries and programs.
Python dependancies can easily be downloaded

`pip install lazy`
`pip install forgi`
`pip install biopython`
`pip install bcbio-gff`

Download and install the ViennaRNA package from 
http://www.tbi.univie.ac.at/RNA/.

If not installed globally, adapt the path to the RNAFold binary in `config.py` 

Further you might need: 
* numpy
and for plotting: 
* matplotlib
* the *seaborn* addon

# Scripts
* `parse_mirbase_structures.py` parses structural information from mirbase and converts it to the dotbracket format
* `mk_mirna_objects.py` reads mirbase annotations and generates Hairpin-Objects. See below in the section "Build miRNA objects"
* `reads_to_stemloop.py` maps reads from sequencing experiments to the sequences of the hairpins. 
* `mk_candidate_hairpins.py` generates CandidateHairpins for a Hairpin that can be used for MachineLearning.
* `mk_candidate_haipins_sungrid.py` submits multiple `mk_candidate_hairpins.py` instances on a qsub system. 
* `mk_feature_hairpins.py` generates CandidateHairpins with a given cleavage site, for feature Analysis
* `mk_feature_hairpins_sungrid.py` submides multiple instances of `mk_feature_hairpins` on a qsub system. 

# Usage

## Get Data
* miRBase database files from ftp://mirbase.org/pub/mirbase/CURRENT/
* genomes 
    Get the corresponding genome (hg38 in my case)

## Build miRNA objects
With this script one can read the fasta-sequences and annotations from mirbase and create Hairpin-Object. The Hairpin-Objects will be extended by mapping the sequences to the genome, such that they have an equal amount of nucleotides beyond the cleavage site. 
`mk_mirna_objects.py`

Args:
    precursor_file: hairpin.fa from mirbase
    constraint_file: additional constraints in fasta/dotbracket format. The constraints
        from mirbase (miRNA.str.gz) can be converted to dotbracket
        with the `parse_mirbase_structures` script. Other than indicated, this is not optional. You can use a file 
        with empty constraints ("." * len(seq)) if you don't want to use the constraints from mirbase. 
    mature_file: mature.fa from mirbase
    coordinate_file: genome.gff from mirbase
    genome_path: path to the full genome, split up in
        chromosomes, s.t. there are files named
        chr1.txt .. chrY.txt, only containing the sequences
        in plain text. You can use fastaget from the
        biosnippets repository (https://github.com/grst/biosnippets)
        to make such files from a fasta-genome. 
    pickle_file: output file, here it will pickle a file containing 
        AnnotatedHairpin Objects
    fasta_file: the same miRNAS, only as fasta sequences. 
    annotation_files: A list of files, containing one mirna-id per line,
        the miRNA will be annotated with the label specified
        in annotation_labels.
    annotation_labels: A list of labels corresponding to the annotation
        files. 

Example:
    `ipython ./mk_mirna_objects.py -- -p ../../data/mirbase21/hsa-mir.fa --constraint ../../data/mirbase21/miRNA.dotbracket.str --mature ../../data/mirbase21/hsa-mature.fa --gff ../../data/mirbase21/hsa.gff3 --genome-path ../../data/genomes/ -e 30 --annotation ../../data/mirbase21/filters/mirtrons.txt ../../data/mirbase21/selections/conserved-hsa.txt ../../data/mirbase21/selections/hsa-confidence-high.txt --annotation_labels mirtron conserved high_confidence --output_pickle ../../results/mirna-objects/hsa-30-2.pickle --output_fasta ../../results/mirna-objects/hsa-30-2.fa`

Further Information:
    Error logs will be written to STDOUT. MiRNAs that fail to be extended (e.g. due to bad
    RNAfold folding will still be included but get the annotaion "extension_fail". 


## Map sequencing data to precursor hairpin sequences
Use `reads-to-stemloop.py` to map reads from deep-sequencing to the extended precursor sequences generated with `mk_mirna_objects.py`

## Generate ExpressedHairpins
Use the data from `reads-to-stemloop.py` (read the data with read_results.read_precursor_results) to generate ExpressedHairpins. 
This class provides cleavage site detection based on different methods (mature/offset reads) and can be used to evaulate the consistency of the data. 
Most of the features (like base-pair probability or entropy) can directly be retrieved from this class. 

## Generate CandidateHairpins
Generating a CandidateHairpin is computationally intensive, since it involves the analysis of many suboptimal secondary structures. The script `mk_candidate_hairpins.py` generates an instance for all possible cleavage site positions of a miRNA hairpin. These instances can be used to build a training/test set for machine learning. The script `mk_feature_hairpins.py` generates only one instance at the cleavage site as specified by the ExpressedHairpin. The resulting instance contains
all possible features and can be used for feature analysis. 

Both scripts come with a wrapper to submit them on a sun grid engine to parallelize the computation. 


