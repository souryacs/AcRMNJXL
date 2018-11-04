AcRNJ(XL) (species tree estimation using Accumulated coalescence Rank and optionally eXcess gene Lineage, with Neighbor Joining)
------------------

AcRNJ(XL) is a python based tool for computing species tree from a set of incongruent gene trees with Incomplete Lineage Sorting (ILS). AcRNJ proposes a novel couplet based distance measure, termed as the "accumulated coalescence rank" (AcR), 
to create a distance matrix, and construct a phylogenetic tree using NJ based species tree estimation method. The method 
is termed as AcRNJ.

A second feature, termed as the excess gene leaf count (XL), is also computed for individual couplets, and can 
be used together with the AcR measure. Using both AcR and XL measures lead to a method named AcRNJXL.

Description
-----------------------

Input
-----------

A collection of gene trees with overlapping set of species (from which the genes are sampled), having topological incongruence 
due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; both AcRNJ and AcRNJXL do not use the branch lengths of the gene trees for species tree estimation.

Input gene trees can be either in NEWICK format or in NEXUS format. However, all the gene trees should have identical input formats. They should be placed in a standard tree list file, according to the syntax of NEXUS or NEWICK formats. Such a tree list text file is to be provided as an input of this executable.

Output
--------

A species tree covering all the taxa of the gene trees. Output species tree is generated in the NEWICK format.

Installation
--------------

AcRNJ / AcRNJXL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems having linux OS (Fedora / Ubuntu). User can download the zipped archieve from GitHub, or clone the current repository. User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: Current version does not support python 3, and a feature release is planned to cover it.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: Future release is planned to cover the latest release of Dendropy, since Dendropy 4.0 (and onwards) made significant 
changes in their repositories.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest version.


Command line options
--------------

The python file "AcRNJXL.py" is the main executable. User should execute the file with the following command:

python AcRNJXL.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

	name of the input file containing gene trees

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

	name of the output file to contain target species tree

-m METHOD_type, --METHOD=METHOD_type

                1 - The method AcRNJ is to be executed. (default)
                2 - The method AcRNJXL is to be executed. (published in AlCOB 2016)

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

	1 - input file format is NEWICK (default)
	2 - input file format is NEXUS

-r TAXON_NAME, --ROOT=TAXON_NAME

	User can specify a taxon name to root the output species tree with that specified taxon. Useful to benchmark the 
	generated tree with respect to a given rooting configuration.

Example of a command (followed for the results published in the manuscript)

python AcRNJXL.py -I source_tree_input.txt -O 'output_file_name' -p1 -m1 -r "target_root_taxon_sample"

The method = 1 means AcRNJ will be used. Else the method AcRNJXL will be invoked.

If the output file name is not specified, a folder of name either "AcRNJ" (if the method = 1) or “AcRNJXL” (if method = 2) 
is created. A file named 'outtree_newick.tre' within the specified output directory contains the output species tree. 
A separate output log file 'Complete_Desription.txt' within the given output folder contains the detailed execution of the 
current method.


Citation
-----------

Upon using the package, user should cite the following papers:

1) Sourya Bhattacharyya, Jayanta Mukhopadhyay, Accumulated Coalescence Rank and Excess Gene Count for Species Tree Inference, proceedings of 3rd International Conference on Algorithms for Computational Biology (AlCoB), Trujillo, Spain, Springer LNBI 9702, pp. 93-105. 

2) Sourya Bhattacharyya, Jayanta Mukhopadhyay, Couplet Supertree based Species Tree Estimation  from Incongruent Gene Trees with Deep Coalescence, proceedings of 11th International Symposium on Bioinformatics Research and Applications (ISBRA), Virginia, USA, June 2015, Springer LNBI 9096, pages 48-59.


For any queries, please contact
---------------------------

1) Sourya Bhattacharyya 

Postdoctoral Researcher

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

<sourya.bhatta@gmail.com>


2) Jayanta Mukherjee

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

<jay@cse.iitkgp.ernet.in>





