*********************************
AcRNJXL
*********************************

******************************
Species tree estimation using mode based accumulated coalescence rank and extra gene leaves information
******************************

AcRNJXL is a python based tool for computing species tree from a set of incongruent gene trees 
with Incomplete Lineage Sorting (ILS). One of the following measures between individual couplets are compurted 
for species tree estimation.

A) Accumulated coalescence rank between individual couplets 

B) Accumulated extra gene leaves count between individual couplets.

These measures are averaged with respect to all input gene trees, and subsequently used to form the distance matrix 
D. It is then applied to the NJ method for species tree construction.


Description
-----------------------

Input
-----------

A collection of gene trees with overlapping taxa set (sampled genes), having topological incongruence 
due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; current species tree estimation 
does not consider the branch lengths of the tree.

Input gene trees can be either in NEWICK format or in NEXUS format. 
However, all the gene trees should have identical input formats. They should be placed in a 
standard tree list file, according to the syntax of NEXUS or NEWICK formats. Such a tree list 
text file is to be provided as an input of this executable.

Output
--------

A species tree covering all the taxa of the gene trees. Output species tree 
is generated in the NEWICK format.

*********************************
Dependencies
*********************************

AcRMNJXL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems 
having linux OS (Fedora / Ubuntu).

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to 
check the correct execution of our code, and optionally needs to upgrade it accordingly.
We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We 
did not upgrade the code for Dendropy 4.0 support, so any user having this new version of Dendropy 
might need to check the functionalities of AcRNJXL and possibly upgrade / replace / edit few 
dendropy related functions. So, we recommend users to use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest 
Numpy module in it. We found that Numpy module in the traditional Apt-Get repository is of lower version.

***************
Command line options
****************

After downloading the current archieve in a zipped format, extracting the archieve reveals a file AcRMNJXL.py which 
is the main executable file. Assuming current directory contains this python file, 
following command is to be executed:

./AcRNJXL [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

                name of the input file containing gene trees

-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

                name of the output file to contain target species tree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                1 - input file format is NEWICK (default)
                2 - input file format is NEXUS

-r TAXON_NAME, --ROOT=TAXON_NAME

		User can specify a taxon name to root the output species tree with that specified taxon.

Example of a command (followed for the results published in the manuscript)

./AcRMNJXL -I source_tree_input.txt -p1

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

In addition, the package contains another option: -O 'output_file_name'

Here, user can specify the output file name containing the derived species tree file.

the species tree can be rooted using a custom provided taxon, by using the option '-r'. 
Although AcRNJXL creates a rooted tree, custom tree rooting can be provided using this option.

If no such option is provided, a directory “AcRNJXL” is created within the same directory containing the input treelist file. 
Within this new created directory, one file 'outtree_newick.tre' is created, which contains the derived species tree. 
Another text file named 'Complete_Desription.txt' is created, which contains execution and timing information 
for the method. 

*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



