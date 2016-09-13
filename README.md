*********************************
AcRNJ and AcRNJXL
*********************************

******************************
Species tree estimation using either one or both of the following parameters:

1) accumulated coalescence rank (primary feature)
2) excess gene leaf count
******************************

AcRNJ / AcRNJXL is a python based tool for computing species tree from a set of incongruent gene trees 
with Incomplete Lineage Sorting (ILS).

AcRNJ uses the accumulated coalescence rank for individual couplets, to construct a distance matrix, 
which is then used in the NJ based method to generate the species tree.

AcRNJXL uses both accumulated coalescence rank and the excess gene leaf count, for  individual 
couplets, to construct two different distance matrices. These matrices are then employed for NJ 
method based species tree generation.

Description
-----------------------

Input
-----------

A collection of gene trees with overlapping set of species (from which the genes are sampled), having topological incongruence 
due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; both AcRNJ and AcRNJXL
do not use the branch lengths of the gene trees for species tree estimation.

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

AcRNJ / AcRNJXL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems 
having linux OS (Fedora / Ubuntu).

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default) 

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to 
check the correct execution of our code, and optionally needs to upgrade it accordingly.
We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) 

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We 
did not upgrade the code for Dendropy 4.0 support, so any user having this new version of Dendropy 
might need to check the functionalities of AcRNJ / AcRNJXL and possibly upgrade / replace / edit few 
Dendropy related functions. So, we recommend users to use the earlier version 
of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest 
Numpy module in it. We found that Numpy module in the traditional apt-get repository is of lower version.

***************
Command line options
****************

After downloading the current archieve in a zipped format, extracting the archieve generates the 
folder "AcRNJXL_exec". Within that folder, there exists an executable file "AcRNJXL" (otherwise, user 
needs to apply the executable permission on that file).

Assuming that the current directory is within the folder AcRNJXL, following command is to be executed:

./AcRNJXL [options]

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

		User can specify a taxon name to root the output species tree with that specified taxon.

Example of a command (followed for the results published in the manuscript)

./AcRNJXL -I source_tree_input.txt -p1 -m1

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

4) The method AcRNJ is to be executed, since the -m option is set as 1.

In addition, the package contains another option: -O 'output_file_name'

Here, user can specify the output file name containing the derived species tree file.

the species tree can be rooted using a custom provided taxon, by using the option '-r'. 
Although AcRNJXL creates a rooted tree, custom tree rooting can be provided using this option.

If the method (-m option) is provided as 1, a folder "AcRNJ" is created within the same directory 
containing the input treelist file. If the method (-m option) is set as 2, a directory “AcRNJXL” is created 
instead. 

Within this new created directory, one file 'outtree_newick.tre' is created, which contains the derived species tree. 
Another text file named 'Complete_Desription.txt' is created, which contains execution and runtime information 
for the method. 

*********************************
For any queries, please contact
*********************************

Sourya Bhattacharyya 
Department of Computer Science and Engineering
Indian Institute of Technology Kharagpur
<sourya.bhatta@gmail.com>



