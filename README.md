# AcRMNJXL
#------------------
# Species tree estimation using mode based accumulated coalescence rank and extra lineage information
#-------------------

Implements species tree construction from incongruent gene trees with Incomplete Lineage Sorting (ILS). One of the following measures between individual couplets are compurted for species tree estimation. 

A) Accumulated coalescence rank between individual couplets

B) Accumulated extra lineage count between individual couplets.

Description
------------

Input
-------

A collection of gene trees with overlapping taxa set (sampled genes), having topological incongruence due to Incomplete Lineage Sorting (ILS). Gene trees may or may not be weighted; species tree estimation does not consider the branch lengths of the tree.

Output
------

A species tree covering all the taxa of the gene trees. The species tree is aimed to be topologically closer to the model species tree (if available), or to the input gene trees.


Methods implemented
---------------------

We have implemented following 3 kinds of methods for species tree estimation.

A) AcRNJ: Proposes accumulated coalescence rank between individual couplets, and computes them using individual gene trees. The average values of these accumulated rank information for individual couplets are used in NJ (Neighbor Joining) based implementation to produce the final species tree. The method is named as Accumulated coalescence Rank based NJ (AcRNJ).

B) AcRMNJ: Here we have used accumulated coalescence rank between individual couplets for species tree construction. However, instead of computing their average statistic over all the input gene trees, we consider only the accumulated rank values (for individual couplets) whose frequency of occurrence (in the gene trees) is at least 50% of the mode (coalescence rank) value. Corresponding subset of accumulated rank values are then averaged to use in NJ based species tree construction. The method is named as Accumulated coalescence Rank with Mode based NJ (AcRMNJ).

E) AcRMNJXL: Here also, we have used the mode statistic of accumulated coalescence rank between individual couplets. In addition, we have computed couplet based average extra lineage count, computed with respect to all the input gene trees. Product of these two measures is used for NJ based species tree construction. The method is named as Accumulated coalescence Rank with Mode based NJ and Extra Lineage information (AcRMNJXL).

Dependencies / Installation Requirements
------------------------------------------

This package is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. It is tested and meant for systems having linux OS (Fedora / Ubuntu).

Development was done using Python 2.7. Note: We plan to support Python 3 environment in some future release.

We have used the phylogenetic library Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ ) for implementation.

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We did not upgrade the code for Dendropy 4.0 support, and plan it as a future work.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

#********** 
User do not need to install anything. Please follow the instructions below to execute and use the implementation 
#****************

Execution
------------

The file AcRMNJXL.py is the main file in this package. It is to be executed with the following command line options, from a terminal. In terminal, go to the directory containing the source codes, and type the following commands:

chmod +x AcRMNJXL.py (To change its permission to make it an executable file)

./AcRMNJXL.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

                        name of the input file containing gene trees
  
-O OUT_FILENAME, --OUTFILE=OUT_FILENAME

                        name of the output file to contain target species tree

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                        1 - input file format is NEWICK (default)
                        2 - input file format is NEXUS

-m METHOD_TYPE, --method=METHOD_TYPE

                        1 - average accumulated coalescence rank of the
                        couplets (AcRNJ) (Default Method)
                        
                        2 - mode based accumulated coalescence rank of the
                        couplets (AcRMNJ)                    
                        
                        3 - product of accumulated coalescence rank between couplets, with 
                        couplet based average extra lineage information, for species tree inference (AcRMNJXL)

Example of a command (followed for the results published in the manuscript)

./AcRMNJXL.py -I source_tree_input.txt -p1 -m3

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, as specified by the option (-p1) (1 stands for newick)

4) -m option is used to specify the species tree construction method. Here 3 is used to denote that AcRMNJXL is employed. The value can vary from 1 to 3.

The output texts are printed at console. User can redirect the output results to any standard text file by using standard redirection operation (>). For example, in the above command, all the detailed results (textual descriptions) are redirected to file out.txt.

As mentioned in the command options, specification of output species tree containing file is not mandatory. In such a case, a folder named identical with the name of the method employed (such as AcRNJ / AcRMNJ / AcRMNJXL) will be created in the directory containing the input gene tree list file. Within that new created directory, a file 'outtree_Newick.tre' will contain the output species tree.

Utilities
---------

All three species tree construction methods associate O(MN^2 + N^3) time complexity and O(N^2) space complexity, for N input taxa and M input trees.

For any queries, please contact
---------------

Sourya Bhattacharyya

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

email: sourya.bhatta@gmail.com


