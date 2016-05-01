#!/usr/bin/env python

import Header
from Header import * 

##-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist

def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
	Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
		preserve_underscores=PRESERVE_UNDERSCORE, default_as_rooted=ROOTED_TREE)
	return Inp_TreeList

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

##--------------------------------------------------
#"""
#this function computes the entropy of a numpy array input
#"""
#def Compute_Entropy(labels, outfile):
	#if 1:	#(DEBUG_LEVEL >= 2):
		#fp = open(outfile, 'a')
	
	#""" 
	#Computes entropy of label distribution. 
	#"""
	#n_labels = len(labels)

	#if n_labels <= 1:
			#return 0

	#counts = numpy.bincount(labels)
	
	#if 1:	#(DEBUG_LEVEL >= 2):
		#fp.write('\n After bincount:   counts: ' + str(counts))
	
	#ent = 0
	
	#for i in range(len(counts)):
		#if (counts[i] > 0):
			#count_frac = ((counts[i] * 1.0) / n_labels)
			#log_count_frac = math.log(count_frac, 2)
			#ent = ent - (count_frac * log_count_frac)
			#if 1:	#(DEBUG_LEVEL >= 2):
				#fp.write('\n Nonzero count index: ' + str(i) + ' count_frac: ' + str(count_frac) + ' log_count_frac: ' + str(log_count_frac))
			
	#if 1:	#(DEBUG_LEVEL >= 2):
		#fp.close()

	#return ent

