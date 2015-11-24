#!/usr/bin/env python

"""
program for species tree estimation
computes four criteria - 
1) Level / branch count information
2) Accumulated Rank statistics 
"""

import dendropy
from optparse import OptionParser
import math
import time
import os
import numpy
import sys

""" this dictionary defines the taxa pair relations
each entry is indexed by two nodes """
TaxaPair_Reln_Dict = dict()

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# variables depicting the method employed for species tree construction
# first two methods use average statistics of the coalescence rank or the branch count
# last two methods employ the mode statistics

# accumulated coalescence rank with average statistics
AcRNJ = 1
# accumulated coalescence rank with mode based statistics
AcRMNJ = 2
# mode based accumulated coalescence rank and extra lineage statistics
AcRMNJXL = 3

# variables used to denote whether we use traditional NJ method
# or use a variant of it, namely the agglomerative clustering
TRADITIONAL_NJ = 1
AGGLO_CLUST = 2

##-----------------------------------------------------
""" 
this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
"""
class Reln_TaxaPair(object):  
  def __init__(self):
    # this is the count of trees for which the couplet is supported
    self.tree_support_count = 0        
    # this list contains the accumulated sum of ranks of the internal nodes
    # present between individual couplets, for all the gene trees
    # corresponds to AcRNJ or AcRMNJ
    self.accumulated_rank_list = []
    # this is the extra lineage count list for this couplet
    self.XL_val_list = []            
            
  # this function adds the count of tree according to the support of 
  # corresponding couplet in the input tree
  def _IncrSupportTreeCount(self):
    self.tree_support_count = self.tree_support_count + 1
          
  # this function returns the number of trees supporting the couplet
  def _GetSupportTreeCount(self):
    return self.tree_support_count        
        
  def _AddXLVal(self, val):
    self.XL_val_list.append(val)

  def _GetAvgXLVal(self):
    return (sum(self.XL_val_list) * 1.0) / self.tree_support_count
          
  def _GetMultiModeXLVal(self):
    candidate_score_sum = 0
    candidate_freq_sum = 0
    curr_arr = numpy.array(self.XL_val_list)
    # returns the counts of individual elements
    # array size: max_elem + 1
    counts = numpy.bincount(curr_arr)
    # remove the zero values 
    values = numpy.nonzero(counts)[0]
    # mode value and corresponding frequency
    mode_val = numpy.argmax(counts)
    mode_count = numpy.max(counts)
    # check for the values having frequency at least half of the maximum frequency
    for v in values:
      if (counts[v] >= 0.5 * mode_count):
	candidate_score_sum = candidate_score_sum + (v * counts[v])
	candidate_freq_sum = candidate_freq_sum + counts[v]
    return (candidate_score_sum * 1.0) / candidate_freq_sum                
                            
  def _AddAccumulatedRank(self, val):
    self.accumulated_rank_list.append(val)
    
  def _GetAvgAccumulatedRank(self):
    return (sum(self.accumulated_rank_list) * 1.0) / self.tree_support_count
             
  def _GetMultiModeAccumulatedRank(self):
    candidate_score_sum = 0
    candidate_freq_sum = 0
    curr_arr = numpy.array(self.accumulated_rank_list)
    # returns the counts of individual elements
    # array size: max_elem + 1
    counts = numpy.bincount(curr_arr)
    # remove the zero values 
    values = numpy.nonzero(counts)[0]
    # mode value and corresponding frequency
    mode_val = numpy.argmax(counts)
    mode_count = numpy.max(counts)
    # check for the values having frequency at least half of the maximum frequency
    for v in values:
      if (counts[v] >= 0.5 * mode_count):
	  candidate_score_sum = candidate_score_sum + (v * counts[v])
	  candidate_freq_sum = candidate_freq_sum + counts[v]
    return (candidate_score_sum * 1.0) / candidate_freq_sum
                                             
  # this function prints information for the current couplet
  def _PrintTaxaPairRelnInfo(self, key, out_text_file, METHOD_USED):
    fp = open(out_text_file, 'a')    
    fp.write('\n taxa pair key: ' + str(key))
    fp.write('\n supporting number of trees: ' + str(self._GetSupportTreeCount()))
    if (METHOD_USED == AcRNJ):
      fp.write('\n *** Accumulated rank list: ' + str(sorted(self.accumulated_rank_list)))
      fp.write('\n *** average accumulated rank : ' + str(self._GetAvgAccumulatedRank()))   
    elif (METHOD_USED == AcRMNJ):
      fp.write('\n *** 50 percent mode accumulated rank : ' + str(self._GetMultiModeAccumulatedRank()))   
    elif (METHOD_USED == AcRMNJXL):
      fp.write('\n *** 50 percent mode accumulated rank : ' + str(self._GetMultiModeAccumulatedRank()))   
      fp.write('\n *** average XL val : ' + str(self._GetAvgXLVal()))   
    fp.close()
    
      
