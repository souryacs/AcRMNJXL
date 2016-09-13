import Header
from Header import *
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
"""
this function is a shortcut to obtain the normalized expression 
used in the agglomerative clustering proposed in this code
as various methods are experimented, corresponding various forms of 
agglomerative clustering is tried
"""
#--------------------------------------------------------
def ObtainNormalizedVal(num, denom1, denom2):
	if ((denom1 + denom2) > 0):
		return (num * 1.0) / (denom1 + denom2)
	else:
		return 0

#---------------------------------------------
""" 
function to print the matrix content
N is the matrix dimension
"""
def PrintMatrixContent(N, TaxaList, inp_data, inp_str, textfile):
	fp = open(textfile, 'a')
	fp.write('\n\n ===>>> printing contents of ' + str(inp_str) + ' ---- ')
	for i in range(N):
		fp.write('\n ' + str(i) + '--' + str(TaxaList[i]) + '--->>')
		if (i > 0):
			for j in range(i):
				fp.write(' ' + str(inp_data[i][j]))
	fp.close()

#----------------------------------------------------
"""
this function fills the distance matrix using accumulated internode count
"""
def Fill_DistMat_CoalRankInfo(DistMat, type_distmat, ntaxa):
	for i in range(ntaxa - 1):
		for j in range(i+1, ntaxa):
			spec1 = COMPLETE_INPUT_TAXA_LIST[i]
			spec2 = COMPLETE_INPUT_TAXA_LIST[j]
			key1 = (spec1, spec2)
			key2 = (spec2, spec1)
			if key1 in TaxaPair_Reln_Dict:
				l = key1
			elif key2 in TaxaPair_Reln_Dict:
				l = key2
			else:
				l = None
			
			if l is not None:
				"""
				there exists a valid entry in the couplet dictionary
				find three different measures
				"""
				if (type_distmat == 1) or (type_distmat == 4) or (type_distmat == 5) or (type_distmat == 7) or (type_distmat == 8):
					avgR = TaxaPair_Reln_Dict[l]._GetAvgAccumulatedRank()
				if (type_distmat == 2) or (type_distmat == 4) or (type_distmat == 5) or (type_distmat == 6):
					medR = TaxaPair_Reln_Dict[l]._GetMedianAccumulatedRank()
				if (type_distmat == 3) or (type_distmat == 5) or (type_distmat == 6) or (type_distmat == 7) or (type_distmat == 8):
					modeR = TaxaPair_Reln_Dict[l]._GetMultiModeAccumulatedRank()

				# add - sourya
				if (type_distmat == 1):
					DistMat[i][j] = avgR
				elif (type_distmat == 2):
					DistMat[i][j] = medR
				elif (type_distmat == 3):
					DistMat[i][j] = modeR
				elif (type_distmat == 4):
					DistMat[i][j] = min(avgR, medR)
				elif (type_distmat == 5):
					DistMat[i][j] = min(avgR, medR, modeR)
				elif (type_distmat == 6):
					DistMat[i][j] = min(medR, modeR)
				elif (type_distmat == 7):
					DistMat[i][j] = (avgR + modeR) / 2.0
				elif (type_distmat == 8):
					DistMat[i][j] = min(avgR, modeR)
				# end add - sourya
				
				# symmetric property
				DistMat[j][i] = DistMat[i][j]
				
			else:
				"""
				there exists no valid couplet
				so, set the distance matrix values as (-1)
				"""
				DistMat[i][j] = -1
				DistMat[j][i] = DistMat[i][j]	# symmetric property
				
	return

#--------------------------------------------------------
"""
if excess gene count information is used, this 
function fills the distance matrix using average excess gene count
"""
def Fill_DistMat_ExcessGeneCount(DistMat, type_distmat, ntaxa):
	for i in range(ntaxa - 1):
		for j in range(i+1, ntaxa):
			spec1 = COMPLETE_INPUT_TAXA_LIST[i]
			spec2 = COMPLETE_INPUT_TAXA_LIST[j]
			key1 = (spec1, spec2)
			key2 = (spec2, spec1)
			if key1 in TaxaPair_Reln_Dict:
				l = key1
			elif key2 in TaxaPair_Reln_Dict:
				l = key2
			else:
				l = None
			
			if l is not None:
				"""
				there exists a valid entry in the couplet dictionary
				find three different measures
				"""
				if (type_distmat == 1) or (type_distmat == 4) or (type_distmat == 5) or (type_distmat == 7) or (type_distmat == 8):
					avgXL = TaxaPair_Reln_Dict[l]._GetAvgXLVal()
				if (type_distmat == 2) or (type_distmat == 4) or (type_distmat == 5) or (type_distmat == 6):
					medXL = TaxaPair_Reln_Dict[l]._MedianXLVal()
				if (type_distmat == 3) or (type_distmat == 5) or (type_distmat == 6) or (type_distmat == 7) or (type_distmat == 8):
					modeXL = TaxaPair_Reln_Dict[l]._GetMultiModeXLVal()
				
				# add - sourya
				if (type_distmat == 1):
					DistMat[i][j] = avgXL
				elif (type_distmat == 2):
					DistMat[i][j] = medXL
				elif (type_distmat == 3):
					DistMat[i][j] = modeXL
				elif (type_distmat == 4):
					DistMat[i][j] = min(avgXL, medXL)
				elif (type_distmat == 5):
					DistMat[i][j] = min(avgXL, medXL, modeXL)
				elif (type_distmat == 6):
					DistMat[i][j] = min(medXL, modeXL)
				elif (type_distmat == 7):
					DistMat[i][j] = (avgXL + modeXL) / 2.0
				elif (type_distmat == 8):
					DistMat[i][j] = min(avgXL, modeXL)
				# end add - sourya
	
				# symmetric property
				DistMat[j][i] = DistMat[i][j]
				
			else:
				"""
				there exists no valid couplet
				so, set the distance matrix values as (-1)
				"""
				DistMat[i][j] = -1
				DistMat[j][i] = DistMat[i][j]	# symmetric property
	
	return

#---------------------------------------------
"""
computing the row wise sum for individual taxca clusters
"""
def ComputeSumRowsDistMat(sum_list, nclust, DistMat, Output_Text_File, inpstr):
	for i in range(nclust):
		t1 = 0
		for j in range(nclust):
			if (DistMat[i][j] >= 0):	# add the condition - sourya
				t1 = t1 + DistMat[i][j]
		sum_list.append(t1)
		
	if (DEBUG_LEVEL > 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n content of ' + str(inpstr) + ' : ' + str(sum_list))
		fp.close()
	
	return

#----------------------------------------------------------
"""
fill the normalized matrix entries for agglomerative clustering based method
"""
def FillAggloClustNormalizeMatrix(NormMat, DistMat, sum_list, nclust):
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			# comment - sourya
			#NormMat[i][j] = ObtainNormalizedVal(DistMat[i][j], sum_list[i], sum_list[j])
			# add - sourya
			NormMat[i][j] = ObtainNormalizedVal(DistMat[i][j], \
				(sum_list[i] - DistMat[i][j]), (sum_list[j] - DistMat[i][j]))
			# end add - sourya
			NormMat[j][i] = NormMat[i][j]

	return

#----------------------------------------------------------
"""
fill the normalized matrix entries for NJ based method
"""
def FillNJNormalizeMatrix(NormMat, DistMat, sum_list, nclust):
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			# comment - sourya
			ri = sum_list[i] / (nclust - 2)
			rj = sum_list[j] / (nclust - 2)
			## add - sourya
			#ri = (sum_list[i] - DistMat[i][j]) / (nclust - 2)
			#rj = (sum_list[j] - DistMat[i][j]) / (nclust - 2)
			## end add - sourya
			NormMat[i][j] = (DistMat[i][j] - ri - rj)
			NormMat[j][i] = NormMat[i][j]

	return

#---------------------------------------------
"""
finds the minimum of the distance matrix
when the product of both accumulated coalescence rank and excess gene count based 
measures are used
"""
def Find_Unique_Min_XL(DistMat_CoalRank, Norm_DistMat_CoalRank, DistMat_XL, \
	Norm_DistMat_XL, nclust, clust_species_list, NJ_RULE_USED):
	
	target_val = (Norm_DistMat_CoalRank[0][1] * Norm_DistMat_XL[0][1])
	min_idx_i = 0
	min_idx_j = 1
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (NJ_RULE_USED == AGGLO_CLUST):
				if ((Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j]) < target_val):
					target_val = (Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j])
					min_idx_i = i
					min_idx_j = j
			else:
				if ((Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j]) > target_val):
					target_val = (Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j])
					min_idx_i = i
					min_idx_j = j
				elif (FlEq((Norm_DistMat_CoalRank[i][j] * Norm_DistMat_XL[i][j]), target_val) == True):
					"""
					# condition add - sourya - 30.08.2016
					equal value with respect to the earlier minimum
					here we agglomerate the clusters according to the following condition:
					1) if DistMat_CoalRank[i][j] is strictly lower then use the new cluster pair
					2) if DistMat_CoalRank[i][j] = DistMat_CoalRank[min_idx_i][min_idx_j], 
					and DistMat_XL[i][j] is strictly lower then use the new cluster pair
					3) if both Acc Coal Rank count and XL measures are identical for the above mentioned clusters, 
					agglomerate the new cluster pair, if they have higher sum of cardinality
					"""
					if (DistMat_CoalRank[i][j] < DistMat_CoalRank[min_idx_i][min_idx_j]):
						min_idx_i = i
						min_idx_j = j
					else:
						if (FlEq(DistMat_CoalRank[i][j], DistMat_CoalRank[min_idx_i][min_idx_j]) == True):
							if (DistMat_XL[i][j] < DistMat_XL[min_idx_i][min_idx_j]):
								min_idx_i = i
								min_idx_j = j
							else:
								if (FlEq(DistMat_XL[i][j], DistMat_XL[min_idx_i][min_idx_j]) == True):
									"""
									higher sum of cardinality of the new cluster pair
									"""
									if (len(clust_species_list[i]) + len(clust_species_list[j])) > (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
										min_idx_i = i
										min_idx_j = j
				
	return min_idx_i, min_idx_j
	
#---------------------------------------------
"""
finds the minimum of the distance matrix
when only branch count based measure is used
"""
def Find_Unique_Min(Norm_DistMat_CoalRank, DistMat_CoalRank, nclust):
	target_val = Norm_DistMat_CoalRank[0][1]
	min_idx_i = 0
	min_idx_j = 1
	for i in range(nclust - 1):
		for j in range(i+1, nclust):
			if (i == j):
				continue
			if (DistMat_CoalRank[i][j] >= 0):	# condition add - sourya
				if (Norm_DistMat_CoalRank[i][j] < target_val):
					target_val = Norm_DistMat_CoalRank[i][j]
					min_idx_i = i
					min_idx_j = j
	
	return min_idx_i, min_idx_j

#-------------------------------------------
"""
checks whether a taxa cluster specified by the input index is a leaf
"""
def IsLeafCluster(clust_species_list, idx):
	if (len(clust_species_list[idx]) == 1):
		return True
	return False

#-------------------------------------------
"""
this function has following parameters:
1) first_cluster_mrca_node: root of 1st subtree 
2) second_cluster_mrca_node: root of 2nd subtree 
3) all_taxa_mrca_node: root of all these trees
4) Curr_tree: Tree containing all these subtrees

It creates one new internal node as a child of all_taxa_mrca_node
and places above mentioned subtrees as its children
"""
def MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
	all_taxa_mrca_node, taxa_list, Output_Text_File):
	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
		fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
		fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
		fp.close()

	# create new internal node 
	newnode = dendropy.Node()
	
	# its parent node will be the previous MRCA node of all the taxa in two clusters
	all_taxa_mrca_node.add_child(newnode)
	newnode.parent_node = all_taxa_mrca_node
	all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = None
	all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = None
	
	# add these individual clusters' MRCA node as its children
	newnode.add_child(first_cluster_mrca_node)
	first_cluster_mrca_node.parent_node = newnode
	newnode.add_child(second_cluster_mrca_node)
	second_cluster_mrca_node.parent_node = newnode
	
	# update splits of the resulting tree
	Curr_tree.update_splits(delete_outdegree_one=False)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
		fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(Curr_tree.mrca(taxon_labels=taxa_list))))  
		fp.close()

	return Curr_tree

#-------------------------------------------
"""
this function merges a pair of clusters whose indices are pointed by the min_idx_i and min_idx_j entries
this is part of the proposed agglomerative clustering
taxa_list is the union of these two clusters (species contents)
"""
def Merge_Cluster_Pair(Curr_tree, clust_species_list, min_idx_i, min_idx_j, taxa_list, Output_Text_File):
	isleaf_clust1 = IsLeafCluster(clust_species_list, min_idx_i)
	isleaf_clust2 = IsLeafCluster(clust_species_list, min_idx_j)
	
	if (isleaf_clust1):
		first_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
	else:
		first_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_i])
	
	if (isleaf_clust2):
		second_cluster_mrca_node = Curr_tree.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
	else:
		second_cluster_mrca_node = Curr_tree.mrca(taxon_labels=clust_species_list[min_idx_j])
	
	all_taxa_mrca_node = Curr_tree.mrca(taxon_labels=taxa_list)
	
	Curr_tree = MergeSubtrees(Curr_tree, first_cluster_mrca_node, second_cluster_mrca_node, \
		all_taxa_mrca_node, taxa_list, Output_Text_File)
	
	return Curr_tree

#---------------------------------------------
"""
this function refines the initial species tree (in terms of a star network) to 
find the true species tree
it does using agglomerative clustering (NJ principle)
the distance metric employed for NJ algorithm can vary depending on experimentation 
"""
def Form_Species_Tree_NJ_Cluster(Curr_tree, COMPLETE_INPUT_TAXA_LIST, METHOD_USED, \
	NJ_RULE_USED, Output_Text_File, XL_DIST_MAT_TYPE, ACC_RANK_DIST_MAT_TYPE):

	"""
	initially we have N of clusters for N taxa, where individual clusters are isolated
	agglomerating technique introduces a bipartition (speciation) which contains two taxa as its children
	"""
	no_of_taxa_clust = len(COMPLETE_INPUT_TAXA_LIST)

	# initialize the taxa clusters
	# copying the taxa list is done since initial clusters contain single species  
	# comment - sourya - we do not just copy ordinarily
	#clust_species_list = COMPLETE_INPUT_TAXA_LIST[:]
	# add - sourya - we enclose individual elements within a list and then copy these single element lists
	clust_species_list = []
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		subl = []
		subl.append(COMPLETE_INPUT_TAXA_LIST[i])
		clust_species_list.append(subl)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n\n\n COMPLETE_INPUT_TAXA_LIST ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp.write('\n\n Initial formed clust_species_list ' + str(clust_species_list))
		fp.close()        

	"""
	allocate a 2D square matrix of no_of_taxa_clust dimension
	for a pair of taxa clusters Cx and Cy, it contains the employed AcR measure for the cluster pairs
	"""
	Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	allocate one new square matrix which will contain the NJ based relative distance values between the 
	(used for minimum finding routine)
	"""
	Norm_Mean_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)
	
	if (METHOD_USED == ProdAcRNJXL):
		"""
		we allocate another matrix which will contain the XL measure among individual couplets
		"""
		XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)
		"""
		similarly, allocate one new square matrix to contain the relative XL based distance measure 
		for individual couplets
		"""
		Norm_XLVal_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	fill input distance matrices using accumulated coalescence rank information
	"""
	Fill_DistMat_CoalRankInfo(Mean_DistMat_ClustPair_NJ, ACC_RANK_DIST_MAT_TYPE, no_of_taxa_clust)

	"""
	for the method ProdAcRNJXL, fill the XL based distance matrix as well
	"""
	if (METHOD_USED == ProdAcRNJXL): 
		Fill_DistMat_ExcessGeneCount(XLVal_DistMat_ClustPair_NJ, XL_DIST_MAT_TYPE, no_of_taxa_clust)

	#--------------------------------------------------------
	# loop to execute the agglomerative clustering
	#--------------------------------------------------------
	while(no_of_taxa_clust > 2): 

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n\n\n iteration start --- number of clusters: ' + str(no_of_taxa_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.close()
			"""
			printing the AcR matrix 
			"""
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Mean_DistMat_ClustPair_NJ, \
				'Mean_DistMat_ClustPair_NJ', Output_Text_File)
			"""
			for the method ProdAcRNJXL, printing the XL matrix 
			"""
			if (METHOD_USED == ProdAcRNJXL):
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, XLVal_DistMat_ClustPair_NJ, \
					'XLVal_DistMat_ClustPair_NJ', Output_Text_File)
					
		"""
		for individual cluster Cx, it contains AcR(Cx, :) - sum of AcR measures for all 
		the cluster pairs (Cx, Cy) for all other clusters Cy
		"""
		sum_DistMat_Clust = []
		ComputeSumRowsDistMat(sum_DistMat_Clust, no_of_taxa_clust, Mean_DistMat_ClustPair_NJ, \
			Output_Text_File, 'sum_DistMat_Clust')
		
		"""
		for the method ProdAcRNJXL, compute for individual cluster Cx, the measure XL(Cx, :) - 
		sum of XL measures for all 
		the cluster pairs (Cx, Cy) for all other clusters Cy
		"""
		if (METHOD_USED == ProdAcRNJXL):
			sum_XLVal_Clust = []
			ComputeSumRowsDistMat(sum_XLVal_Clust, no_of_taxa_clust, XLVal_DistMat_ClustPair_NJ, \
				Output_Text_File, 'sum_XLVal_Clust')
		
		"""
		fill the normalized (relative) distance matrix entries for NJ based computation
		we have used two different mechanisms for computing such relative distance matrix
		1) Traditional NJ based relative distance matrix
		2) Agglomerative clustering based distance matrix
		"""
		if (NJ_RULE_USED == AGGLO_CLUST):
			FillAggloClustNormalizeMatrix(Norm_Mean_DistMat_ClustPair_NJ, \
				Mean_DistMat_ClustPair_NJ, sum_DistMat_Clust, no_of_taxa_clust)
			if (METHOD_USED == ProdAcRNJXL):
				FillAggloClustNormalizeMatrix(Norm_XLVal_DistMat_ClustPair_NJ, \
					XLVal_DistMat_ClustPair_NJ, sum_XLVal_Clust, no_of_taxa_clust)
		else:
			FillNJNormalizeMatrix(Norm_Mean_DistMat_ClustPair_NJ, Mean_DistMat_ClustPair_NJ, \
				sum_DistMat_Clust, no_of_taxa_clust)
			if (METHOD_USED == ProdAcRNJXL):
				FillNJNormalizeMatrix(Norm_XLVal_DistMat_ClustPair_NJ, XLVal_DistMat_ClustPair_NJ, \
					sum_XLVal_Clust, no_of_taxa_clust)

		"""
		printing the NJ based relative distance measures
		"""
		if (DEBUG_LEVEL >= 2):
			PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_Mean_DistMat_ClustPair_NJ, \
				'Norm_Mean_DistMat_ClustPair_NJ', Output_Text_File)
			if (METHOD_USED == ProdAcRNJXL): 
				PrintMatrixContent(no_of_taxa_clust, clust_species_list, Norm_XLVal_DistMat_ClustPair_NJ, \
					'Norm_XLVal_DistMat_ClustPair_NJ', Output_Text_File)
			
		#--------------------------------------------------
		"""
		find the cluster pairs having minimum distance values
		"""
		if (METHOD_USED == ProdAcRNJXL):
			"""
			this is earlier product of XL and Rank and minimum finding based method
			"""
			min_idx_i, min_idx_j = Find_Unique_Min_XL(Mean_DistMat_ClustPair_NJ, Norm_Mean_DistMat_ClustPair_NJ, \
				XLVal_DistMat_ClustPair_NJ, Norm_XLVal_DistMat_ClustPair_NJ, \
					no_of_taxa_clust, clust_species_list, NJ_RULE_USED)
		else:
			min_idx_i, min_idx_j = Find_Unique_Min(Norm_Mean_DistMat_ClustPair_NJ, Mean_DistMat_ClustPair_NJ, no_of_taxa_clust)
		
		#---------------------------------------------------------------  
		"""
		note down the taxa list in these two indices (min_idx_i and min_idx_j) 
		of the clust_species_list
		"""
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete taxa list (union) ' + str(taxa_list))
			fp.close()
		#----------------------------------------------------
		"""
		now we merge the pair of clusters pointed by these indices
		"""
		Curr_tree = Merge_Cluster_Pair(Curr_tree, clust_species_list, \
			min_idx_i, min_idx_j, taxa_list, Output_Text_File)
		#---------------------------------------------------------------------
		"""
		adjust the Mean_DistMat_ClustPair_NJ by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		Mean_DistMat_ClustPair_NJ = numpy.vstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		# then append one column
		Mean_DistMat_ClustPair_NJ = numpy.hstack((Mean_DistMat_ClustPair_NJ, \
			numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		# now reshape the distance matrix
		Mean_DistMat_ClustPair_NJ = numpy.reshape(Mean_DistMat_ClustPair_NJ, \
			((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		"""
		if the method used is ProdAcRNJXL
		apply these operations on the XL based distance matrix as well
		"""
		if (METHOD_USED == ProdAcRNJXL):
			XLVal_DistMat_ClustPair_NJ = numpy.vstack((XLVal_DistMat_ClustPair_NJ, \
				numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
			XLVal_DistMat_ClustPair_NJ = numpy.hstack((XLVal_DistMat_ClustPair_NJ, \
				numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
			XLVal_DistMat_ClustPair_NJ = numpy.reshape(XLVal_DistMat_ClustPair_NJ, \
				((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		"""
		add taxa_list as a new element of clust_species_list
		"""
		clust_species_list.append(taxa_list)          
		
		"""
		update the distance values with respect to this merged new cluster
		to all other clusters
		according to the NJ based method
		"""
		for m in range(no_of_taxa_clust):
			if (m == min_idx_i) or (m == min_idx_j):
				continue
			#-------------------------------
			"""
			we have used simple averaging to approximate the accumulated coalescence rank measure
			this averaging produces the best performance
			"""
			# modified - sourya
			if (Mean_DistMat_ClustPair_NJ[min_idx_i][m] >= 0) and (Mean_DistMat_ClustPair_NJ[min_idx_j][m] < 0):
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = Mean_DistMat_ClustPair_NJ[min_idx_i][m]
			elif (Mean_DistMat_ClustPair_NJ[min_idx_i][m] < 0) and (Mean_DistMat_ClustPair_NJ[min_idx_j][m] >= 0):
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = Mean_DistMat_ClustPair_NJ[min_idx_j][m]
			else:
				# this was used at the latest code
				Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
					(Mean_DistMat_ClustPair_NJ[min_idx_i][m] + Mean_DistMat_ClustPair_NJ[min_idx_j][m]) / 2.0
			# end modification - sourya
			
			# symmetric property
			Mean_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = Mean_DistMat_ClustPair_NJ[no_of_taxa_clust][m]
			#-------------------------------
			if (METHOD_USED == ProdAcRNJXL): 
				#----------------------------------------------
				# update the matrix containing excess gene count
				#----------------------------------------------
				# modified - sourya
				if (XLVal_DistMat_ClustPair_NJ[min_idx_i][m] >= 0) and (XLVal_DistMat_ClustPair_NJ[min_idx_j][m] < 0):
					XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = XLVal_DistMat_ClustPair_NJ[min_idx_i][m]
				elif (XLVal_DistMat_ClustPair_NJ[min_idx_i][m] < 0) and (XLVal_DistMat_ClustPair_NJ[min_idx_j][m] >= 0):
					XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = XLVal_DistMat_ClustPair_NJ[min_idx_j][m]
				else:
					XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m] = \
						(XLVal_DistMat_ClustPair_NJ[min_idx_i][m] + XLVal_DistMat_ClustPair_NJ[min_idx_j][m]) / 2.0
				# end modification - sourya
				
				# symmetric property
				XLVal_DistMat_ClustPair_NJ[m][no_of_taxa_clust] = XLVal_DistMat_ClustPair_NJ[no_of_taxa_clust][m]
				#-------------------------------
				
		"""
		now remove the rows and columns corresponding to min_idx_i and min_idx_j
		"""
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		Mean_DistMat_ClustPair_NJ = numpy.delete(Mean_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		if (METHOD_USED == ProdAcRNJXL):
			"""
			remove the entries in the XL based matrix as well
			"""
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=0)	# delete the row
			XLVal_DistMat_ClustPair_NJ = numpy.delete(XLVal_DistMat_ClustPair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		"""
		clear the entries in the NJ based relative distance matrices as well
		"""
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_Mean_DistMat_ClustPair_NJ = numpy.delete(Norm_Mean_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_Mean_DistMat_ClustPair_NJ.fill(0)
		
		if (METHOD_USED == ProdAcRNJXL):
			"""
			clear the entries in the NJ based relative distance matrices as well
			"""
			Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
			Norm_XLVal_DistMat_ClustPair_NJ = numpy.delete(Norm_XLVal_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
			Norm_XLVal_DistMat_ClustPair_NJ.fill(0)    
		
		"""
		remove individual clusters' taxa information from the clust_species_list
		"""
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		"""
		decrement the number of clusters considered
		"""
		no_of_taxa_clust = no_of_taxa_clust - 1

	return
