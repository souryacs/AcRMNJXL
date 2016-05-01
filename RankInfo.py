import Header
from Header import *  
import UtilFunc
from UtilFunc import *

#--------------------------------------------------------
# this function defines accumulated coalescence rank between two nodes
# via the path through their LCA node
def DefineAccCoalRankXL(xl_val, no_of_taxa, lca_node_rank, \
	node1_rank, node2_rank, node1, node2, METHOD_USED, FRACT_ACC_RANK):
  #--------------------------------------------------
	""" 
	derive the accumulated coalescence rank information
	it is the sum of coalescence rank measure between node1 to LCA node 
	and between node2 to the LCA node
	"""
	sum_acc_rank = 0
	if (node1_rank < lca_node_rank):
		sum_acc_rank = sum_acc_rank + (((lca_node_rank - node1_rank) * (lca_node_rank + node1_rank - 1)) / 2)
	if (node2_rank < lca_node_rank):
		sum_acc_rank = sum_acc_rank + (((lca_node_rank - node2_rank) * (lca_node_rank + node2_rank - 1)) / 2)
	# add the rank of the LCA node
	sum_acc_rank = sum_acc_rank + lca_node_rank
	
	##--------------------------------------------------
	## add - sourya
	#if (RANK_WRITE_DEBUG == 1):
		#AcR_Complete_List.append(sum_acc_rank)
		#Coal_Rank_Complete_List.append(lca_node_rank)
	## end add - sourya
	##--------------------------------------------------
	# add - sourya
	"""
	If fractional rank information is sought, 
	we divide the accumulated coalescence rank by square of the ntaxa
	and also normalize the LCA coalescence rank with ntaxa
	currently we use this fractional coalescence rank measure
	"""
	if (FRACT_ACC_RANK == 1):
		sum_acc_rank = sum_acc_rank * (1.0 / (no_of_taxa * no_of_taxa))
		lca_node_rank = (lca_node_rank * 1.0) / no_of_taxa	# add - sourya
	# end add - sourya
	
	key1 = (node1.taxon.label, node2.taxon.label)
	key2 = (node2.taxon.label, node1.taxon.label)
	
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		if (METHOD_USED == AcRNJ) or (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL):
			TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(sum_acc_rank)
		else:
			TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(lca_node_rank)
		if (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL) or \
			(METHOD_USED == ProdRNJXL) or (METHOD_USED == RankRNJXL):
			TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
		if (METHOD_USED == AcRNJ) or (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL):
			TaxaPair_Reln_Dict[key2]._AddAccumulatedRank(sum_acc_rank)
		else:
			TaxaPair_Reln_Dict[key2]._AddAccumulatedRank(lca_node_rank)
		if (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL) or \
			(METHOD_USED == ProdRNJXL) or (METHOD_USED == RankRNJXL):
			TaxaPair_Reln_Dict[key2]._AddXLVal(xl_val)
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		if (METHOD_USED == AcRNJ) or (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL):
			TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(sum_acc_rank)
		else:
			TaxaPair_Reln_Dict[key1]._AddAccumulatedRank(lca_node_rank)
		if (METHOD_USED == ProdAcRNJXL) or (METHOD_USED == RankAcRNJXL) or \
			(METHOD_USED == ProdRNJXL) or (METHOD_USED == RankRNJXL):
			TaxaPair_Reln_Dict[key1]._AddXLVal(xl_val)

	return

#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree, METHOD_USED, FRACT_ACC_RANK):
	
	# number of taxa in the current tree
	no_of_taxa = len(Curr_tree.infer_taxa().labels())
		
	# traverse the internal nodes of the tree in postorder fashion
	for curr_node in Curr_tree.postorder_internal_node_iter():
		# compute the rank associated with this node
		curr_node_rank = no_of_taxa - curr_node.level()
		curr_node_level = curr_node.level()
		# add - sourya
		"""
		for fractional excess gene leaf measure
		normalize it with the number of taxa of the input tree
		"""
		if (FRACT_ACC_RANK == 0):
			xl_val = len(curr_node.leaf_nodes()) - 2
		else:
			xl_val = (len(curr_node.leaf_nodes()) - 2) * (1.0 / no_of_taxa)
		# end add - sourya
		
		# list the leaf and internal children of the current node
		curr_node_child_leaf_nodes = []
		curr_node_child_internal_nodes = []
		for x in curr_node.child_nodes():
			if (x.is_leaf() == True):
				curr_node_child_leaf_nodes.append(x)
			else:
				curr_node_child_internal_nodes.append(x)
		
		# pair of leaf nodes will be related by sibling relations
		if (len(curr_node_child_leaf_nodes) > 1):
			for i in range(len(curr_node_child_leaf_nodes) - 1):
				for j in range(i+1, len(curr_node_child_leaf_nodes)):
					node1_rank = no_of_taxa - curr_node_child_leaf_nodes[i].parent_node.level()
					node2_rank = no_of_taxa - curr_node_child_leaf_nodes[j].parent_node.level()
					DefineAccCoalRankXL(xl_val, no_of_taxa, curr_node_rank, node1_rank, node2_rank, \
						curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], METHOD_USED, FRACT_ACC_RANK)
		
		# one leaf node (direct descendant) and another leaf node (under one internal node)
		# will be related by ancestor / descendant relations
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						node1_rank = no_of_taxa - p.parent_node.level()
						node2_rank = no_of_taxa - r.parent_node.level()
						DefineAccCoalRankXL(xl_val, no_of_taxa, curr_node_rank, \
							node1_rank, node2_rank, p, r, METHOD_USED, FRACT_ACC_RANK)      
				
		# finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							node1_rank = no_of_taxa - p.parent_node.level()
							node2_rank = no_of_taxa - q.parent_node.level()
							DefineAccCoalRankXL(xl_val, no_of_taxa, curr_node_rank, \
								node1_rank, node2_rank, p, q, METHOD_USED, FRACT_ACC_RANK)    
		

	return
