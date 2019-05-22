#!/bin/python

import sys
import numpy as np
import copy
from tqdm import *
from GraphClass import *


def ShortestPath(newvNodes, vEdges, Labels, v1_id, head1, v2_id, head2):
	component_id = Labels[v1_id]
	# the newvNodes is modified for removing some of the connections in HeadEdges and TailEdges, but the edges are not removed in vEdges
	# for the shortest path, only access edges via HeadEdges and TailEdges, not from vEdges
	# vEdges is for collecting the other information, for example the other node the current node connects to, and the corresponding orientation
	dist = np.ones(2*len(newvNodes)) * np.inf
	prev = np.ones(2*len(newvNodes)) * np.nan
	# initialze the correponding head / tail of v1_id to have dist 0
	ind = v1_id * 2 + (1-head1)
	dist[ind] = 0
	Q = sum([[2*i, 2*i + 1] for i in range(len(newvNodes)) if Labels[i] == component_id], [])
	while len(Q) != 0:
		min_in_Q = np.argmin(dist[np.array(Q)])
		mindist_arg = Q[min_in_Q]
		v3_id = int(mindist_arg/2)
		head3 = (mindist_arg%2 == 0)
		# remove this node from Q
		Q = Q[:min_in_Q] + Q[(min_in_Q+1):]
		# update the distance and orientation of the other end of this node
		if dist[mindist_arg] + 1 < dist[v3_id * 2 + head3]:
			dist[v3_id * 2 + head3] = dist[mindist_arg] + 1
			prev[v3_id * 2 + head3] = v3_id * 2 + (1-head3)
		# update the distance and orientation of the nodes connecting to this end of v3_id
		edge_ids = newvNodes[v3_id].HeadEdges if head3 else newvNodes[v3_id].TailEdges
		for e_id in edge_ids:
			# collect the information of the other node and orientation
			v4_id = vEdges[e_id].Ind1 if vEdges[e_id].Ind1 != v3_id else vEdges[e_id].Ind2
			head4 = vEdges[e_id].Head1 if vEdges[e_id].Ind1 != v3_id else vEdges[e_id].Head2
			# compare with existing distance and update
			if dist[mindist_arg] + 1 < dist[v4_id * 2 + (1-head4)]:
				dist[v4_id * 2 + (1-head4)] = dist[mindist_arg] + 1
				prev[v4_id * 2 + (1-head4)] = v3_id * 2 + (1-head3)
	passed_nodes = []
	if not np.isinf(dist[v2_id * 2 + (1-head2)]):
		passed_nodes.append( v2_id * 2 + (1-head2))
		while passed_nodes[-1] != v1_id * 2 + 1-head1:
			x = passed_nodes[-1]
			passed_nodes.append( int(prev[x]) )
	# convert passed node to segment id and orientation
	passed_nodes = [(int(x/2), (x%2==0)) for x in passed_nodes]
	return dist[v2_id * 2 + (1-head2)], passed_nodes


def Find_Other_Path(vNodes, vEdges, Labels, e_id):
	v1_id = vEdges[e_id].Ind1
	v2_id = vEdges[e_id].Ind2
	component_id = Labels[v1_id]
	assert( component_id == Labels[v2_id] )
	assert( IsDiscordant(vNodes, vEdges, e_id) )
	# if the discordant edge appears in conflict, there must a loop passing over the same pair of nodes
	# picking e_id as the starting edge of the loop, the rest path must start / end with at least one different node compared to e_id
	is_in_conflict = False
	for head1 in [True, False]: # head1 represent the start orientation of the rest of loop
		for head2 in [True, False]: # head2 represent the end orientatino of the rest of loop
			# check whether at least one of them is different from the orientation of e_id
			if (head1, head2) != (vEdges[e_id].Head1, vEdges[e_id].Head2):
				# remove the edges connectng to not head1 in v1_id
				id_to_removed = [e_id]
				if head1:
					id_to_removed += vNodes[v1_id].TailEdges
				else:
					id_to_removed += vNodes[v1_id].HeadEdges
				if head2:
					id_to_removed += vNodes[v2_id].TailEdges
				else:
					id_to_removed += vNodes[v2_id].HeadEdges
			else:
				id_to_removed = [e_id]
			id_to_removed = set(id_to_removed)
			newvNodes = []
			for v_id in range(len(vNodes)):
				v = vNodes[v_id]
				if Labels[v_id] == component_id:
					newv = copy.deepcopy(v)
					newv.HeadEdges = [x for x in v.HeadEdges if not (x in id_to_removed)]
					newv.TailEdges = [x for x in v.TailEdges if not (x in id_to_removed)]
					newvNodes.append(newv)
				else:
					newvNodes.append(vNodes[v_id])
			dist, passed_nodes = ShortestPath(newvNodes, vEdges, Labels, v1_id, head1, v2_id, head2)
			# count number of duplicate
			if len(passed_nodes) != 0:
				map_count = {}
				for i in range(1, len(passed_nodes) - 1):
					if passed_nodes[i][0] in map_count:
						map_count[passed_nodes[i][0]] += 1
					else:
						map_count[passed_nodes[i][0]] = 1
				count = len([c for c in map_count.values() if c != 2])
				if passed_nodes[0][1] == vEdges[e_id].Head2:
					count += 1
				if passed_nodes[-1][1] == vEdges[e_id].Head1:
					count += 1
			if (not np.isinf(dist)) and count != 2:
				is_in_conflict = True
				break
	return is_in_conflict


def WriteDiscordantEdgesInConflict(outputfile, vNodes, vEdges, Labels):
	fp = open(outputfile, 'w')
	fp.write("# e_id\tInd1\tHead1\tInd2\tHead2\tWeight\n")
	with tqdm.tqdm(total=len(vEdges)) as pbar:
		for e_id in range(len(vEdges)):
			e = vEdges[e_id]
			if IsDiscordant(vNodes, vEdges, e_id):
				is_in_conflict = Find_Other_Path(vNodes, vEdges, Labels, e_id)
				if is_in_conflict:
					fp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(e_id, e.Ind1, e.Head1, e.Ind2, e.Head2, e.Weight))
			pbar.update(1)
	fp.close()


# vNodes, vEdges = ReadGraphFile("/pghbio/cure/congma/TCGAData/blcashort/TCGA-BL-A0C8/TCGA-BL-A0C8-01/StarAlign/sqnew_graph.txt")
vNodes, vEdges = ReadGraphFile(sys.argv[1])
Labels = ConnectedComponent(vNodes, vEdges)
WriteDiscordantEdgesInConflict(sys.argv[2], vNodes, vEdges, Labels)

# # test example
# v1 = Node_t(0, 0, 100)
# v2 = Node_t(0, 100, 200)
# v3 = Node_t(0, 200, 300)
# v4 = Node_t(0, 300, 400)
# vNodes = [v1, v2, v3, v4]
# e1 = Edge_t(0, 1, False, False, 1, 1)
# e2 = Edge_t(0, 2, False, False, 1, 1)
# e3 = Edge_t(1, 2, False, True, 1, 1)
# e4 = Edge_t(1, 3, True, True, 1, 1)
# vEdges = [e1, e2, e3, e4]

# vNodes, vEdges = UpdateNodeLink(vNodes, vEdges)
# Labels = ConnectedComponent(vNodes, vEdges)
# WriteDiscordantEdgesInConflict("test_cs_example1.txt", vNodes, vEdges, Labels)