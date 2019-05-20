import networkx as nx
from collections import defaultdict
import pandas as pd
import os, sys
import glob
import time

def isdiscordant(n1, n2):
    if (n1[-1] == 'H' and n2[-1] == 'T') and (int(n1[:-1]) < int(n2[:-1])):
        return True
    elif n1[-1] == n2[-1]:
        return True
    else:
        return False

def equal(cycle1, cycle2):
    if len(cycle1) != len(cycle2):
        return False
    if set(cycle1) != set(cycle2):
        return False
    idx = cycle2.index(cycle1[0])
    cycle2_rot = cycle2[idx:]+cycle2[:idx]
    cycle2_rev = cycle2[::-1]
    idx_rev = cycle2_rev.index(cycle1[0])
    cycle2_rev_rot = cycle2_rev[idx_rev:]+cycle2_rev[:idx_rev]
    if cycle1 == cycle2_rev_rot or cycle2 == cycle2_rot:
        return True
    else:
        return False

def cycle_to_edge(cycle):
    return [(cycle[i], cycle[i+1]) for i in range(len(cycle)-1)]

# path is a list of edges
def test_conflict(path):
    dis = 0
    for e in path:
        if isdiscordant(e[0],e[1]):
            dis +=1
    if dis == 0:
        return False
    inter = []
    for e in path:
        if e[1][:-1] == e[0][:-1]:
            inter.append(0)
        else:
            inter.append(1)
    adj = 0
    for i in range(len(inter)-1):
        if inter[i] == 1 and inter[i+1] == 1:
            adj +=1
    if inter[-1] == 1 and inter[0] == 1:
        adj+=1
    if adj == 2:
        return False
    else:
        return True

def get_edges(graph_file):
    graph_df = pd.read_csv(graph_file,skiprows=2,sep="\t", header=None)
    edge_df = graph_df[graph_df[0]=='edge'].dropna(axis=1).reset_index()
    del edge_df['index']
    edge_df.columns=['edge','id','Ind1','Head1','Ind2','Head2','Weight']
    return edge_df

def build_graph(edge_df):
    G = nx.Graph()
    a = edge_df.apply(add_edges, G = G, axis=1)
    for node in set([n[:-1] for n in G.nodes()]):
        G.add_edge(node+"H", node+"T", weight=-1)
    return G
        
def add_edges(row, G):
    G.add_edge(str(row['Ind1'])+row['Head1'], str(row['Ind2'])+row['Head2'],weight=row['Weight'])
    
def get_conf_structures(G):
    # find all cycles in G
    G = G.to_directed()
    #all_cycles = [c for c in nx.simple_cycles(G) if len(c)>2]
    counter = 0
    unique_cycles = []
    conflict_structures = []
    added = 1
    for c in nx.simple_cycles(G):
        if len(c) > 2:
            added = 1
            for c2 in unique_cycles:
                if equal(c2, c):
                    added = 0
            if added == 1:
                unique_cycles.append(c)
                if test_conflict(cycle_to_edge(c+c[:1])):
                    conflict_structures.append(cycle_to_edge(c+c[:1]))
                    counter +=1
                    #if counter %100 == 0:
                    #    print (counter, len(unique_cycles))
        if counter > 5000:
            print("Upper limit of 5000 cycles reached.")
            break

    # Get all unique cycles
    #unique = [1 for i in all_cycles]
    #for i in range(len(all_cycles)):
    #    for j in range(i+1, len(all_cycles)):
    #        if unique[i] == 0:
    #            continue
    #        if equal(all_cycles[i], all_cycles[j]):
    #            unique[j] = 0
    #unique_cycles = [cycle_to_edge(all_cycles[i]+all_cycles[i][:1]) for i in range(len(all_cycles)) if unique[i] == 1]

#    conflict_structures = []
#    for c in unique_cycles:
#        if test_conflict(c):
#            conflict_structures.append(c)

    conflict_edges = set()
    for c in conflict_structures:
        for e in c:
            conflict_edges.add(e)
    return conflict_structures, conflict_edges

def get_weights(conflict_edges, edge_df):
    weights = []
    for e in conflict_edges:
        n1 = int(e[0][:-1])
        n2 = int(e[1][:-1])
        h1 = e[0][-1]
        h2 = e[1][-1]
        if n1 != n2:
            weight_df = edge_df[((edge_df['Ind1']==n1) & (edge_df['Ind2']==n2) & (edge_df['Head1']==h1) & (edge_df['Head2']==h2)) | ((edge_df['Ind1']==n2) & (edge_df['Ind2']==n1) & (edge_df['Head1']==h2) & (edge_df['Head2']==h1))]
            weights.append(weight_df['Weight'].tolist()[0])
    return weights

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("!!!! USAGE: python find_conflict_structures.py <input_graph> <outdir>")
		exit()
	fname = sys.argv[1]
	outdir= sys.argv[2]
	print(fname)
	edge_df = get_edges(fname)
	G = build_graph(edge_df)
	start = time.time()
	conflict_structures, conflict_edges = get_conf_structures(G)
	end = time.time()

	# write to file
	weights = nx.get_edge_attributes(G, 'weight')
	outfname = outdir+"/"+os.path.basename(fname).split(".")[0]+".confstruct.csv"
	print("Wrting output to "+ outfname)
	conf_file = open(outfname, 'w')
	conf_file.write("\t".join(["struct_id","n1","n2","weight"])+"\n")
	for i,struct in enumerate(conflict_structures):
		for e in struct:
			w = ""
			if e in weights:
				w = str(weights[e])
			elif (e[1],e[0]) in weights:
				w = str(weights[(e[1],e[0])])
			else:
				print("Error: edges (" + ",".join(e) + ") not found")
				exit() 
			conf_file.write("\t".join([str(i), e[0], e[1], w])+"\n")
	conf_file.close()
