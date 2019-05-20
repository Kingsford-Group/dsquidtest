import pandas as pd
from collections import defaultdict
import sys

def isdiscordant(n1, n2):
    if (n1[-1] == 'H' and n2[-1] == 'T') and (int(n1[:-1]) < int(n2[:-1])):
        return True
    elif n1[-1] == n2[-1]:
        return True
    else:
        return False

def node2Edge(ind1, ind2):
    n1 = ""
    n2 = ""
    if ind1<0:
        n1=str(abs(ind1))+"H"
    else:
        n1=str(abs(ind1))+"T"
    if ind1<0:
        n2=str(abs(ind2))+"T"
    else:
        n2=str(abs(ind2))+"H"
    return n1+","+n2

def read_conf_struct(fname):
    df = pd.read_csv(fname, sep="\t")
    all_struct = defaultdict(list)
    all_struct_dict = defaultdict(int)
    weights = defaultdict(int)
    for i in df.index:
        all_struct[df.loc[i,'struct_id']].append((df.loc[i,'n1'], df.loc[i,'n2']))
        all_struct_dict[df.loc[i, 'n1']+","+df.loc[i,'n2']] = df.loc[i,'struct_id']
        weights[df.loc[i, 'n1']+","+df.loc[i,'n2']] = df.loc[i,'weight']
    return all_struct, all_struct_dict, weights

def read_edge(fname):
    df = pd.read_csv(fname, sep="\t")
    all_edges = []
    for i in df.index:
        ind1 = df.loc[i, 'ind1']
        ind2 = df.loc[i, 'ind2']
        all_edges.append(node2Edge(ind1, ind2))
    return all_edges

if __name__ == "__main__":
    fname = sys.argv[1]
    all_struct, all_struct_dict, weights = read_conf_struct(fname)
    disc_set = set()
    portions = []
    # count all discordant edge
    for edge in all_struct_dict:
        n1 = edge.split(",")[0]
        n2 = edge.split(",")[1]
        if isdiscordant(n1,n2):
            disc_set.add(edge)
    #print(len(disc_set))

    # count portions of discordant edges
    avg_portion = 0
    for struct in all_struct:
        disc_set_2 = set()
        total_edge = 1
        for e in all_struct[struct]:
            if e[0][:-1] == e[1][:-1]:
                continue
            if isdiscordant(e[0], e[1]):
    #           print(e[0],e[1])
                disc_set_2.add(edge)
            total_edge +=1

        #print(len(disc_set_2)/len(all_struct[struct]))
    #    print(len(disc_set_2))
        avg_portion+= len(disc_set_2)/total_edge
        print(len(disc_set_2)/total_edge)
    #print(str(len(disc_set))+"\t"+str(avg_portion/len(all_struct)))
