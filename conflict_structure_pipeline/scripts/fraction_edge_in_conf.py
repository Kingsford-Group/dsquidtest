import pandas as pd
import os
import sys

def enumerate_edge_in_conf(cd_file, conf_file):
    # Concordant edges that were discordant in original genome
    edge_Cd_df = pd.read_csv(cd_file, sep="\t")
    # edges in conflict structures
    edge_conf_df = pd.read_csv(conf_file, sep="\t", skiprows=1, names=["Ind1","Head1","Ind2","Head2"], header=None)

    fraction = 0
    total = 0
    for i in edge_Cd_df.index:
        id1 = abs(int(edge_Cd_df.loc[i, 'ind1']))-1
        id2 = abs(int(edge_Cd_df.loc[i, 'ind2']))-1
        conf_edges = edge_conf_df[((edge_conf_df['Ind1']==id1) & (edge_conf_df['Ind2']==id2)) | ((edge_conf_df['Ind1']==id2) & (edge_conf_df['Ind2']==id1))]
        if len(conf_edges.index.tolist()) >0:
            fraction +=1
        total+=1
    fraction_in_Cd = fraction/total
    fraction_in_conf = fraction/len(edge_conf_df.index.tolist())
    return fraction_in_Cd, fraction_in_conf

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("!!! USAGE: python fraction_edge_in_conf.py <inputdir> <outname>")
        exit()

    directory = sys.argv[1]
    outname = sys.argv[2]
    edge_names_d = []
    edge_names_s = []
    conf_names = []
    samples = []
    for l in os.listdir(directory):
        if "diploid" in l and "C_Discordant" in l:
            edge_names_d.append(l)
            samples.append("_".join(l.split(".")[0].split("_")[3:]))
        elif "squid" in l and "C_Discordant" in l:
            edge_names_s.append(l)
        elif "edge.txt" in l:
            conf_names.append(l)
    edge_names_d = sorted(edge_names_d)
    edge_names_s = sorted(edge_names_s)
    conf_names = sorted(conf_names)

    fraction_dir = {'sample':samples,'diploid_frac_Cd':[], 'diploid_frac_conf':[], 'squid_frac_Cd':[], 'squid_frac_conf':[]}
    edge_names = [edge_names_d, edge_names_s]
    arr = ['diploid_frac_Cd','diploid_frac_conf','squid_frac_Cd', 'squid_frac_conf']
    for i in [0,1]:
        fractions_Cd = []
        fractions_conf = []
        for e, c in zip(edge_names[i], conf_names):
            fraction_Cd, fraction_conf = enumerate_edge_in_conf(directory+"/"+e,directory+"/"+c)
            fractions_Cd.append(fraction_Cd)
            fractions_conf.append(fraction_conf)
        fraction_dir[arr[i*2]] = fractions_Cd
        fraction_dir[arr[i*2+1]] = fractions_conf
    fraction_df = pd.DataFrame(fraction_dir)
    fraction_df.to_csv(outname,index=None, sep="\t")
