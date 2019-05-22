import numpy as np
from collections import defaultdict
import pandas as pd
import math
import sys

def extract_bed(gtf_fname, out_name):
    f = open(gtf_fname)
    gene_region = defaultdict(list)
    for l in f:
        line = l.rstrip("\n").split("\t")
        transcript = line[-1].split(";")[0].split(" ")[-1][1:-1]
        gene = line[-1].split(";")[1].split(" ")[-1][1:-1]
        gene_region['chr'].append(line[0])
        gene_region['start'].append(line[3])
        gene_region['end'].append(line[4])
        gene_region['gene_name'].append(gene)
    gene_df = pd.DataFrame(gene_region)
#     gene_df.to_csv(out_name, index=None)
    return gene_df

def get_random_regions(gene_df, num_regions, out_name, dir_name, num):
    gene_df = gene_df.sort_values(by=['chr','start'])
    names = gene_df['gene_name'].unique()
    chosen = np.random.choice(names, num_regions)
    regions = defaultdict(list)
    for c in chosen:
        startrow = gene_df[gene_df['gene_name'] == c].head(1)
        endrow =gene_df[gene_df['gene_name'] == c].tail(1)
        regions['chr'].append(startrow['chr'].tolist()[0])
        regions['start'].append(int(startrow['start'].tolist()[0])-(num/100)*80000)
        regions['end'].append(int(endrow['end'].tolist()[0])+(num/100)*80000)
        regions['gene_name'].append(c)
    regions_df = pd.DataFrame(regions)
    regions_df.to_csv(dirname+"/"+out_name, index=None)
    return regions_df

def read_chrom_size(fname):
    chrom_dict = {}
    for l in open(fname):
        chrom = l.split("\t")[0]
        size = int(l.split("\t")[-1])
        chrom_dict[chrom] = size
    return chrom_dict

def SVconflict(start, size, chrom, chrom_size, all_loc):
    end = start + size
    if end >= chrom_size:
        return True
    for i in range(len(all_loc['chr'])):
        if all_loc['chr'][i] != chrom:
            continue
        if (start >= all_loc['start'][i] and start <= all_loc['end'][i]) or (end >= all_loc['start'][i] and end <= all_loc['end'][i]):
#             print(start, size, all_loc['start'][i], all_loc['end'][i])
            return True
    return False

# def add_dup(dup_st, chrom_dict, dup_size, all_loc, max_dup, dup_loc):
#     chrom = np.random.choice(list(chrom_dict.keys()),1)[0]
#     chrom_size = chrom_dict[chrom]

#     # initialize start position
#     dup_st = np.random.randint(1, chrom_size-1)
#     dup_size = dup_size * 
    
#     # check conflict
#     counter = 0
#     if 'chr' in all_loc:
#         while SVconflict(dup_st, dup_size, chrom, chrom_size, all_loc):
#             if counter > 10000000:
#                 print("Deletion simulation failed. Cannot find a position")
#                 print(chrom, dup_st, dup_size)
#             dup_st = np.random.randint(1, chrom_size)
#             counter += 1

#     dup_loc['chr'].append(chrom)
#     dup_loc['start'].append(dup_st)
#     dup_loc['end'].append(dup_st+dup_size)

#     all_loc['chr'].append(chrom)
#     all_loc['start'].append(dup_st)
#     all_loc['end'].append(dup_st+dup_size)
#     return dup_loc, all_loc

def add_tra(chrom_dict, all_loc, tra_loc):
# randomly choose 5' or 3' end, and randomly choose location based on the first/last SV
    chrom = np.random.choice(list(chrom_dict.keys()),1)[0]
    chrom_size = chrom_dict[chrom]
    chrom2 = np.random.choice(list(chrom_dict.keys()),1)[0]
    chrom_size2 = chrom_dict[chrom2]
    
    all_loc_df = pd.DataFrame(all_loc)
    first_1 = min(all_loc_df[all_loc_df['chr']==chrom]['start'].tolist())
    last_1 = max(all_loc_df[all_loc_df['chr']==chrom]['end'].tolist())
    first_2 = min(all_loc_df[all_loc_df['chr']==chrom2]['start'].tolist())
    last_2 = max(all_loc_df[all_loc_df['chr']==chrom2]['end'].tolist())
    while (first_1 == 1) or (first_2 == 1) or (last_1 >= chrom_size) or (last_2 >= chrom_size2):
        chrom = np.random.choice(list(chrom_dict.keys()),1)[0]
        chrom_size = chrom_dict[chrom]
        chrom2 = np.random.choice(list(chrom_dict.keys()),1)[0]
        chrom_size2 = chrom_dict[chrom2]
        first_1 = min(all_loc_df[all_loc_df['chr']==chrom]['start'].tolist())
        last_1 = max(all_loc_df[all_loc_df['chr']==chrom]['end'].tolist())
        first_2 = min(all_loc_df[all_loc_df['chr']==chrom2]['start'].tolist())
        last_2 = max(all_loc_df[all_loc_df['chr']==chrom2]['end'].tolist())
    
#     first = min(first_1, first_2)
#     last = max(last_1, last_2)
    
    start_1, end_1, start_2, end_2, size_1, size_2 = [0]*6
    
    head1 = np.random.randint(0,2)
    if head1 == 0:
        start_1 = 1
        end_1 = np.random.randint(2000, first_1)
#         print("End 1: " +str(end_1))
        size_1 = end_1 - start_1
    
    if head1 == 1:
        start_1 = np.random.randint(last_1, chrom_size)
        end_1 = chrom_size
        size = end_1 - start_1
    
    head2 = np.random.randint(0,2)
    if head2 == 0:
        start_2 = 1
        end_2 = np.random.randint(2000, first_2)
#         print("End 1: " +str(end_1))
        size_2 = end_2 - start_2
    if head2 == 1:
        start_2 = np.random.randint(last_2, chrom_size2)
        end_2 = chrom_size2
        size_2 = end_2 - start_2

    tra_loc['ChrA'].append(chrom)
    tra_loc['StartA'].append(start_1)
    tra_loc['EndA'].append(end_1)
    tra_loc['SizeA'].append(size_1)
    tra_loc['ChrB'].append(chrom2)
    tra_loc['StartB'].append(start_2)
    tra_loc['EndB'].append(end_2)
    tra_loc['SizeB'].append(size_2)
    
    all_loc['chr'].append(chrom)
    all_loc['start'].append(start_1)
    all_loc['end'].append(end_1)
    
    all_loc['chr'].append(chrom2)
    all_loc['start'].append(start_2)
    all_loc['end'].append(end_2)
#     print(chrom, start_1, end_1, chrom2, start_2, end_2, size_2)
    return tra_loc, all_loc

def add_del(chrom_dict, del_size, all_loc, del_loc, begin=0, end=0, chrom=0):
    if begin !=0 and end!=0 and chrom!=0:
        del_st = np.random.randint(begin-del_size, end)
        counter = 0
        while SVconflict(del_st, del_size, chrom, chrom_dict[chrom], all_loc):
            if counter > 100000:
                print("Deletion simulation failed in gene region: " + str(begin) +"," + str(end))
                print(del_size)
                return None
            del_st = np.random.randint(begin-del_size, end)
            counter +=1
    else:
        chrom = np.random.choice(list(chrom_dict.keys()),1)[0]
        chrom_size = chrom_dict[chrom]

        # initialize start position
        del_st = np.random.randint(5000, chrom_size-5000)

        while (del_st+del_size >= chrom_size):
            del_st = np.random.randint(5000, chrom_size-5000)

        # check conflict
        counter = 0
        if 'chr' in all_loc:
            while SVconflict(del_st, del_size, chrom, chrom_size, all_loc) or (del_st+del_size >= chrom_size):
                if counter > 100000:
                    print("Deletion simulation failed. Cannot find a position")
                    print(chrom, del_st, del_size)
                del_st = np.random.randint(5000, chrom_size-5000)
                counter += 1

        if del_st >= chrom_size or del_st+del_size >= chrom_size:
            print("?")

    del_loc['Chr'].append(chrom)
    del_loc['Start'].append(del_st)
    del_loc['End'].append(del_st+del_size)

    all_loc['chr'].append(chrom)
    all_loc['start'].append(del_st)
    all_loc['end'].append(del_st+del_size)
    return del_loc, all_loc

        
def add_ins(chrom_dict, ins_size, all_loc, ins_loc, begin=0, end=0, chrom=0):
    ins_st, ins_st2 = 0,0
    chrom2 = chrom
    if begin !=0 and end!=0 and chrom!=0:
        ins_st = np.random.randint(begin, end-del_size)
        ins_st2 = np.random.randint(begin, end-del_size)
        counter = 0
        while SVconflict(ins_st, ins_size, chrom, chrom_dict[chrom], all_loc):
            if counter > 100000:
                print("Insertion/Inversion simulation failed in gene region: " + str(begin) +"," + str(end))
                print(ins_size)
                return None
            
            ins_st = np.random.randint(begin-ins_size, end)
            counter +=1
            
        counter = 0
        while SVconflict(ins_st2, ins_size, chrom, chrom_dict[chrom], all_loc) or (ins_st2 > ins_st and ins_st2 <= ins_st+ins_size):
            if counter > 100000:
                print("Insertion/Inversion simulation failed in gene region: " + str(begin) +"," + str(end))
                print(chrom, ins_st2, ins_size)
                return None
            ins_st2 = np.random.randint(begin-ins_size, end)
            counter += 1
    
    
    else:
        chrom = np.random.choice(list(chrom_dict.keys()),1)[0]
        chrom_size = chrom_dict[chrom]
        chrom2 = np.random.choice(list(chrom_dict.keys()),1)[0]
        chrom_size2 = chrom_dict[chrom2]

        # initialize start position
        ins_st = np.random.randint(5000, chrom_size-5000)
        ins_st2 = np.random.randint(5000, chrom_size2-5000)

        while (ins_st+ins_size >= chrom_size) or (ins_st2+ins_size >= chrom_size):
            ins_st = np.random.randint(5000, chrom_size-5000)
            ins_st2 = np.random.randint(5000, chrom_size-5000)

        # check conflict
        counter = 0
        if 'chr' in all_loc:
            while SVconflict(ins_st, ins_size, chrom, chrom_size, all_loc) or (ins_st+ins_size >= chrom_size):
                if counter > 10000000:
                    print("Insertion/Inversion simulation failed. Cannot find a position")
                    print(chrom, ins_st, ins_size)
                ins_st = np.random.randint(5000, chrom_size-5000)
                counter += 1
        counter = 0
        if 'chr' in all_loc:
            while SVconflict(ins_st2, ins_size, chrom2, chrom_size2, all_loc) or (ins_st2 > ins_st and ins_st2 <= ins_st+ins_size) or (ins_st2+ins_size >= chrom_size2):
                if counter > 10000000:
                    print("Insertion/Inversion simulation failed. Cannot find a position")
                    print(chrom, ins_st2, ins_size)
                ins_st2 = np.random.randint(5000, chrom_size-5000)
                counter += 1

        if ins_st+ins_size >= chrom_size or ins_st2+ins_size >= chrom_size:
            print("?")

    ins_loc['ChrA'].append(chrom)
    ins_loc['StartA'].append(ins_st)
    ins_loc['EndA'].append(ins_st+ins_size)
    
    ins_loc['ChrB'].append(chrom2)
    ins_loc['StartB'].append(ins_st2)
    ins_loc['EndB'].append(ins_st2+ins_size)
    
    all_loc['chr'].append(chrom)
    all_loc['start'].append(ins_st)
    all_loc['end'].append(ins_st+ins_size)
    
    all_loc['chr'].append(chrom2)
    all_loc['start'].append(ins_st2)
    all_loc['end'].append(ins_st2+ins_size)
    return ins_loc, all_loc

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("!!! USAGE: python simulate_breakpoint.py <gtf> <chrom-sizes> <out-dir>")

    gtf_fname = sys.argv[1]
    genes_out_name = "genes_clean.bed"
    selected_out_name = "regions.bed"
    chrom_size = sys.argv[2]
    chrom_dict = read_chrom_size(chrom_size)
    num_regions = 6
    outdir=sys.argv[3]

    gene_df = extract_bed(gtf_fname, genes_out_name)

    for num in [100,250,400]:

        total_sv = num
        #region_df = get_random_regions(gene_df, num_regions, "SVRNA"+str(num)+"_rep1_"+selected_out_name, out_dir,num)

        for rep in range(1,7):
            if rep %2 ==1:
            	# ensure that every two replicates share the same random regions to induce heterogeneity
                region_df = get_random_regions(gene_df, num_regions, "SVRNA"+str(num)+"_rep"+str(rep)+"_"+selected_out_name, out_dir, num)
            print("Generating SVRNA"+str(total_sv)+"_"+str(rep))
            sv_size_df = pd.read_csv("SV_sizes/sizes_"+str(total_sv*2)+"_"+str(rep)+".csv", sep="\t")
            sv_size_df = sv_size_df.sample(frac=1).reset_index(drop=True)
            maxdup = 2

            all_loc = defaultdict(list)
            ins_loc = defaultdict(list)
            inv_loc = defaultdict(list)
            del_loc = defaultdict(list)
            dup_loc = defaultdict(list)
            tra_loc = defaultdict(list)

            # simulate background SVs
            for i in range(math.floor(total_sv * (1 - 0.05*num_regions))):
                ins_size = sv_size_df.loc[i, 'insSizes']
                inv_size = sv_size_df.loc[i, 'invSizes']
                del_size = sv_size_df.loc[i, 'delSizes']
                dup_size = sv_size_df.loc[i, 'dupSizes']

                # find a position on the chromosome such that the simulated SV does not overlap with any other SVs
            #     ins_st, ins_st2, inv_st, inv_st2, del_st, dup_st = np.random.randint(1, chrom_size-1,6)
                ins_loc, all_loc = add_ins(chrom_dict, ins_size, all_loc, ins_loc)
                inv_loc, all_loc = add_ins(chrom_dict, inv_size, all_loc, inv_loc)
                del_loc, all_loc = add_del(chrom_dict, del_size, all_loc, del_loc)

            #     dup_loc, all_loc = add_dup(dup_st, chrom, chrom_size, dup_size, all_loc, max_dup, dup_loc)
            for i in region_df.index:
                begin = region_df.loc[i,'start']
                end = region_df.loc[i,'end']
                chrom = region_df.loc[i,'chr']

                for j in range(math.ceil(total_sv * 0.05)):
                    index = math.floor(total_sv * (1 - 0.05*num_regions))
                    ins_size = sv_size_df.loc[index+j, 'insSizes']
                    inv_size = sv_size_df.loc[index+j, 'invSizes']
                    del_size = sv_size_df.loc[index+j, 'delSizes']
                    dup_size = sv_size_df.loc[index+j, 'dupSizes']

                    ins_loc, all_loc = add_ins(chrom_dict, ins_size, all_loc, ins_loc, begin, end, chrom)
                    inv_loc, all_loc = add_ins(chrom_dict, inv_size, all_loc, inv_loc, begin, end, chrom)
                    del_loc, all_loc = add_del(chrom_dict, del_size, all_loc, del_loc, begin, end, chrom)
            tra_loc, all_loc = add_tra(chrom_dict, all_loc, tra_loc)
            tra_loc, all_loc = add_tra(chrom_dict, all_loc, tra_loc)        
        
            ins_df = pd.DataFrame(ins_loc)
            ins_df['Name'] = ["insertion_"+str(i) for i in range(1, len(ins_df.index)+1)]
            cols = ins_df.columns.tolist()
            cols = [cols[-1]]+cols[:-1]
            ins_df.to_csv(out_dir+"/SVRNA"+str(total_sv)+"_"+str(rep)+"/insertions.csv", sep="\t", index=None, columns=cols)
            
            inv_df = pd.DataFrame(inv_loc)
            inv_df['Name'] = ["inversion"+str(i) for i in range(1, len(inv_df.index)+1)]
            cols = inv_df.columns.tolist()
            cols = [cols[-1]]+cols[:-1]
            inv_df.to_csv(out_dir+"/SVRNA"+str(total_sv)+"_"+str(rep)+"/inversions.csv", sep="\t", index=None, columns=cols)
            
            del_df = pd.DataFrame(del_loc)
            del_df['Name'] = ["deletion"+str(i) for i in range(1, len(del_df.index)+1)]
            cols = del_df.columns.tolist()
            cols = [cols[-1]]+cols[:-1]
            del_df.to_csv(out_dir+"/SVRNA"+str(total_sv)+"_"+str(rep)+"/deletions.csv", sep="\t", index=None, columns=cols)
            
            tra_df = pd.DataFrame(tra_loc)  
            tra_df['Name'] = ["translocation_"+str(i) for i in range(1, len(tra_df.index)+1)]
            cols = tra_df.columns.tolist()
            cols = [cols[-1]]+cols[:-1]
            tra_df.to_csv(out_dir+"/SVRNA"+str(total_sv)+"_"+str(rep)+"/translocations.csv", sep="\t", index=None, columns=cols)
