#!/bin/bash

#### CHANGE paths below ####
diploid_path=~/SQUID/diploid-squid
squid_path=~/SQUID/squid
####

graph_dir=$1
edge_dir=$2
out_prefix=$3

if [ "$#" -ne 3 ]; then
    echo "!!!! USAGE: bash Find_conflict.sh <graphdir> <edgedir> <outprefix>"
    exit
fi

# For each discordant edge, determine if it is involved in a conflict structure
for f in ${in_graph_dir}/*.txt
do
    prefix=`basename $f | cut -f 1,2 -d '_'`
    echo $prefix
    python3.6 find_conf_disc_all_paths.py $f ${edge_dir}/${prefix}_allpaths_edge.txt
done

file=${out_prefix}"_num_disc_conf.txt"

echo -e "type\tsample\tnum_disc_conf\ttotal_disc\tfraction" > $file
for f in ${edge_dir}/*allpaths_edge.txt
do
    type=`basename $f | cut -f 1 -d '_'`
    sample=`basename $f | cut -f 2 -d '_' `
    disc_conf=$(grep -v TOTAL $f | wc -l)
    total=$(head -1 $f | cut -f 2)
    echo -e $type'\t'$sample'\t'$disc_conf'\t'$total'\t'$(bc -l <<< "${disc_conf}/${total}")>> $file
done

# Get fraction of edges
file=${out_prefix}"_frac_disc_conf.txt"
for graphfile in ${graphdir}/*txt
do
    prefix=`basename ${graphfile} | cut -f 1,2 -d '_'`
    echo ${prefix}
    # Useful output file: C_Discordant_squid/diploid_${prefix}.txt
    ${diploid_path}/test -D ${graphfile} ${prefix} ${outdir} > ${outdir}/diploid_${prefix}.log.txt
    ${squid_path}/test ${graphfile} ${prefix} ${outdir} > ${outdir}/squid_${prefix}.log.txt
done

python fraction_edge_in_conf.py ${outdir} ${file}
