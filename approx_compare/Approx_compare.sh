#!/bin/bash 

#### CHANGE paths below ####
SQUID=~/SQUID/squid
DSQUID=~/SQUID/diploid-squid
####

indir=$1
outdir=$2
logdir=$3
outprefix=$4

time_file=${outprefix}"_runtime.txt"
rm ${time_file}
for graph in ${indir}/*graph.txt
do
    prefix=`basename ${graph} | cut -f 1,2 -d "_"`
    for op in D
    do
        start=`date +%s.%N`
        ${DSQUID}/test -${op} ${graph} ${prefix}_${op} ${outdir} > ${logdir}/${prefix}_${op}.log.txt | head -1 | cut -f 2
        end=`date +%s.%N`
        duration=`bc <<< ${end}-${start}`
        echo -e ${prefix}_${op}"\t"${duration} >> ${time_file}
    done
done

out_name=${outprefix}"_fraction.txt"
rm ${out_name}
#echo -e "name\tobj\tfrac"> $out_name
for f in ${logdir}/*log*txt
do
    obj=`grep OBJ $f |  awk -F "\t" '{SUM+=$2}END{print SUM}'`
    all=`tail -3 ${f} | head -1 | awk -F ": " '{print $2}'`
    res=`tail -3 ${f} | tail -1 | awk -F ": " '{print $2}'`
    frac=`bc -l <<< ${res}/${all}`
    prefix=`basename $f | cut -f 1 -d "."`
    echo -e "${prefix}\t${obj}\t${frac}" >> $out_name
done

sort -k 1 -V ${time_file} > ${outprefix}"_runtime.sorted.txt"
sort -k 1 -V ${out_name} > ${outprefix}"_fraction.sorted.txt"

join -j 1 -a 2 ${outprefix}"_runtime.sorted.txt" ${outprefix}"_fraction.sorted.txt"
