SQUID=~/SQUID/squid
DSQUID=~/SQUID/diploid-squid

indir=$1
outdir=$2
logdir=$3

for graph in ${indir}/*graph.txt
do
    prefix=`basename ${graph} | cut -f 1,2 -d "_"`
    for op in D
    do
        start=`date +%s.%N`
        ${DSQUID}/test -${op} ${graph} ${prefix}_${op} ${outdir} > ${logdir}/${prefix}_${op}.log.txt | head -1 | cut -f 2
        end=`date +%s.%N`
        duration=`bc <<< ${end}-${start}`
        echo -e ${prefix}_${op}"\t"${duration}
    done
done
