dir=$1

for f in ${dir}/squid*graph.txt
do
    echo ${f}
    prefix=`basename $f | cut -f 2,3 -d '_'`
    echo ${prefix}
    mkdir -p ${dir}/conc_edges
    ~/SQUID/squid/test ${f} ${prefix} ${dir}/conc_edges
    ~/SQUID/diploid-squid/test ${f} ${prefix} ${dir}/conc_edges
done
