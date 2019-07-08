diploid_path=~/SQUID/diploid-squid
squid_path=~/SQUID/squid
graphdir=$1
outdir=$2
outname=$3

if [ "$#" -ne 3 ]; then
    echo "!!!! USAGE: bash get_fraction.sh <graphdir> <outdir> <outname>"
    exit
fi

for graphfile in ${graphdir}/*txt
do
    prefix=`basename ${graphfile} | cut -f 1,2 -d '_'`
    echo ${prefix}
    # Useful output file: C_Discordant_squid/diploid_${prefix}.txt
    ${diploid_path}/test -D ${graphfile} ${prefix} ${outdir} > ${outdir}/diploid_${prefix}.log.txt
    ${squid_path}/test ${graphfile} ${prefix} ${outdir} > ${outdir}/squid_${prefix}.log.txt
done

# ${outdir} should also contain files produced by find_conf_disc_all_paths.py
# ${outname} will be a tab-delimited file. Header is "${prefix} diploid_frac squid_frac"
python fraction_edge_in_conf.py ${outdir} ${outname}
