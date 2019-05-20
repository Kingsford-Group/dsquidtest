dir=$1
simdir=$2
script_dir=~/SQUID/dsquidtest/simulation_pipeline
for f in ${simdir}/SVRNA*
do
    prefix=`basename $f`
    echo ${prefix}
    echo $f
    echo "Getting Rearranged Transcripts from SV breakpoints..."
    python3.6 ${script_dir}/RearrangedTranscripts_2.py $f genome_original.fa genes_clean.gtf ${dir}/${prefix}
done

