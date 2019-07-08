
final_out_name=$1

for f in TCGA_graphs/*.txt
do
    prefix=`basename $f | cut -f 1,2 -d '_'`
    echo $prefix
    python3.6 find_conf_disc_all_paths.py $f TCGA_edges/${prefix}_allpaths_edge.txt
done

bash find_fraction.sh TCGA_graphs TCGA_edges $1
