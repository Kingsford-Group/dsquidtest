file="tcga_num_disc_conf_allpaths.txt"

echo -e "type\tsample\tnum_disc_conf\ttotal_disc\tfraction" > $file
for f in *TCGA*.txt
do
    type=`basename $f | cut -f 1 -d '_'`
    sample=`basename $f | cut -f 2 -d '_' `
    disc_conf=$(grep -v TOTAL $f | wc -l)
    total=$(head -1 $f | cut -f 2)
    echo -e $type'\t'$sample'\t'$disc_conf'\t'$total'\t'$(bc -l <<< "${disc_conf}/${total}")>> $file
done

