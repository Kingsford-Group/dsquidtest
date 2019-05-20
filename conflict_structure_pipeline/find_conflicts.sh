dir=$1
for f in ${dir}/squid*graph.txt
do
    python3.6 /home/yutongq/savanna/SQUID-simulation/true-SV-sim/scripts/find_conflict_structures.py ${f} ./conflict_structures
    python3.6 get_conf_discedges.py ${some_conf_structure}
done
