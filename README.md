The scripts in this git repository helps reproduce results in the following paper:

Detecting Transcriptomic Structural Variants in Heterogeneous Contexts via the Multiple Compatible Arrangements Problem 
[[WABI version]](http://drops.dagstuhl.de/opus/volltexte/2019/11048/)

You can install Diploid-SQUID from [here](https://github.com/Kingsford-Group/diploidsquid).

Please find the following scripts in their corresponding directories. Note that you might need to change the relative paths in the scripts.
### HCC dataset
- script:`HCC_pipeline/scripts/HCC_pipeline.sh`
- *_ref stores annotated fusion and non-fusion gene events
- To verify detected SVs, run `python3 VerifyFusionGene.py 1 <TrueSV.bedpe> <prefix_pred.txt>`

### Find conflict structures
- script: `conflict_structure_pipeline/Find_conflict.sh`

### Evaluate approximation algorithms
- script: `Approx_compare.sh`

