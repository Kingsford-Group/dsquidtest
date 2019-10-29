The scripts in this git repository helps reproduce results in the following paper:

Detecting Transcriptomic Structural Variants in Heterogeneous Contexts via the Multiple Compatible Arrangements Problem 

### HCC dataset
- script:`HCC_pipeline/scripts/HCC_pipeline.sh`
- Scripts for downloading, preprocessing and running SQUID and D-SQUID on HCC are in HCC_pipeline/scripts
- *_ref stores annotated fusion and non-fusion gene events
- To verify detected SVs, run `python3 VerifyFusionGene.py 1 <TrueSV.bedpe> <prefix_pred.txt>`

### Find conflict structures
- script: `conflict_structure_pipeline/Find_conflict.sh`

### Evaluate approximation algorithms
- script: `Approx_compare.sh`

