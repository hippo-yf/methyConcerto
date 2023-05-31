# methyConcerto

**Refined data analyses for single-cell bisulfite-based sequencing**


## Empirical Bayesian Methylation Caller (EBMC)

**Determinate discrete methyation status for each cytosine from C/T read counts**


## BS-SNV-Caller

**SNV calling with bisuifite-converted sequencing data**

### usage

`zcat ./data/sample.ATCG.gz | awk '$6+$7+$8+$9+$11+$12+$13+$14>=10' | Rscript ./source/BS-SNV-Caller.R | gzip > sample.bssnv.gz`

### parameters

tune the parameters in the code

```{R}
# mutation rate
pm = 1/1000/3

# error rate
# in GV oocyte (2n 4c) samples, error rate is set larger
# pe = 1/100/3
pe = 3/100/3

# total mis rate
p = pm + pe

# methylation rate/proportion 
pr.cg = 0.6        # CG context
pr.ncg = 1/100     # non-CG context
```

## Biologically consistent ASM CpG (BCA CpG) caller

**call BCA SNVs among a group of cells upon analysing allele-specific methylation (ASM) in each single cells**

usage

`source ./source/BCA-ASM-Caller.R`


**input files** :




SNV file: `merge_sample.bssnv.simple`

 a merged file containing detected SNV sites in each samples, four columns of chromosome, reference base, postion, and sample, no header line

```
1       T       1554    AF1_12
1       C       15023   AF1_12
1       G       15088   AF1_12
1       G       15091   AF1_12
1       T       15103   AF1_12
```


ASM file: `ASM_site_methpipe_allelicmeth_DP_merge.asm`

 a merged file containing detected ASM CpGs in each cells, 12 columns of 'chr', 'pos', 'strand', 'CpG', 'p_value', 'dp', 'MM', 'MU', 'UM', 'UU', and 'sample', no header line. The first 11 clomuns are returned by `methpipe::allelicmeth`


```
1       23736   +       CpG     1       11      11      0       0       0       AF1_12
1       50772   +       CpG     1       13      0       0       0       13      AF1_12
1       50809   +       CpG     1       13      0       0       0       13      AF1_12
1       69991   +       CpG     1       12      0       0       0       12      AF1_12
1       99975   +       CpG     1       16      16      0       0       0       AF1_12
1       100868  +       CpG     1       16      16      0       0       0       AF1_12
```


**output file:**

`./output/asm.consensus.tsv`

