# methyConcerto

Refined data analyses for **single-cell bisulfite-based sequencing** in DNA methylation researches, including

- Empirical Bayesian Methylation Caller (EBMC): Determinate discrete methyation status for each cytosine
- BS-SNV-Caller: SNV calling with bisuifite-converted sequencing data
- Biologically consistent ASM CpG (BCA CpG) caller: consistent allele-specific methylation (ASM) among a group of samples (cells), rather than occuring in individual cells occasionally




## Empirical Bayesian Methylation Caller (EBMC)

Determinate discrete methyation status for each cytosine from C/T read counts. The methylation state of a cytosine is discrete: methylated (1), unmethylated (0), and hemimethyalted (0.5). EBMC gives the predicted state with corresponding posterior probability.

### usage

process CG and nonCG sites seperately

`zcat ./data/sample.CGmap.gz | cut -f 7,8 | Rscript ./source/EBMC.R --output ./output/sample.prob --summary ./output/sample.methycall --CG TRUE`

### parameters

|param | discription|
|  ----  | ----  |
|--CG  |boolean, TRUE or FALSE, CG site or not |
|--prior| prior probabilities, defaults: `(0.59, 0.4, 0.01)` for CG and `(0.01, 0.98, 0.01)` for nonCG |
|--error_rate| prior rate of C/T mis-detction, default: 0.01|
|--Poisson_approx| threshold of depth the Poisson approximation applied, default: 20|
|--rounds|rounds of iteration, default: 2|


### input 

`.CGmap` file for a sample/cell:

 ./data/sample.CGmap.gz, returned by `baseeker2/3` or `cgmaptools`

### output

`.prob` probability file with same \#rows of input: two columns

|column | discription|
|  ----  | ----  |
|1| predicted methylation state with maximum posterior probability|
|2| posterior probability of the predicted state|

```
0       0.99911
0       0.990017
0       0.996467
0       0.996467
1       0.99932
0.5     0.800168
0.5     0.800168
0       0.999775
0       0.999775
0       0.996467
0       0.990017
1       0.999995
1       0.999995
```


`.methycall` summary file

|var | discription|
|  ----  | ----  |
|Rows|number of input rows/sites|
|methylated(1)| proportion of methylated cytosines  |
|unmethylated(0)| proportion of unmethylated cytosines  |
|epi-heterozygous(0.5)| proportion of hemimethylated cytosines  |
|error rate of measurement (0->1 or 1->0)| error rate of unexpected reads, errors due to libraray construction, sequencing, alignment etc.|


```
Rows:	10000
methylated(1):	0.188000
unmethylated(0):	0.807900
epi-heterozygous(0.5):	0.004100
error rate of measurement (0->1 or 1->0):	0.002497
```


## BS-SNV-Caller

SNV calling with bisuifite-converted sequencing data

### usage

`zcat ./data/sample.ATCGmap.gz | awk '$6+$7+$8+$9+$11+$12+$13+$14>=10' | Rscript ./source/BS-SNV-Caller.R | gzip > ./output/sample.bssnv.gz`

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

### input

 a `.ATCG` file, returned by `bsseeker2/3` or `cgmaptools`

### output

13 columns, no header line

|column | discription|
|  ----  | ----  |
|1-5| chromosome, reference base, position, CG/CHG/CHH, dinucleoties. They are same with first 5 columns of `.ATCG` file|
|6| ***p*-value** testing SNV (different from reference) which means the posterior probability of reference base, the smaller, the more likely to be a SNV, homozygous or heterzygous |
|7-10| estimated allele frequencies of A,T,C, and G respectively |
|11-12| coverage depths of Watson and Crick strands|
|13| probability of homozygote|

```
1	C	1550	CHG	CT	9.79052e-01	7.46622e-07	1.05278e-02	9.89471e-01	7.46622e-07	14	0	0.97916
1	T	1551	--	--	9.78959e-01	7.46623e-07	9.89423e-01	1.05750e-02	7.46623e-07	14	0	0.97907
1	G	1552	CHG	CA	9.99995e-01	7.55074e-07	7.55074e-07	7.55074e-07	9.99998e-01	14	0	1.00000
1	C	1553	CHG	CT	9.79052e-01	7.46622e-07	1.05278e-02	9.89471e-01	7.46622e-07	14	0	0.97916
1	T	1554	--	--	2.13408e-24	3.49345e-05	6.55191e-03	3.49345e-05	9.93378e-01	14	0	0.98676
1	G	1555	CHG	CA	9.99995e-01	7.55074e-07	7.55074e-07	7.55074e-07	9.99998e-01	14	0	1.00000
1	C	1556	CHH	CT	9.79052e-01	7.46622e-07	1.05278e-02	9.89471e-01	7.46622e-07	14	0	0.97916
1	T	1557	--	--	9.78959e-01	7.46623e-07	9.89423e-01	1.05750e-02	7.46623e-07	14	0	0.97907
1	C	1558	CHH	CT	9.79052e-01	7.46622e-07	1.05278e-02	9.89471e-01	7.46622e-07	14	0	0.97916
1	T	1559	--	--	9.78959e-01	7.46623e-07	9.89423e-01	1.05750e-02	7.46623e-07	14	0	0.97907
```

## Biologically consistent ASM CpG (BCA CpG) caller

Call BCA allele-specific methylation (ASM) among a group of cells upon the ASM analysis in each single cells. The BCA ASM is ASM consistent among a group of samples (cells). The *p*-value is evaluated under the null hypothesis that ASM CpGs are distributed randomly in the whole genome and independently among cells.


### usage

`source ./source/BCA-ASM-Caller.R`

### input files

SNV file: `./data/snv-10K-merge.simple.gz`

 a merged file containing detected SNV sites in each samples, four columns of chromosome, reference base, postion, and sample. no header line

```
1       T       1554    AF1_12
1       C       15023   AF1_12
1       G       15088   AF1_12
1       G       15091   AF1_12
1       T       15103   AF1_12
```


ASM file: `./data/methpipe-10K-merge.asm.gz`

a merged file containing detected ASM CpGs in each cells, 12 columns of 'chr', 'pos', 'strand', 'CpG', 'p_value', 'dp', 'MM', 'MU', 'UM', 'UU', and 'sample'. no header line. The first 11 clomuns are returned by `methpipe::allelicmeth`


```
1       23736   +       CpG     1       11      11      0       0       0       AF1_12
1       50772   +       CpG     1       13      0       0       0       13      AF1_12
1       50809   +       CpG     1       13      0       0       0       13      AF1_12
1       69991   +       CpG     1       12      0       0       0       12      AF1_12
1       99975   +       CpG     1       16      16      0       0       0       AF1_12
1       100868  +       CpG     1       16      16      0       0       0       AF1_12
```


### output file 

`./output/asm.consensus.tsv`


|variable | discription|
|  ----  | ----  |
|chr| chromosome|
|pos| position|
|count.cov | count of cells with this CpG effectively covered|
|count.asm| count of cells indicating ASM in this CpG|
|chi.sq| Chi-square ($\chi^2$) statistic|
|p.chisq| *p*-value of $\chi^2$ test, minor method|
|lambda| lambda ($\lambda$) statistic|
|p.pois|***p*-value of Poisson test, major method**|


```
chr	pos	count.cov	count.asm	count.snv	chi.sq	p.chisq	lambda	p.pois
16	8868	48	47	2	NA	NA	26.758615249594396	4.508674530477948e-4
16	9276	49	47	4	Inf	0	26.067905464137255	8.393246565554962e-4
16	9284	50	47	0	Inf	0	29.110167311050276	8.200382170259079e-4
16	9507	50	47	1	Inf	0	28.493378989882395	9.142007577857896e-4
8	1661	48	47	0	Inf	0	28.361876688204955	4.803623038606877e-4
```

## Citations

- Xiaolong Yuan, Na Chen, Yance Feng, et al., Single-cell multi-omics profiling reveals key regulatory mechanisms that poise germinal vesicle oocytes for meiotic resumption in pigs. *Cellular and Molecular Life Sciences*, 2023, under review.

