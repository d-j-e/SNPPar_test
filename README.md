# SNPPar_test
 Data and scripts for testing [SNPPar](https://github.com/d-j-e/SNPPar) with either simulated or empirical datasets

# 1. Explanation of test data

## 1.1 Three levels of data 

* Vietnam (Ho Chi Minh City) lineage 2 - 'HCMC_L2'
* Global lineage 2 - 'Global_L2'
* Global lineages 1, 2 and 4 - 'Global_L124'

## 1.2 Four dataset sizes each

HCMC_L2 and Global_L2 (n=821 and n=940 respectively)

Random samples of 10, 20 and 50 percent of isolates, plus all isolates

	Size	Labels	HCMC_L2	Global_L2
	10%	'_r10p'	82	94
	20%	'_r20p'	164	188 
	50%	'_r50p'	410	470 
	100%	'_all'	821	940

Global_L124 (n=2965)
Random samples of (n=)100, 500, 1000 and 2000 isolates

	size	label
	100	'_r100'
	500	'_r500'
	1000	'_r1000'
	2000	'_r2000'


Examples of each step using 'HCMC_L2_r10p' are given below. The relevant data will be found in the first folder for the step each is introduced.

# 2. Obtaining random samples

Generated from one of three different lists;

'step2/HCMC_L2.txt', 'step2/Global_L2.txt', and 'step2/Global_L124.txt'.

Note: All of the python scripts provided use python 3, often along with the packages biopython and ete3 (same as SNPPar).

### Command:

python getRandomSet.py -l step2/HCMC_L2.txt -n 82 -o step2/HCMC_L2_r10p.txt

The randomly generated lists used for testing are provided (as I forgot a 'seed') in the folder 'step2'

# 3. Obtaining initial SNP tables

Note: 'step3/NC_000962_alleles_3786strains_var.csv.zip' needs to be unzipped before the next step 
Warning: The SNP table 'NC_000962_alleles_3786strains_var.csv' is a 900 MB file after decompression

The three modules used in this run are 'filter', which uses 'TB_excluded_regions_with_repeat.txt' to remove SNPs that occur in previously identified repeat regions, 'cons'(erve), which removes those SNPs with more than (1-cons)\*100\*isolates missing calls, and 'aln' (alignment), which outputs the SNP table in MFASTA format for Tree generation using RAxML.

### Command:

python parseSNPtable.py -s NC_000962_alleles_3786strains_var.csv -x step3/TB_excluded_regions_with_repeat.txt -r step3/NC_00962_3_1.gbk -m filter,cons,aln -l step2/HCMC_L2_r10p.txt -p HCMC_L2_r10p_ -c 0.95

Conservation Note:

	dataset			conservation level
	r10p, r20p, r100	-c 0.95
	r50p, r500		-c 0.99
	all, r1000, r2000	-c 0.995

Outputs (of importance to us)

SNP table (CSV) - 'X' rename to 'HCMC_L2_r10p_alleles_82strains_var_rF_c95.csv'
MFASTA alignment - 'X' rename to 'HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta'

Resulting SNP tables for all twelve data sets are in 'step3/SNP_tables.zip'. Resulting MFASTA files are in 'step3/MFASTA.zip'.

# 4. Obtaining the 'best' Tree

Note: RAxML was run on a computational cluster (Massive LINK) to make use of the multiple-threaded version. Each of up to 16 threads were allocated 4 Gb. Each data set was run five times using different seed values for each run - these seeds were 'recycled' for each data set.

		Seeds	Run
		 	1		2		3		4		5
		'-p'	8876251475	1885216875	8282263475	1118374265	8144432745
		'-x'	31566		315646		3444666		15623245	553324366

### Version used: 8.2.12 (multi-threaded)

### Command:

raxml -T 7 -s HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta -n HCMC_L2_r10p_t01.tre -f a -m ASC_GTRGAMMA --asc-corr=lewis -p 8876251475 -x 31566 -N 100

Information from each run was collected from the resulting 'RAxML_info' file (e.g. 'set4/RAxML_info.HCMC_L2_r10p_t01.out') and collated in 'set4/RAxML_runs.tsv'. Information collected include the final maximum likelihood Optimization Likelihood (OL), the alpha distribution of GAMMA estimated from the data set along with the individual mutation rates between base pairs (e.g. A<->T) relative to the rate of G<->T, and the base frequencies. The latter three were used as estimates for modelling of sequence evolution in SeqGen (see below), whilst the first value was used to pick the 'best' Tree from the five replicates (i.e. ml estimate closest to zero). A edited version of 'data/RAxML_runs.tsv' with only the statistics for the 'best' Tree was also produced ('set4/best_RAxML_runs.tsv').

In each case, the tree used was the 'RAxML_bipartitions' version. The 'best' RAxML Tree for each data set was opened in FigTree (Link), where each was midpoint rooted and nodes set to 'descending' order. Each tree was then exported in Newick format. The folder 'data/trees' contains a reformatted 'best' Tree for each of the twelve data sets (eg. 'set4/trees/HCMC_L2_r10p.tre').

# 5. Run SeqGen

# 6. Extract variable sites

# 7. Extract parallel mutation events (expected)

# 8. Run simulated data with SNPPar (observed)

# 9. Compare expected and observed results

# 10. Further statistical analysis

# 11. Real Datasets

## 11.1 _Elizabethkingia anophelis_

## 11.2 _Burkholderia dolosa_

## 11.3 _Mycobacterium tuberculosis_

