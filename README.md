# SNPPar_test
 Data and scripts for testing SNPPar with simulated sequence data


1. Explanation of test data

1.1 3 levels of data 

Viet Nam (Ho Chi Minh City) lineage 2 - 'HCMC_L2'
Global lineage 2 - 'Global_L2'
Global lineages 1, 2 and 4 - 'Global_L124'

1.2 4 data set sizes each

HCMC_L2 and Global_L2 (n=821 and n=940 respectively)

Random samples of 10, 20 and 50 percent of isolates, plus all isolates

	Size	Labels	HCMC_L2 	Global_L2
	10% 	'_r10p'	n=82    	n=94 
	20% 	'_r20p'	n=164   	n=188 
	50% 	'_r50p'	n=410   	n=470 
	100%	'_all'	n=821   	n=940

Global_L124 (n=2965)
Random samples of (n=)100, 500, 1000 and 2000 isolates

'_r100' 
'_r500'
'_r1000' 
'_r2000'


Examples of each step using 'HCMC_L2_r10p' are given below.

2. Obtaining random samples

Generated from one of three different lists;

'HCMC_L2.txt', 'Global_L2.txt', and 'Global_L124.txt'.

Note: Most scripts use python 3, often along with the packages biopython and ete3 (same as SNPPar), unless otherwise noted.

Command:

python getRandomSet.py -l HCMC_L2.txt -n 82 -o HCMC_L2_r10p.txt

3. Obtaining initial SNP tables

Note: 'NC_000962_3_NC_000962_alleles_var.csv.zip' needs to be unzipped before the next step 
Warning: The SNP table 'NC_000962_3_NC_000962_alleles_var.csv' is a 3 Gb file after decompression

Note: 'parseSNPtable.py' uses python 2 and not 3. It also requires ete2 and biopython

The three modules used in this run are 'filter', which uses 'TB_excluded_regions_with_repeat.txt' to remove SNPs that occur in previously identified repeat regions, 'cons'(erve), which removes those SNPs with more than (1-cons)*100*isolates missing calls, and 'aln' (alignment), which outputs the SNP table in MFASTA format for Tree generation using RAxML.

Command:

python parseSNPtable.py -s NC_000962_3_NC_000962_alleles_var.csv -x TB_excluded_regions_with_repeat.txt -r NC_00962_3_1.gbk -m filter,cons,aln -l HCMC_L2_r10p.txt -p HCMC_L2_r10p_ -c 0.95

Conservation Note:

					r10p, r20p, r100 	-c 0.95
					r50p, r500			-c 0.99
					all, r1000, r2000	-c 0.995

Outputs (of importance to us)

SNP table (CSV) - 'X' rename to 'HCMC_L2_r10p_alleles_82strains_var_rF_c95.csv'
MFASTA alignment - 'X' rename to 'HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta'

Resulting SNP tables for all twelve data sets are in 'SNP_tables.zip'. Resulting MFASTA files are in 'MFASTA.zip'.

4. Obtaining the 'best' Tree

Note: RAxML was run on a computational cluster (Massive LINK) to make use of the multiple-threaded version. Each of up to 16 threads were allocated 4 Gb. Each data set was run five times using different seed values for each run - these seeds were 'recycled' for each data set.

		Seeds		Run
		 		1			2			3			4			5
		'-p'	8876251475	1885216875	8282263475	1118374265	8144432745
		'-x'	31566		315646		3444666		15623245	553324366

Version used: 8.2.12

Command:

raxml -T 7 -s HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta -n HCMC_L2_r10p_t01.tre -f a -m ASC_GTRGAMMA --asc-corr=lewis -p 8876251475 -x 31566 -N 100

Information from each run was collected from the resulting 'RAxML_info' file (e.g. 'data/RAxML_info.HCMC_L2_r10p_t01.out') and collated in 'data/RAxML_runs.tsv'. Information collected include the final maximum likelihood Optimization Likelihood (OL), the alpha distribution of GAMMA estimated from the data set along with the individual mutation rates between base pairs (e.g. A<->T) relative to the rate of G<->T, and the base frequencies. The latter three were used as estimates for modelling of sequence evolution in SeqGen (see below), whilst the first value was used to pick the 'best' Tree from the five replicates (i.e. ml estimate closest to zero). A edited version of 'data/RAxML_runs.tsv' with only the statistics for the 'best' Tree was also produced ('data/best_RAxML_runs.tsv').

In each case, the tree used was the 'RAxML_bipartitions' version. The 'best' RAxML Tree for each data set was opened in FigTree (Link), where each was midpoint rooted and nodes set to 'descending' order. Each tree was then exported in Newick format. The folder 'data/trees' contains a reformatted 'best' Tree for each of the twelve data sets (eg. 'data/trees/HCMC_L2_r10p.tre').

5. Run SeqGen

6. Extract variable sites

7. Extract parallel mutation events (expected)

8. Run simulated data with SNPPar (observed)

9. Compare expected and observed results

10. Further statistical analysis
