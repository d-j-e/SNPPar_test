# SNPPar_test
 Data and scripts for testing [SNPPar](https://github.com/d-j-e/SNPPar) with either simulated or empirical datasets

# 1. Explanation of simulated data

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

'Step2/HCMC_L2.txt', 'Step2/Global_L2.txt', and 'Step2/Global_L124.txt'.

Note: All of the python scripts provided use python 3, often along with the packages biopython and ete3 (same as SNPPar).

### Command:

python getRandomSet.py -l HCMC_L2.txt -n 82 -o HCMC_L2_r10p.txt

The randomly generated lists used for testing are provided (as I forgot a 'seed') in the folder 'Step2'

# 3. Obtaining initial SNP tables

Note: 'data/Step3/NC_000962_alleles_3786strains_var.csv.zip' needs to be unzipped before the next step 
Warning: The SNP table 'NC_000962_alleles_3786strains_var.csv' is a 900 MB file after decompression

The three modules used in this run are 'filter', which uses 'TB_excluded_regions_with_repeat.txt' to remove SNPs that occur in previously identified repeat regions and PE/PPE genes, 'cons'(erve), which removes those SNPs with more than (1-cons)\*100\*isolates missing calls, and 'aln' (alignment), which outputs the SNP table in MFASTA format for Tree generation using RAxML.

### Command:

python parseSNPtable.py -s NC_000962_alleles_3786strains_var.csv -x TB_excluded_regions_with_repeat.txt -r NC_00962_3_1.gbk -m filter,cons,aln -l HCMC_L2_r10p.txt -p HCMC_L2_r10p_ -c 0.95

**Conservation Note:**

	dataset			conservation level
	r10p, r20p, r100	-c 0.95
	r50p, r500		-c 0.99
	all, r1000, r2000	-c 0.995

Outputs (of importance to us)

SNP table (CSV) - 'HCMC_L2_r10p_alleles_82strains_var_rF_c95.csv'
MFASTA alignment - HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta'
(Note: these have been renamed to reduce the name lengths)

Resulting SNP tables for all twelve data sets are in 'dat/Step3/SNP_tables.zip'. Resulting MFASTA files are in 'data/Step3/MFASTA.zip'.

# 4. Obtaining the 'best' tree using RAxML

Note: [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) was run on a computational cluster ([MASSIVE](https://www.monash.edu/research/infrastructure/platforms-pages/massive)) to make use of the multiple-threaded version. Each of up to 16 threads were allocated 4 Gb. Each data set was run five times using different seed values for each run - these seeds were 'recycled' for each data set.

	Seeds	Run
	 	1		2		3		4		5
	'-p'	8876251475	1885216875	8282263475	1118374265	8144432745
	'-x'	31566		315646		3444666		15623245	553324366

### Version used: 8.2.12 (multi-threaded)

### Command:

raxml -T 16 -s HCMC_L2_r10p_alleles_82strains_var_rF_c95.mfasta -n HCMC_L2_r10p_t01.tre -f a -m ASC_GTRGAMMA --asc-corr=lewis -p 8876251475 -x 31566 -N 100

Information from each run was collected from the resulting 'RAxML_info' file (e.g. 'set4/RAxML_info.HCMC_L2_r10p_t01.out') and collated in 'set4/RAxML_runs.tsv'. Information collected include the final maximum likelihood Optimization Likelihood (OL), the alpha distribution of GAMMA estimated from the data set along with the individual mutation rates between base pairs (e.g. A<->T) relative to the rate of G<->T, and the base frequencies. The latter three were used as estimates for modelling of sequence evolution in SeqGen (see below), whilst the first value was used to pick the 'best' tree from the five replicates (i.e. ML estimate closest to zero). A edited version of 'data/RAxML_runs.tsv' with only the statistics for the 'best' Tree was also produced ('step4/best_RAxML_runs.tsv').

In each case, the tree used was the 'RAxML_bipartitions' version. The 'best' RAxML tree for each data set was opened in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), where each was midpoint rooted and nodes set to 'descending' order. Each tree was then exported in Newick format. The folder 'data/Step4/trees' contains a reformatted 'best' tree for each of the twelve data sets (eg. 'Step4/trees/HCMC_L2_r10p.tre').

# 5. Run SeqGen

Note: [SeqGen](http://tree.bio.ed.ac.uk/software/seqgen/) was also run on a computational cluster ([MASSIVE](https://www.monash.edu/research/infrastructure/platforms-pages/massive)) to allow for the output files (\~1.3 GB each simulation). The processing of these output files (Command 2) were also completed on the same cluster environment.

As the amount of simulated SNPs was found to be much higher than expected (due to a poor choice of ascertainment bias correction used in RAxML), for each of the 12 datatsets an internal correction was used where branch lengths were multiplied by the expected number of SNPs (the number of SNPs in the real dataset that created the tree) divided by the observed simulated SNPs in the first run of each of the twelve populations/samplesize combinations ('E/O', see table below for values used here).

Note: as the seqgen output files are so large, these have been excluded from the repository.

The variable sites for each replicate were extracted using processSeqGen.py in the 'scripts' folder, and separated into those calls for internal nodes and tips. It produces the tips as both FASTA and SNP table output.

### Command 1, run 1:

seqgen -m GTR -a [alpha] -f [#A],[#C],[#G],[#T] -r [A<>C],[A<>G],[A<>T],[C<>G],[C<>T],[G<>T] -of -wa -k1 [treefile]

where values for -a -f and -r are found in 'set4/best_RAxML_runs.tsv'
and treefile is the 'Phylip' format input for SeqGen containing the seed sequence (H37Rv, shortened to account for filtered sites). These can be found in 'data/Step5/phy_files'. The shortened sequence length is 4,022,248 bp (4,411,532 total genome length - 389,284 excluded bases).

For HCMC_L2_r10p pre-run 1, the SeqGen command was:

seqgen -m GTR -a 2.76 -f 0.15,0.35,0.34,0.16 -r 1.06,3.29,0.35,0.54,3.11,1 -of -wa -k1 mtb_short_seq_HCMC_L2_r10p.phy

### Command 1, run 2 

(note: number of SNPs established by running processSeqGen.py - see command 2):

seqgen -m GTR -a [alpha] -f [#A],[#C],[#G],[#T] -r [A<>C],[A<>G],[A<>T],[C<>G],[C<>T],[G<>T] -of -wa -k1 [treefile] -t [E/O]

For HCMC_L2_r10p run 1, the updated SeqGen command was:

seqgen -m GTR -a 2.76 -f 0.15,0.35,0.34,0.16 -r 1.06,3.29,0.35,0.54,3.11,1 -of -wa -k1 mtb_short_seq_HCMC_L2_r10p.phy **-t 0.044**

### Command 2:

python processSeqGen.py -i [seqgen_output] -p [output_prefix]

For HCMC_L2_r10p run 1:
python processSeqGen.py -i [seqgen_output_example] -p HCMC_L2_r10p_r1

**E/O corrections**

    replicate 		expected	observed	E/O
    HCMC_L2_r10p		5561		125150		0.044
    HCMC_L2_r20p		9027		188805		0.048
    HCMC_L2_r50p		16824		219728		0.077
    HCMC_L2_all		25547		289022		0.088
    Global_L2_r10p		4303		161383		0.027
    Global_L2_r20p		7956		201160		0.04
    Global_L2_r50p		17971		250597		0.052
    Global_L2_all		20014		256046		0.078
    Global_L124_r100	12911		139064		0.093
    Global_L124_r500	29442		213650		0.138
    Global_L124_r1000	43038		281296		0.153
    Global_L124_r2000	63451		398169		0.159

In the subfolder 'data/Step5/simulated/HCMC_L2)r10p' is the folder with the results for run 1 to 10 of HCMC_L2_r10p (*sans* SeqGen output as discussed above). The rest of the outputs for the 11 other datasets are either zipped into a single file (e.g. 'HCMC_L2_r20p.zip' in 'Step5/simulated') or split into separate zipped files (<100 MB for GitHub, e.g. 'Global_L124_r1000_r1_to_r3.zip' in 'Step5/simulated/Global_L124_r1000').

# 6. Extract expected mutation events and homoplasies

### Command:

# 7. Run simulated data with SNPPar (observed)

### Command:

# 8. Compare expected and observed results

### Command:

# 9. Further statistical analysis

### R-markdown script

# 10. Real Datasets

## 10.1 _Elizabethkingia anophelis_

### Reference genome:

### SNPtable:

### Tree:

### Command:

### Post-run analysis:

## 10.2 _Burkholderia dolosa_

### Reference genomes:

### SNPtable:

### Tree:

### Commands:

### Post-run analysis:

## 10.3 _Mycobacterium tuberculosis_

### Reference genome:

### SNPtable:

### Tree:

### Command:

### Post-run analysis:

### R-markdown script:
