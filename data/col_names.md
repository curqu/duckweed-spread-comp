# Key for data files and column names 

### Data summary 
data_summary_clean.csv
 - data_summarycores_clean.csv
 - data_summaryedges_clean.csv\
by replicate summary of genetic and trait data
parent file divied into edges and cores for easier processing -- column names are consistent except some only in edge file (noted below)

Columns L-R\
Treatment\
&ensp; C - control ie. -_Spirodela_\
&ensp; S - spirodela ie. +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
BenchID\
&ensp; ID of the greenhouse bench replicate was positioned on, range is 1-5\
final_ext\
&ensp; extent of the expansion on day 30\
position\
&ensp; CORE - summarizes the samples taken from the core of the population\
&ensp; EDGE - summarizes the samples taken from the edge of the population\
n\
&ensp; number of samples obtained\
prop_LJ01 -- prop_LM06\
&ensp; these columns give the relative proportion of each of 10 genotypes\
raft.z\
&ensp; z score of genotype-weighted raft size, ie. number of fronds per raft\
SLA.z\
&ensp; z score of genotype-weighted specific leaf area\
root.z\
&ensp; z score of genotype-weighted root length\
growth.z\
&ensp; z score of genotype-weighted low density growth rate\
LJ01_change -- LM06_change\
&ensp; these columns give the difference in frequency between day 30 and day 1\
GenNo\
&ensp; The number of distinct genotypes present in the sample\
Simpson\
&ensp; Simpson's diversity index\
InvSimpson\
&ensp; The inverse of the Simpson's diversity index\
Abundance\
&ensp; ???\
CI_freq\
&ensp; Upper bound of 95% confidence interval from simulated samples of size n with equal genotype frequencies\
CIup_diff / CIlow_diff\
&ensp; upper and lower bounds of 95% confidence intervals from simulated pairwise differences in genotype frequency for samples of size n with equal genotype frequencies\
sig_LJ01 -- sig_LM06\
&ensp; A binary variable for whether the sample falls outside the 95% confidence estimate for equal genotype frequencies\
CI_raftUP / CI_raftLOW\
&ensp; upper and lower bounds for 95% confidence intervals of raft size given equal genotype frequencies\
CI_SLAUP / CI_SLALOW\
&ensp; upper and lower bounds for 95% confidence intervals of specific leaf area given equal genotype frequencies\
CI_rootUP / CI_rootLOW\
&ensp; upper and lower bounds for 95% confidence intervals of root length given equal genotype frequencies\
CI_growthUP / CI_growthLOW\
&ensp; upper and lower bounds for 95% confidence intervals of low density populaton growth rate given equal genotype frequencies\
CIdiff_raftUP / CIdiff_raftLOW\
&ensp; upper and lower bounds of 95% confidence intervals of differences in raft size from simulated pairwise samples of size n with equal genotype frequencies\
CIdiff_SLAUP / CIdiff_SLALOW\
&ensp; upper and lower bounds of 95% confidence intervals of differences in specific leaf area from simulated pairwise samples of size n with equal genotype frequencies\
CIdiff_rootUP / CIdiff_rootLOW\
&ensp; upper and lower bounds of 95% confidence intervals of differences in root length from simulated pairwise samples of size n with equal genotype frequencies\
CIdiff_growthUP / CIdiff_growthLOW\
&ensp; upper and lower bounds of 95% confidence intervals of differences in low density population growth rate from simulated pairwise samples of size n with equal genotype frequencies\
