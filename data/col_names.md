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
GenNo\
&ensp; The number of distinct genotypes present in the sample\
InvSimpson\
&ensp; The inverse of the Simpson's diversity index\
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
&ensp; upper and lower bounds of 95% confidence intervals of differences in low density population growth rate from simulated pairwise samples of size n with equal genotype frequencies
#### Columns only in data_summaryedges
LJ01_diff -- LM06_diff\
&ensp; these columns give the difference in frequency between edge and core samples in a given replicate\
LJ01_sigdiff -- LM06_sigdiff\
&ensp; A binary variable for whether the difference between edge and core genotype frequency sample falls outside the 95% confidence estimate for pairwise comparisons with equal genotype frequencies\
SLA_diff\
&ensp; The difference in specific leaf area z-score between edge and core samples for the given replicate.\
root_diff\
&ensp; The difference in root length z-score between edge and core samples for the given replicate.\
raft_diff\
&ensp; The difference in raft size z-score between edge and core samples for the given replicate.\
growth_diff\
&ensp; The difference in low density population growth rate z-score between edge and core samples for the given replicate.\
div_diff\
&ensp; The difference in Inverse Simpson's  between edge and core samples for the given replicate.
### Extent
extent-over-time_clean.csv

Columns L-R\
Treatment\
&ensp; C - control i.e. -_Spirodela_\
&ensp; S - Spirodela i.e. +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
BenchID\
&ensp; ID of the greenhouse bench replicate was positioned on, range is 1-5\
day\
&ensp; Day of experiment (0-30)\
back\
&ensp; Furthest back occupied cm of the gutter\
front\
&ensp; Furthest forward occupied cm of the gutter\
gap_back\
&ensp; Furthest back cm of any large gap within the expansion extent; NA otherwise\
gap_front\
&ensp; Furthest forward cm of any large gap within the expansion extent; NA otherwise\
algae\
&ensp; Level of algae present in the gutter, on a visual scale from 0 (none) to 3 (most water surface covered)\
extent\
&ensp; expansion extent i.e. difference between front and back.\
### Genotype data
genotypes_clean.csv

Columns L-R\
Treatment
&ensp; CTRL - -_Spirodela_\
&ensp; SPIRO - +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
Position\
&ensp; CORE - sample taken from population core\
&ensp; EDGE - sample taken from population edge\
Table\
&ensp; Greenhouse bench where given replicate was located (1-5)\
SampleID\
&ensp;Unique ID given to DNA sample\
GenotypeID\
&ensp; Genotype determined from microsatellite analysis
### Density data
density_data_clean.csv

Columns L-R\
cam\
&ensp; camera image was taken with\
img_id\
&ensp; file name of corresponding image\
treatment\
&ensp; C - control ie. -_Spirodela_\
&ensp; S - spirodela ie. +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
cm\
&ensp; centimeter of the gutter corresponding to the transect counted along\
gutter_width\
&ensp; width of the gutter at the water surface along the transect
frond_no\
&ensp; Number of fronds counted along the transect\
LE_cm\
&ensp; cm corresponding to the furthest forward transect counted along\
cm_plus\
&ensp; corrected distances so that highest cm value always corresponds to the leading edge, removes negative cm induced in expansions that spread from the far end of their gutter\
fromLE\
&ensp; The distance in cm of the given transect from the furthest forward transect counted along.
### Wave Steepness
 - steepness_clean.csv

Columns from L-R
treatment\
&ensp; C - control ie. -_Spirodela_\
&ensp; S - spirodela ie. +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
peak\
&ensp;Location (according to cm_plus column in densities) of the largest population density\
fifty_pct_edge\
&ensp;Location (according to cm_plus column in densities) of the furthest forward transect to reach at least 50% of the peak population density\
edge\
&ensp; Location of the furthest forward transect counted (according to cm_plus column in densities)\
peak_edge_dist\
&ensp; Difference between edge and peak\ 
fifty_pct_edge_dist\
&ensp; Difference between edge and fifty_pct_edge

- peak_edge_pop.csv\
  
Columns from L-R\
treatment\
&ensp; C - control ie. -_Spirodela_\
&ensp; S - spirodela ie. +_Spirodela_\
Rep\
&ensp; replicate number 1-25 for each treatment\
peak\
&ensp;Population size at the location (according to cm_plus column in densities) of the largest population density\
edge\
&ensp; Population size at the location of the furthest forward transect counted (according to cm_plus column in densities)

### Traits
traits_clean.csv\

columns L-R\
genotypeID\
&ensp; identifier for each unique genotype\
root_mean\
&ensp; mean of root length (mm) across measured individuals\
root_sd\
&ensp; standard deviation of root length\
root_z\
&ensp;z-score of root length compared to average across genotypes\
SLA_mean\
&ensp; mean of specific leaf area across measured individuals\
SLA_sd\
&ensp; standard deviation of specific leaf area\
SLA_z\
&ensp;z-score of specific leaf area compared to average across genotypes\
rafts_mean\
&ensp; mean of raft size (fronds per raft) across measured individuals\
rafts_sd\
&ensp; standard deviation of raft size\
rafts_z\
&ensp;z-score of raft size compared to average across genotypes\
growth_mean\
&ensp; mean number of fronds produced in 7 days by a single frond across measured individuals (i.e. low density population growth rate)\
growth_sd\
&ensp; standard deviation of low density population growth rate\
growth_z\
&ensp;z-score of low density population growth rate compared to average across genotypes\
