# duckweed-spread-comp

# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code
CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication

## 2. Data derived from other sources

NA

## 3. Recommended citation for this data/code archive



# DATA & CODE FILE OVERVIEW

This data repository consist of 9 data files, 2 code scripts, and this README document, with the following data and code filenames and variables

## Data files and variables

See col_names.md in data folder

## Code scripts and workflow

1. plots_and_analyses.R - generates all statistical analyses and figures in "Eco-evolutionary dynamics are shaped by competition in experimental range expansions". To run, download whole repo and set working directory to parent directory (i.e. download location).

2. simulations.R - demonstrates how simulated confidence intervals were generated. Note this is for review purposes, all data files have relevant columns added in advance, so this file does not need to be run in order to reproduce results.

# SOFTWARE VERSIONS

All analyses were conducted in R v. 4.5.1. Required packages are:
1. dplyr 2.5.0
2. permute 0.9-7
3. readxl 1.4.5
4. tidyR 1.3.1
5. vegan 2.7-1
6. boot 1.3-32
