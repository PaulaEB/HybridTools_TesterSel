################################################################################
# Script S.5: Tester profiles
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
##
# Evaluates diagnostic metrics for optimising tester selection considering the 
# tester profiles in different breeding schemes based on long-term population 
# improvement vs. short-term product development. This script considers three 
# distinct profiles:
#
# 1. The Classic Tester: 
#    Maximises the selection accuracy of candidate lines by amplifying GCA 
#    variance (high positive GCA-SCA alignment, high b1).
# 2. The Parent Tester: 
#    Maximises per se GCA to produce high-performing commercial hybrids, 
#    even if it masks the variance of candidate lines (low or negative b1).
# 3. The Representative Tester: 
#    Captures the allele frequencies of the Target Pool of Genotypes (TPG), 
#    ensuring the development of reciprocal complementarity between pools.
# 
# Key outputs:
# - Calculation of all metrics (b1, candidates GCA variance, unaligned var, 
#   representativeness scores).
# - Visualization of Biological Trade-offs (Product vs Resolution/Accuracy).
# - Consolidated quantitative summary of tester profiles.
################################################################################

rm(list = ls())

source("./functions_TesterSel.R")
load("./data/S1_data.RData")

## ---------------------------------------------------------
## 1) Calculate key metrics for each of the tester profiles 
## ---------------------------------------------------------

# Between-pool metrics (GCA, SCA, b1, Noise, Candidate Variance)

# Candidate lines observed GCA variance per tester.
# Total variance of the candidate lines GCA is var_gca2
var_gca2_obs <- apply(H_centre, 2, pop_var)
 
# Testers own GCA per line 
gca1

# GCA-SCA alignment (b1) for each tester
cov_gca2_sca <- apply(sca, 2, function(sca_col) pop_cov(sca_col, gca2))
b1 <- cov_gca2_sca/var_gca2

# Within-pool metric (Representativeness)
x_bar_l <- colMeans(X1)

# Assign some weights

