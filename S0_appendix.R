################################################################################
# Script S.0: Appendix
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
# Here there are demonstrations for the subsections of the Appendix of the manuscript.
#
# 1. Heterosis is not the same as specific combining ability (SCA):
#    
#   
# 2. Genetic distances do not necessarily measure complementarity: 
#    
#    
# 3. Maximise testcross variance at one locus and GCA rankings:
#
#
# 
# Key outputs:
# - Calculation of within-pool diagnostic metrics (representativeness).
# - Calculation of between-pool diagnostic metrics (GCA-SCA alignment, 
#   candidate lines GCA variance, and per se GCA of testers).
# - Identification of biological trade-offs between profiles (e.g., demonstrating
#   how high-performing testers often mask candidate line variance).
#
################################################################################

Load data from S1 (Global Parameters) and S3 (Testcross Parameters)
load("./testersel/data/S1_data.RData")
load("./testersel/data/S3_testcrosses_data.RData")

