################################################################################
# Script S.3: Heterogeneous variance in combining ability components
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
#
# Demonstration of the heterogeneous variance in combining ability components and
# its impact on the accuracy of GCA estimation. In this context, SCA is not just 
# uniform noise; if aligned to GCA, it amplifies the GCA discrimination signal. 
# This script simulates a testcross with a single-tester
#
# Assumptions:
# - GCA and SCA are not independent (non-zero covariances exist).
# - SCA variance is heteroscedastic across testers.
# 
# Key outputs:
# - Demonstration of variable SCA variance and covariance(GCA, SCA) != 0.
# - Comparison of expected accuracy under homogeneous vs heterogeneous assumptions.
# - Identification of testers that amplify (High b1), mask (Low b1), or 
#   contribute no ranking information (Zero b1) to GCA estimation.
#
################################################################################
rm(list = ls())

source("./functions_TesterSel.R")
load("./data/S1_data.RData")

## ---------------------------------------------------------------
## 1) Demonstrate Heterogeneous Variance and Non-Independence
## ---------------------------------------------------------------

# Define population variance and covariance functions 
pop_var <- function(x) { mean((x - mean(x))^2) }
pop_cov <- function(x, y) { 
  if (is.matrix(x)) {
    apply(x, 2, function(col) mean((col - mean(col)) * (y - mean(y))))
  } else {
    mean((x - mean(x)) * (y - mean(y)))
  }
}

# Obtain the SCA variance in crosses with each tester. This demonstrates that
# SCA variance is heteroscedastic across testers (sigma^2_s1j != sigma^2_s1k)
var_sca_testers <- apply(sca, 2, pop_var)
barplot(sort(var_sca_testers), las=2, ylab = "SCA Variance")

# Demonstrate that GCA of the candidate lines and SCA has a covariance different
# from 0. This breaks the standard independent model assumption
cov_gca2_sca <- apply(sca, 2, function(sca_col) pop_cov(sca_col, gca2))
hist(cov_gca2_sca); abline(v = 0, lty=2, lwd=3)

## -----------------------------------------------------------------------------
## 2) Expressing GCA-SCA Alignment as a Linear Regression Scope (b1)
## -----------------------------------------------------------------------------

# Knowing that the covariance between GCA and SCA is different to 0, let's see  
# how the correlation varies. This demonstrates that there is variation in this  
# correlation reflecting the heterogeneity in the SCA expected from the crosses 
# with different testers
r_gca2_sca <- cor(gca2,sca)
hist(r_gca2_sca); abline(v = 0, lty=2, lwd=3)

# Express this in linear regression terms where a higher slope (b1) reflects a 
# higher positive alignment of the SCA with the candidate lines GCA effects
b1 <- cov_gca2_sca/var_gca2
hist(b1); abline(v = 0, lty=2, lwd=3)

# See which testers have the highest alignment (amplifiers) and which have the 
# lowest alignment (maskers)
barplot(sort(drop(b1)), ylim = c(-0.6, 0.6), las = 2, ylab = "Alignment GCA-SCA (b1)")

## -----------------------------------------------------------------------------
## 3) Accuracy of GCA Selection Driven by the Alignment (r)
## -----------------------------------------------------------------------------

# Observed GCA Selection Accuracy defined as the correlation between the 
# true candidate lines GCA and the observed testcross performance
r_gca_obs <- cor(H_centre, gca2)

# Expected accuracy assuming heterogeneous variance (allowing covariances > 0)
b1_vec <- drop(b1)
names(r_gca_obs) <- names(b1)

aligned_var <- ((1 + b1)^2) * as.vector(var_gca2)
unaligned_var <- var_sca_testers - (b1^2 * as.vector(var_gca2)) #error variance see report 
total_var <- ((1 + 2 * b1) * as.vector(var_gca2)) + var_sca_testers

# Check math:
all.equal((aligned_var + unaligned_var), total_var) #TRUE
r_gca_exp_hetg <- sqrt(aligned_var/total_var)

# Expected accuracy assuming homogeneous variance (ignoring b1)
r_gca_exp_homg <- sqrt(as.vector(var_gca2) / (as.vector(var_gca2) + var_sca_testers/1))

# Plot check: The heterogeneous model explains observed accuracy better 
plot(r_gca_obs, r_gca_exp_homg, 
     pch = 19, col = adjustcolor("black", 0.5),
     xlab = "Observed r (testcross performance, true GCA)", 
     ylab = "Expected Accuracy",
     xlim = c(0.3,1), ylim = c(0.3,1))

points(r_gca_obs, r_gca_exp_hetg, pch = 19, col = "blue")
lines(x = c(0.3, 1), y = c(0.3, 1), lwd = 1, col = "grey30")
legend("topleft", legend = c("Heterogeneous variance", "Homogeneous variance"), 
       col = c("blue", "black"), pch = 19, bty = "n", cex = 0.8)

# Check if the observed variance matches the expected variance under heterogeneous
# assumption 
all.equal(unname(apply(H_centre, 2, pop_var)), unname(total_var)) #TRUE

## -----------------------------------------------------------------------------
## 4) Accuracy driven by the slope (b1) 
## -----------------------------------------------------------------------------

# After these results a sensible hypothesis is that testers with higher slope (b1)
# deliver more accuracy. Then, we would want to identify distinct testers based on b1:

# Highest b1 = high alignment GCA-SCA (amplifies GCA)
tester_high_b1 <- names(which.max(drop(b1)))  
# Lowest (negative) b1 = low alignment GCA-SCA (masks GCA)
tester_low_b1 <- names(which.min(drop(b1))) 
# Close to independent assumption (SCA as noise)
tester_zero_b1 <- names(which.min(abs(drop(b1))))  

# Calculate the observed variance of the candidate lines when crossed those single
# testers
var_tc_high <- pop_var(H_centre[, tester_high_b1])
var_tc_low  <- pop_var(H_centre[, tester_low_b1])
var_tc_zero <- pop_var(H_centre[, tester_zero_b1])

# Calculate the expected testcross variance of candidate lines GCA for these three
# testers using the Eq. 5: Var(g_1) + Var(s_1j) + 2*Cov(g_1, s_1j) --REPORT-- 
# because we allow covariances to exist
var_exp_high <- var_gca2 + pop_var(sca[, tester_high_b1]) + 2 * pop_cov(gca2, sca[, tester_high_b1])
var_exp_low  <- var_gca2 + pop_var(sca[, tester_low_b1])  + 2 * pop_cov(gca2, sca[, tester_low_b1])
var_exp_zero <- var_gca2 + pop_var(sca[, tester_zero_b1]) + 2 * pop_cov(gca2, sca[, tester_zero_b1])
var_homog_exp <- as.vector(var_gca2) + (var_sca/1)

# Check if the expected variance matches the variance calculated under the 
# heterogeneous variance assumption (Eq. 5) for the three testers
all.equal(as.numeric(total_var[tester_high_b1]), as.numeric(var_exp_high)) #TRUE
all.equal(as.numeric(total_var[tester_low_b1]),  as.numeric(var_exp_low)) #TRUE
all.equal(as.numeric(total_var[tester_zero_b1]), as.numeric(var_exp_zero)) #TRUE

# Barplot comparison showing that different testers render different candidate 
# lines GCA variances where the tester with a slope of 0 approximates to the 
# homogeneous variance assumption, the tester with low slope has the lower GCA
# variance and the one with the highest slope has the highest GCA variance
barplot_var <- c("Exp. Homog"  = var_homog_exp,
                 "Exp (Low b1)" = var_exp_low,
                 "TC (Low b1)"  = var_tc_low,
                 "Exp (Zero b1)" = var_exp_zero,
                 "TC (Zero b1)" = var_tc_zero,
                 "Exp (High b1)" = var_exp_high,
                 "TC (High b1)" = var_tc_high)
barplot(barplot_var, ylab = "GCA candidate lines variance", las = 1, 
        col = c("black","grey", "darkred", "grey","royalblue","grey", "darkgreen"))
abline(h = var_gca2, col = "black", lwd =1, lty = 2)

# Plot to see how the tester's slope (b1) aligns with the observed accuracy 
plot(b1, r_gca_obs, 
     pch = 19, col = "black", cex=0.5, ylim = c(0.4,1), xlim = c(-0.6,0.6),
     xlab = "Tester Slope (b1)", 
     ylab = "Observed r (Hybrid, true GCA)")

abline(lm(r_gca_obs ~ b1), col = "grey60", lty=2) 
points(b1[tester_high_b1], r_gca_obs[tester_high_b1], col="darkgreen", pch=8, cex=1, lwd=2)
points(b1[tester_zero_b1], r_gca_obs[tester_zero_b1], col="grey40", pch=8, cex=1, lwd=2)
points(b1[tester_low_b1],  r_gca_obs[tester_low_b1], col="darkred", pch=8, cex=1, lwd=2)

legend("bottomright", legend = c("High b1", "Zero b1", "Low b1"), 
       col = c("darkgreen", "grey40", "darkred"), pch = 8, bty = "n", cex = 0.8)

## -----------------------------------------------------------------------------
## 4) The noise of the tester: Unaligned variance
## -----------------------------------------------------------------------------
# A tester should amplify the GCA of the candidate lines (High b1) with the less  
# noise (unaligned variance) as possible. One may think that simply picking a tester 
# with low SCA variance would avoid noise, but that drops b1 to zero. 

# Let's see which tester has the lowest unaligned variance (lowest noise) and 
# how that relates to the observed accuracy.
tester_low_e <- names(which.min(abs(drop(unaligned_var)))) 
plot(unaligned_var, r_gca_obs, 
     pch = 21, bg = ifelse(b1 > 0, "black", "grey70"), col = "white", cex = 1,
     ylim = c(0.4,1),
     xlab = "Unaligned Variance", 
     ylab = "Observed r (Hybrid, true GCA)")

points(unaligned_var[tester_high_b1], r_gca_obs[tester_high_b1], col="darkgreen", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_zero_b1], r_gca_obs[tester_zero_b1], col="royalblue", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_low_b1], r_gca_obs[tester_low_b1], col="darkred", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_low_e], r_gca_obs[tester_low_e], col="orange", pch=8, cex=1, lwd=2)
 
legend("bottomleft", legend = c("Testers with Positive b1", 
                                "Testers with Negative b1",
                                "Highest b1", "Zero b1", "Lowest b1","Lowest e"), 
       pt.bg = c("black", "grey70", NA, NA, NA,NA), 
       col = c("white", "white", "darkgreen","royalblue","darkred","orange"), 
       pch = c(21, 21, 8,8,8,8), bty = "n")

# This shows that there are some other testers with low unaligned variance (orange dot) 
# but let's see how that relates to the total variance of the candidate lines GCA.
plot(total_var, r_gca_obs, 
     pch = 21, bg = ifelse(b1 > 0, "black", "grey70"), col = "white", cex = 1,
     ylim = c(0.4,1),
     xlab = "Total Variance", 
     ylab = "Observed r (Hybrid, true GCA)")

points(total_var[tester_high_b1], r_gca_obs[tester_high_b1], col="darkgreen", pch=8, cex=1, lwd=2)
points(total_var[tester_zero_b1], r_gca_obs[tester_zero_b1], col="royalblue", pch=8, cex=1, lwd=2)
points(total_var[tester_low_b1], r_gca_obs[tester_low_b1], col="darkred", pch=8, cex=1, lwd=2)
points(total_var[tester_low_e], r_gca_obs[tester_low_e], col="orange", pch=8, cex=1, lwd=2)

legend("bottomleft", legend = c("Testers with Positive b1", 
                                "Testers with Negative b1",
                                "Highest b1", "Zero b1", "Lowest b1","Lowest e"), 
       pt.bg = c("black", "grey70", NA, NA, NA,NA), 
       col = c("white", "white", "darkgreen","royalblue","darkred","orange"), 
       pch = c(21, 21, 8,8,8,8), bty = "n")

# Now let's calculate the magnitude of aligned variance (useful signal) for
# high_b1 and tester_low_e. As our goal is to maximise the aligned variance, we 
# want the higher value.
c(High_b1 = unname((r_gca_obs[tester_high_b1]^2) * total_var[tester_high_b1])) #154.9575 
c(Low_e   = unname((r_gca_obs[tester_low_e]^2) * total_var[tester_low_e])) #119.951

# Finally, let's see how the tester with the lowest unaligned variance (orange dot)
# aligns GCA-SCA (b1) and how that relates to the observed accuracy.
plot(unaligned_var, b1,
     pch = 19, col = "black", cex=0.5,
     xlim = c(0,50), ylim=c(-0.8,0.6), 
     xlab = "Unaligned Variance", 
     ylab = "Tester Slope (b1)")

points(unaligned_var[tester_high_b1], b1[tester_high_b1], col="darkgreen", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_zero_b1], b1[tester_zero_b1], col="grey40", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_low_b1],  b1[tester_low_b1],  col="darkred", pch=8, cex=1, lwd=2)
points(unaligned_var[tester_low_e],   b1[tester_low_e],   col="orange", pch=8, cex=1, lwd=2)

legend("bottomleft", legend = c("Testers with Positive b1", 
                                "Testers with Negative b1",
                                "Highest b1", "Zero b1", "Lowest b1","Lowest e"), 
       pt.bg = c("black", "grey70", NA, NA, NA,NA), 
       col = c("white", "white", "darkgreen","grey40","darkred","orange"), 
       pch = c(21, 21, 8,8,8,8), bty = "n")



## -----------------------------
## 5) Save testcross parameters
## -----------------------------
if(!dir.exists("./testersel/data")) {
  dir.create("./testersel/data", recursive = TRUE)
}

save(b1,r_gca_obs,var_sca_testers,unaligned_var,
     file = "./testersel/data/S3_testcrosses_data.RData")

# #
# 
# hist(b1 <- cov(sca, gca2)/var(gca2))
# hist(cor(H_centre,gca2))
# 
# hist(cor(sca,gca2)) #selection for pool 1 testers
# cor(sca,gca2)[1]
# cor(sca[,1], gca2)
# # -0.6254651
# 
# cor(H_centre,gca2)[1]
# cor(H_centre[,1], gca2)
# 
# plot(cor(H_centre,gca2), cor(sca,gca2))
# # 0.8134829
