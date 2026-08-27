################################################################################
# Script S.4: Representativeness
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
# This script evaluates the concept of representativeness to identify a tester 
# that represents its own pool. Thus, it reflects how close a single tester is to
# the mean of its own heterotic pool across all the sites or the QTL, aiming to 
# minimise the mean squared distance (Euclidean distance) between the tester and 
# its reference pool.
# 
# So it contrasts structural genetic distance (unweighted) against functional 
# genetic distance (weighted by QTL). 
#
# Assessed Tester Profiles (identified in S.3):
# 1. The Classic Tester: High b1 (Amplifies GCA).
# 2. The Parent Tester: Low/Negative b1 (Masks GCA via SCA noise).
# 3. The Representative Tester: Lowest weighted Euclidean distance to the pool mean.
# 4. The Low Noise Tester: Lowest unaligned variance (Low e).
# 
# Key outputs:
# - Calculation of weighted (functional) and unweighted (structural) Euclidean 
#   distances to define tester representativeness.
# - Global PCA of both heterotic pools projecting the centroid shifts to demonstrate
#   how different testers drive the candidate population's heterotic divergence
################################################################################

rm(list = ls())

source("./functions_TesterSel.R")
load("./data/S1_data.RData")
# load("./testersel/data/S3_testcrosses_data.RData")

## --------------------------------------------------------------------
## 1) Within pool metrics: unweighted vs weighted representativeness
## ---------------------------------------------------------------------
all_testers <- rownames(X1)

# Calculate the mean allele frequencies of the tester pool at all the sites
x_bar_l <- colMeans(X1) 

# QTL weights: total "importance" at QTL
w_qtl <- (a^2 + d^2) 

w_qtl_norm <- w_qtl / sum(w_qtl) 

# Function to calculate the representativeness score for all the testers 
rep_score <- function(tester_id, weights = NULL) {
  x_jl <- X1[tester_id, ]
  if (!is.null(weights)) {
    # Weighted Euclidean distance
    sqrt(sum(weights * (x_jl - x_bar_l)^2))
  } else {
    # Standard Euclidean distance (all loci treated equally)
    sqrt(sum((x_jl - x_bar_l)^2))
  }
}

# Calculate Representativeness for all single testers 
rep_scores <- sapply(all_testers, rep_score)
rep_scores_w <- sapply(all_testers, function(id) rep_score(id, weights = w_qtl_norm))

# Identification of best/worst profiles based on minimum distance
best_rep_un     <- names(which.min(rep_scores))  # p1_29
low_rep_un    <- names(which.max(rep_scores))  # p1_10

best_rep_w  <- names(which.min(rep_scores_w)) # p1_18
low_rep_w <- names(which.max(rep_scores_w)) # p1_10

# Identify the testers with the b1 and unaligned variances from Script S.3 

cov_gca2_sca <- apply(sca, 2, function(sca_col) pop_cov(sca_col, gca2))
b1 <- cov_gca2_sca / var_gca2

var_sca_testers <- apply(sca, 2, pop_var)
unaligned_var <- var_sca_testers - (b1^2 * var_gca2)

high_b1 <- names(which.max(b1)) # p1_6
low_b1 <- names(which.min(b1)) # p1_10
zero_b1 <- names(which.min(abs(b1))) # p1_3
low_e <- names(which.min(abs(drop(unaligned_var)))) # p1_18

# For weighted distance to work properly, we scale the matrix by the square root
# of the normalized weights.
X1_weighted <- sweep(X1, 2, sqrt(w_qtl_norm), "*")

# Compute tester pool PCA
pca_un <- prcomp(X1, center = TRUE, scale. = FALSE)
pca_w  <- prcomp(X1_weighted, center = TRUE, scale. = FALSE)

# Variance explained by each principal component 
round(summary(pca_un)$importance[2, ] * 100, 2)
round(summary(pca_w)$importance[2, ] * 100, 2)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
barplot(summary(pca_un)$importance[2, ] * 100,
        main = "PCA (No weights)",las=2)
barplot(summary(pca_w)$importance[2, ] * 100,
        main = "PCA (Weights)",las=2)
par(mfrow = c(1, 1)) 

# Plot tester pool PCA with the best/worst profiles highlighted
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
# 1. PCA unweighted
plot(pca_un$x[,1], pca_un$x[,2], pch = 16, col = "grey80",
     main = "PCA (No weights)")
points(0, 0, pch = 3, cex = 1.5, lwd = 3, col = "black")
points(pca_un$x[best_rep_un, 1], pca_un$x[best_rep_un, 2], col = "purple", pch = 19, cex = 1.5)
points(pca_un$x[low_rep_un, 1], pca_un$x[low_rep_un, 2], col = "red", pch = 19, cex = 1.5)
points(pca_un$x[low_e, 1], pca_un$x[low_e, 2], col = "orange", pch = 4, cex = 1.5, lwd = 2)
points(pca_un$x[high_b1, 1], pca_un$x[high_b1, 2], col = "darkgreen", pch = 4, cex = 1.5, lwd = 2)
points(pca_un$x[zero_b1, 1], pca_un$x[zero_b1, 2], col = "royalblue", pch = 4, cex = 1.5, lwd = 2)
points(pca_un$x[low_b1, 1], pca_un$x[low_b1, 2], col = "darkred", pch = 4, cex = 1.5, lwd = 2)

legend("topleft", legend = c("TPG Centroid", "Best Rep", "Low Rep","Low e", "High b1", "Zero b1", "Low b1"),
       col = c("black", "purple", "red","orange", "darkgreen", "royalblue", "darkred"), 
       pch = c(3, 19, 19, 4, 4, 4, 4), pt.cex = 1.2, bty = "n", cex = 0.7)

# 2. PCA weighted
plot(pca_w$x[,1], pca_w$x[,2], pch = 16, col = "grey80",
     main = "PCA (Weights)")
points(0, 0, pch = 3, cex = 1.5, lwd = 3, col = "black")
points(pca_w$x[best_rep_w, 1], pca_w$x[best_rep_w, 2], col = "purple", pch = 19, cex = 1.5)
points(pca_w$x[low_rep_w, 1], pca_w$x[low_rep_w, 2], col = "red", pch = 19, cex = 1.5)
points(pca_w$x[low_e, 1], pca_w$x[low_e, 2], col = "orange", pch = 4, cex = 1.5, lwd = 2)
points(pca_w$x[high_b1, 1], pca_w$x[high_b1, 2], col = "darkgreen", pch = 4, cex = 1.5, lwd = 2)
points(pca_w$x[zero_b1, 1], pca_w$x[zero_b1, 2], col = "royalblue", pch = 4, cex = 1.5, lwd = 2)
points(pca_w$x[low_b1, 1], pca_w$x[low_b1, 2], col = "darkred", pch = 4, cex = 1.5, lwd = 2)
par(mfrow = c(1, 1)) 

## -----------------------------------------------------------------------------
## 2) Simulating Evolutionary Direction: Joint PCA and Selection Shift
## -----------------------------------------------------------------------------

# Impact on the candidate pool: We simulate testcross selection to demonstrate 
# how different tester profiles drive the allele frequencies of the next generation
# of candidates (heterotic divergence).

sel_intensity <- 0.50 
n_select <- max(1, round(npool2 * sel_intensity))

# Global PCA unweighted 
Xall <- rbind(X1, X2)
pca_all <- prcomp(Xall, center = TRUE, scale. = FALSE)
pca_X1 <- pca_all$x[1:npool1, 1:2]
pca_X2 <- pca_all$x[(npool1 + 1):(npool1 + npool2), 1:2]
centroid_X1 <- colMeans(pca_X1)
centroid_X2 <- colMeans(pca_X2)

# Global PCA with weights 
Xall_w <- sweep(Xall, 2, sqrt(w_qtl_norm), "*")
pca_all_w <- prcomp(Xall_w, center = TRUE, scale. = FALSE)
pca_X1_w <- pca_all_w$x[1:npool1, 1:2]
pca_X2_w <- pca_all_w$x[(npool1 + 1):(npool1 + npool2), 1:2]
centroid_X1_w <- colMeans(pca_X1_w)
centroid_X2_w <- colMeans(pca_X2_w)

#Function to plot the shift in the centroid of the selected candidates in PCA space
plot_centroid_shift <- function(tester_id, color_arrow, label = NULL, weighted = FALSE) {
  # Extract tester column in H_centre
  tc_performance <- H_centre[, tester_id] 
  
  # Identify the candidates for selection in X2 using this tester
  selected_idx <- order(tc_performance, decreasing = TRUE)[1:n_select]
  
  # Select the PCA and centroid based on whether weighted or unweighted
  if (weighted) {
    target_pca <- pca_X2_w
    origin_centroid <- centroid_X2_w
  } else {
    target_pca <- pca_X2
    origin_centroid <- centroid_X2
  }
  
  # Get the centroid of the selected candidates in PCA space
  sel_centroid <- colMeans(target_pca[selected_idx, , drop = FALSE])
  
  # Plot the arrow from the original centroid to the selected centroid
  arrows(x0 = origin_centroid[1], y0 = origin_centroid[2], 
         x1 = sel_centroid[1], y1 = sel_centroid[2], 
         col = color_arrow, lwd = 1.5, length = 0.1)
  
  # Include the text if required
  if (!is.null(label)) {
    text(sel_centroid[1], sel_centroid[2], labels = label, 
         col = color_arrow, pos = 3, cex = 0.8, font = 2)
  }
}

# -------------------------------------------------------
# Global PCA visualizations (Shift of candidate lines)
# -------------------------------------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Unweighted
plot(pca_all$x[,1], pca_all$x[,2], pch = 16, col = "grey80",
     main = "PCA all no weights")
points(pca_X1[,1], pca_X1[,2], col = "purple", pch = 1,cex=0.5)
points(pca_X2[,1], pca_X2[,2], col = "black", pch = 1,cex=0.5)

points(centroid_X1[1], centroid_X1[2], pch = 3, cex = 1.5, lwd = 3, col = "purple")
points(centroid_X2[1], centroid_X2[2], pch = 3, cex = 1.5, lwd = 3, col = "black")

points(pca_all$x[best_rep_un, 1], pca_all$x[best_rep_un, 2], col = "purple", pch = 19, cex = 1.5)
points(pca_all$x[low_rep_un, 1], pca_all$x[low_rep_un, 2], col = "red", pch = 19, cex = 1.5)
points(pca_all$x[low_e, 1], pca_all$x[low_e, 2], col = "orange", pch = 4, cex = 1.5, lwd = 2)
points(pca_all$x[high_b1, 1], pca_all$x[high_b1, 2], col = "darkgreen", pch = 4, cex = 1.5, lwd = 2)
points(pca_all$x[zero_b1, 1], pca_all$x[zero_b1, 2], col = "royalblue", pch = 4, cex = 1.5, lwd = 2)
points(pca_all$x[low_b1, 1], pca_all$x[low_b1, 2], col = "darkred", pch = 4, cex = 1.5, lwd = 2)

plot_centroid_shift(best_rep_un, "purple")
plot_centroid_shift(low_rep_un, "red")
plot_centroid_shift(low_e, "orange")
plot_centroid_shift(high_b1, "darkgreen")
plot_centroid_shift(zero_b1, "royalblue")
plot_centroid_shift(low_b1, "darkred")

# Weighted
plot(pca_all_w$x[,1], pca_all_w$x[,2], pch = 16, col = "grey80",
     main = "PCA weights")
points(pca_X1_w[,1], pca_X1_w[,2], col = "purple", pch = 1,cex=0.5)
points(pca_X2_w[,1], pca_X2_w[,2], col = "black", pch = 1,cex=0.5)

points(centroid_X1_w[1], centroid_X1_w[2], pch = 3, cex = 1.5, lwd = 3, col = "purple")
points(centroid_X2_w[1], centroid_X2_w[2], pch = 3, cex = 1.5, lwd = 3, col = "black")

points(pca_all_w$x[best_rep_w, 1], pca_all_w$x[best_rep_w, 2], col = "purple", pch = 19, cex = 1.5)
points(pca_all_w$x[low_rep_w, 1], pca_all_w$x[low_rep_w, 2], col = "red", pch = 19, cex = 1.5)
points(pca_all_w$x[low_e, 1], pca_all_w$x[low_e, 2], col = "orange", pch = 4, cex = 1.5, lwd = 2)
points(pca_all_w$x[high_b1, 1], pca_all_w$x[high_b1, 2], col = "darkgreen", pch = 4, cex = 1.5, lwd = 2)
points(pca_all_w$x[zero_b1, 1], pca_all_w$x[zero_b1, 2], col = "royalblue", pch = 4, cex = 1.5, lwd = 2)
points(pca_all_w$x[low_b1, 1], pca_all_w$x[low_b1, 2], col = "darkred", pch = 4, cex = 1.5, lwd = 2)

plot_centroid_shift(best_rep_w, "purple", weighted = TRUE)
plot_centroid_shift(low_rep_w, "red", weighted = TRUE)
plot_centroid_shift(low_e, "orange", weighted = TRUE)
plot_centroid_shift(high_b1, "darkgreen", weighted = TRUE)
plot_centroid_shift(zero_b1, "royalblue", weighted = TRUE)
plot_centroid_shift(low_b1, "darkred", weighted = TRUE)
par(mfrow = c(1, 1))
par(mfrow = c(1, 1)) 
