################################################################################
# Functions for Tester Selection 
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
#
# This script provides the core functions to demonstrate key concepts in tester 
# selection for hybrid breeding programs:
# 
# 1. Simulation of LD profile using AR(1) spatial model (S1_GlobalParameters)
# 2. Simulation of inbred pools considering recombination rates (S1_GlobalParameters)
# 3. Random single tester or set of testers selection (S2_homogeneousCAvariances)
# 
# Key outputs:
# - Latent vector representing LD correlation across genomic sites
# - Fully homozygous genotype matrices (0/2 coded) for heterotic pools
# - Lists and tracking dataframes specifying random combinations of testers 
#
################################################################################

## --------------------------------------------------
## 0) Function to calculate the population variance
## --------------------------------------------------

pop_var <- function(x) { mean((x - mean(x))^2) }
pop_cov <- function(x, y) { 
  if (is.matrix(x)) {
    apply(x, 2, function(col) mean((col - mean(col)) * (y - mean(y))))
  } else {
    mean((x - mean(x)) * (y - mean(y)))
  }
}
## -------------------------------------
## 1) Function to simulate LD profile
## -------------------------------------

# This function uses the first-order auto-regressive AR(1) spatial model to generate
# a latent vector. Once scaled, this vector is used to introduce heterogeneity in 
# baseline allele frequencies along the sites, simulating LD where adjacent loci
# share similar frequencies.

ld_profile <- function(nSNP, rho) {
  z <- numeric(nSNP); z[1] <- rnorm(1)
  for(k in 2:nSNP) z[k] <- rho*z[k-1] + sqrt(1-rho^2)*rnorm(1)
  z
}

## -------------------------------------
## 2) Function to simulate inbred pools
## -------------------------------------

# This function simulates inbred haplotypes with LD via a simple Markov chain of 
# alleles along sites. Each inbred genotype is 0/2 coded (fully homozygous at each site)
# - rho_switch simulates the recombination rate, i.e., probability to "stay" with 
#   previous allele state, with 92% by default ~ 8% recombination rate

inbred_pool <- function(n, p_vec, rho_switch=0.92) { 
  X <- matrix(0, n, length(p_vec))
  for(i in 1:n){
    a <- rbinom(1, 1, p_vec[1])
    X[i,1] <- 2*a
    for(k in 2:length(p_vec)){
      if(runif(1) < rho_switch){
        a <- X[i,k-1]/2
      } else {
        a <- rbinom(1, 1, p_vec[k])
      }
      X[i,k] <- 2*a
    }
  }
  X
}

## ---------------------------------------
## 3) Function to select testers randomly
## ---------------------------------------

# This function dictates the sampling strategy for selecting different numbers of 
# testers (t). For t=1 and t=2, it computes all combinations using combn(). 
# For larger t values yielding combinations beyond 'max_iter', it shifts to sample
# 5000 combinations to preserve computational efficiency.
#
# Arguments:
# - t: numeric vector listing the sizes of tester sets to generate (e.g., c(1,2,3...))
# - max_iter: threshold to restrict combinations, set to 5000 by default
# - tester_pool: character string defining the target pool ("pool1" or "pool2")
# - df: if TRUE, returns a list containing both the raw sets and a visual dataframe

randomsel_testers <- function(t, max_iter = 5000, tester_pool = "pool1", df = FALSE) {
  
  # Define the tester pool
  if (tester_pool == "pool1") {
    pool_matrix <- X1
  } else if (tester_pool == "pool2") {
    pool_matrix <- X2
  } else {
    stop("Check tester_pool, it must be pool1 or pool2")
  }
  
  tester_idx <- rownames(pool_matrix)
  n_total    <- length(tester_idx)
  
  # Iterate over t values
  all_sets_list <- do.call(c, lapply(t, function(current_t) {
    sets_list <- list() 
    total_combos <- choose(n_total, current_t)
    
    if (total_combos <= max_iter) {
      cat("  -> Generating:", total_combos, "combinations for t =", current_t, "\n")
      all_combos <- combn(tester_idx, current_t) 
      for (i in 1:ncol(all_combos)) {
        sets_list[[i]] <- all_combos[, i]
      }
    } else {
      cat("  -> Sampling:", max_iter, "combinations from", total_combos, "for t =", current_t, "\n")
      for (i in 1:max_iter) {
        sets_list[[i]] <- sample(tester_idx, current_t, replace = FALSE)
      }
    }
    return(sets_list)
  }))
  
  # Return a dataframe with the testers if df = T
  if (df == TRUE) {
    # Visual dataframe 
    visual_df <- data.frame(
      n_testers  = sapply(all_sets_list, length),
      tester_set = sapply(all_sets_list, paste, collapse = ","),
      stringsAsFactors = FALSE
    )
    
    # And the list
    return(list(
      sets = all_sets_list,
      df   = visual_df
    ))
    
  } else {
    # If df = FALSE return only the list to do everything below
    return(all_sets_list)
  }
}

## ---------------------------------------
## 4) Function to evaluate testcross sets
## ---------------------------------------

# This function evaluates the genetic parameters of candidate lines when crossed
# with different tester sets. It simulates a testcross design for each tester set
# generated by the 'randomsel_testers' function, and extracts both the true 
# parameters of the full factorial and the testcross components to allow for accuracy
# assessments.
#
# Key calculations per testcross:
# - GCA estimates for testers and candidates 
# - SCA estimates (if t > 1).
# - Trait performances (Mean and Maximum) to assess hybrids potential.
# - Expected vs Observed Accuracy (Correlations) and Variances for GCA and SCA.
#
# Note: For sets with a single tester (t = 1), SCA variance and alignment are 
# technically unestimable (assigned as NA), as GCA and SCA are perfectly confounded.

evaluate_tc_sets <- function(tester_sets_list, H_mat, gca2_true, sca_true, 
                             var_gca1, var_gca2, var_sca) {
  if (is.list(tester_sets_list) && "sets" %in% names(tester_sets_list)) {
    tester_sets_list <- tester_sets_list$sets
  }
  cat("\nEvaluating", length(tester_sets_list), "testcrosses...\n")
  
  comb_ability_tc <- lapply(tester_sets_list, function(testers_idx) {
    
    t <- length(testers_idx)
    tester_set <- paste(sort(testers_idx), collapse = ",")
    
    H_centre_tc <- H_mat[, testers_idx, drop = FALSE]
    
    gca2_tc <- data.frame(p2_id = rownames(H_centre_tc),
                          GCA_tc = rowMeans(H_centre_tc))
    mean_gca1_tc <- mean(gca1[testers_idx])
    hybrid_perf_tc <- mean(H_centre_tc)
    maxhybrid_perf_tc <- max(H_centre_tc)
    
    var_sca_tc_val <- NA
    sca_align_val  <- NA
    var_gca1_tc_val <- NA
    sca_acc_expected <- NA
    if (t > 1) {
      gca1_tc <- data.frame(p1_id = colnames(H_centre_tc),
                            GCA_tc = colMeans(H_centre_tc))
      
      ######### --PAULA-- Check this logic #############
      sca_obs <- sweep(sweep(H_centre_tc, 1, gca2_tc$GCA_tc, "-"), 2, gca1_tc$GCA_tc, "-")
      true_sca_tc <- sca_true[, testers_idx, drop = FALSE]
      
      var_sca_tc_val <- popVar(as.vector(sca_obs))
      sca_align_val  <- cor(as.vector(true_sca_tc), as.vector(sca_obs))
      var_gca1_tc_val <- popVar(gca1_tc$GCA_tc)
      sca_acc_expected <- sqrt(abs(1 - 1 / t))
    }
    
    return(data.frame(
      n_testers = t,
      tester_set = tester_set,
      tester_mean_gca1 = mean_gca1_tc,
      hybrid_mean_perf = hybrid_perf_tc,
      maxhybrid_perf_tc = maxhybrid_perf_tc,
      true_var_gca2 = var_gca2,
      var_gca2_expected = var_gca2 + var_sca / t,
      true_var_sca  = var_sca,
      var_sca_expected  = var_sca * (1 - 1 / t),
      true_var_gca1 = var_gca1,
      var_gca1_tc   = var_gca1_tc_val,
      var_gca2_tc   = popVar(gca2_tc$GCA_tc),
      var_sca_tc    = var_sca_tc_val,
      gca2_acc_expected = sqrt(var_gca2 / (var_gca2 + var_sca / t)),
      gca_align     = cor(gca2_true, gca2_tc$GCA_tc),
      sca_acc_expected  = sca_acc_expected,
      sca_align     = sca_align_val,
      stringsAsFactors = FALSE
    ))
  })
  return(do.call(rbind, comb_ability_tc))
}

## -------------------------------------------------------------
## 5) Function to calculate the representativeness score for a 
##    single tester 
## -------------------------------------------------------------

# This function calculates the representativeness score for a single tester based 
# on its distance to the average allele frequencies of the target pool. The score
# can be computed using either a weighted or unweighted Euclidean distance, depending 
# on the user's weight to certain loci. 

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

# -------------------------------------------------------------
# 6) Function to calculate the representativeness score for a
#   set of testers
# Representativeness: weighted squared distance between allele freqs in chosen 
# testers and the true base allele frequencies of the pool (p1). Lower distance 
# is better, thus closer to 0 is perfectly representative.
# rep_scores_1 <- sapply(all_testers, function(id) {
#   tester_freq <- X1[id, ] / 2
#   -sum(w_qtl_norm * (tester_freq - target_pool_freqs)^2)
# })
# 
# represent_score <- function(weights, tester_idx, tester_pool = "pool1") {
# if (tester_pool == "pool1") {
#   pool_matrix <- X1
#   target_p <- p1
# } else if (tester_pool == "pool2") {
#   pool_matrix <- X2
#   target_p <- p2
# } else {
#   stop("Check tester_pool, it must be 'pool1' or 'pool2'")
# }
# 
# pT <- colMeans(pool_matrix[tester_idx, , drop = FALSE] / 2)
# -sum(weights * (pT - target_p)^2)
# }

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
