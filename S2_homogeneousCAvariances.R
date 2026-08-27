################################################################################
# Script S.2: Homogeneous variance in combining ability components
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
#
# Demonstration for increasing the number of testers in a testcross increase the
# accuracy of the general combining ability (GCA) for candidate lines, the GCA 
# variance and the hybrid value (H_centre), while approaching the true proportion  
# of specific combining ability (SCA) in the total genetic variance:
# 
# - Assumes that GCA and SCA are independent with homogeneous variance
# - Number of testers to evaluate: 1, 2, 3, 4, 5
# 
# Key outputs:
# - Comparison of true and expected variance components for GCA and SCA
# - Variation among different sets of testers for GCA and SCA variance components
#   in the testcrosses
# - Alignment of GCA and SCA in the testcrosses with the true GCA and SCA
# - Comparison of the mean and maximum hybrid value in the testcrosses with the 
#   true, and the effect of tester's GCA
#
################################################################################
rm(list = ls())

source("./functions_TesterSel.R")
load("./data/S1_data.RData")

## ---------------------------------------------------------
## 1) Simulate testcrosses with different number of testers
##    from pool 1 randomly selected   
## ---------------------------------------------------------

###################### Function to choose testers ##########################

randomsel_testers <- function(t, max_iter = 5000, tester_pool = "pool1", df = FALSE) {
  
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

############################################################################

# Simulate testcrosses with 1, 2, 3, 4 and 5 testers randomly selected. For 1 and 
# 2 testers it samples all the possible sets. For more than 2 testers it samples 
# maximum 5000 sets

ntesters <- c(1,2,3,4,5)
tester_sets <- randomsel_testers(t = ntesters, max_iter = 5000, 
                                 tester_pool = "pool1", df = T)

## -------------------------------------------
## 2) Extract combining ability components 
## -------------------------------------------

#################### Function to simulate testcrosses ########################
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
    
    pop_var <- function(x) { mean((x - mean(x))^2) }
    var_sca_tc_val <- NA
    sca_align_val  <- NA
    var_gca1_tc_val <- NA
    sca_acc_expected <- NA
    if (t > 1) {
      gca1_tc <- data.frame(p1_id = colnames(H_centre_tc),
                            GCA_tc = colMeans(H_centre_tc))
      
      sca_obs <- sweep(sweep(H_centre_tc, 1, gca2_tc$GCA_tc, "-"), 2, gca1_tc$GCA_tc, "-")
      true_sca_tc <- sca_true[, testers_idx, drop = FALSE]
      
      var_sca_tc_val <- pop_var(as.vector(sca_obs))
      sca_align_val  <- cor(as.vector(true_sca_tc), as.vector(sca_obs))
      var_gca1_tc_val <- pop_var(gca1_tc$GCA_tc)
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
      var_sca_expected  = var_sca * (1 - 1 / t), #this assumes homogeneous variance
      true_var_gca1 = var_gca1,
      var_gca1_tc   = var_gca1_tc_val,
      var_gca2_tc   = pop_var(gca2_tc$GCA_tc),
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

############################################################################

# Use the tester sets to obtain combining abilities 

comb_ability_df <- evaluate_tc_sets(tester_sets_list = tester_sets$sets, 
                                    H_mat = H_centre, 
                                    gca2_true = gca2, 
                                    sca_true = sca, 
                                    var_gca1 = var_gca1, 
                                    var_gca2 = var_gca2, 
                                    var_sca = var_sca)

#### Plots ##########

# GCA of candidate lines in pool2 and SCA
plot(comb_ability_df$n_testers, comb_ability_df$true_var_gca2, col="blue");points(comb_ability_df$n_testers,comb_ability_df$var_gca2_expected, col="black")
plot(comb_ability_df$n_testers, comb_ability_df$true_var_sca, col="red", ylim=c(0,30));points(comb_ability_df$n_testers,comb_ability_df$var_sca_expected, col="black")

boxplot(var_gca2_tc ~ n_testers, data=comb_ability_df);abline(h=comb_ability_df$true_var_gca2, col = "blue")
boxplot(var_sca_tc ~ n_testers, data=comb_ability_df);abline(h=comb_ability_df$true_var_sca, col = "red")

boxplot(gca_align ~ n_testers, data=comb_ability_df);points(comb_ability_df$n_testers, comb_ability_df$gca2_acc_expected, col="green")
boxplot(sca_align ~ n_testers, data=comb_ability_df);points(comb_ability_df$n_testers-1, comb_ability_df$sca_acc_expected, col="green")

# Checking SCA to total variance ratio in tester sets
gca1_subset<-subset(comb_ability_df, n_testers > 1)
gca1_subset$ratio_sca_total <- gca1_subset$var_sca_tc / (gca1_subset$var_gca1_tc + gca1_subset$var_gca2_tc + gca1_subset$var_sca_tc)
boxplot(ratio_sca_total ~ n_testers, data=gca1_subset); abline(h = ratio, col = "red")

# Checking testers GCA and hybrid performance
tc1<-subset(comb_ability_df, n_testers == 1)
barplot(tc1$tester_mean_gca1, names.arg = tc1$tester_set,las=2, ylab = "Tester GCA")
boxplot(hybrid_mean_perf ~ n_testers, data = comb_ability_df);abline(h = max(comb_ability_df$hybrid_mean_perf[comb_ability_df$n_testers==1]), col = "darkgreen", lwd = 2, lty = 3)
boxplot(maxhybrid_perf_tc ~ n_testers, data = comb_ability_df);abline(h = max(comb_ability_df$maxhybrid_perf_tc[comb_ability_df$n_testers==1]), col = "darkgreen", lwd = 2, lty = 3)

# ## -----------------------------
# ## 3) Save testcross parameters
# ## -----------------------------
# if(!dir.exists("./testersel/data")) {
#   dir.create("./testersel/data", recursive = TRUE)
# }
# 
# save(comb_ability_df,
#      file = "./testersel/data/S2_testcrosses_data.RData")

