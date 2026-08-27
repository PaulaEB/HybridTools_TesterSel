################################################################################
# Script S.1: Global parameters, simulated data
#
# Authors: Paula Espitia-Buitrago and Daniel Tolhurst 
#
#
# Simulates testcross data parameters from a full factorial cross between
# two complementary heterotic groups:
# 
# - Number of individuals in each heterotic group (pool 1, pool 2): 40
# - Number of biallelic sites with LD along genome: 100
# - Number of QTL with additive (a) + dominance (d): 20 
# - Additive effects: directional and mostly positive
# - Dominance effects: correlated with additive effects and positive dominance 
#   with 50% QTL partial/complete and 50% QTL overdominance
# - True hybrid value centred H_centre (add + dom, no noise)
# - True combining ability components gca_1 (pool 1), gca_2 (pool 2), and sca

# Key outputs:
# - Matrix of complementary heterotic pools (X1, X2)
# - Full factorial hybrid matrix (H_centre)
# - Combining ability components true variance (gca1, gca2, sca)
#
################################################################################

set.seed(1)
source("./functions_TesterSel.R")

## --------------------------
## 1) Simulate LD genotypes
## --------------------------

npool1  <- npool2 <- 40    
nSNP    <- 1000   
nQTL    <- 100
rho_ld  <- 0.8   # LD "strength" along genome (AR1-ish)
p1_base <- 0.3   # baseline allele frequency in pool 1
p2_base <- 0.7   # baseline allele frequency in pool 2
p_sd    <- 0.10  # heterogeneity in allele frequency along sites

# AR(1) latent for site-specific allele frequency variation + mild drift between pools
ld_profile <- function(nSNP, rho) {
  z <- numeric(nSNP); z[1] <- rnorm(1)
  for(k in 2:nSNP) z[k] <- rho*z[k-1] + sqrt(1-rho^2)*rnorm(1)
  z
}
z <- ld_profile(nSNP, rho_ld)   

# Site-specific frequencies for each pool: use the base allele frequencies allowing 
# p_sd heterogeneity along sites and adding some noise. The frequencies are clipped
# to [0.02,0.98]
p1 <- pmin(pmax(p1_base + p_sd*scale(z)[,1] + 2*rnorm(nSNP,0,0.02), 0.02), 0.98)
p2 <- pmin(pmax(p2_base + p_sd*scale(z)[,1] + 2*rnorm(nSNP,0,0.02), 0.02), 0.98)
plot(p1, 1-p2);abline(a=0, b=1); legend("topleft", legend = sprintf("cor = %.3f", cor(p1, 1-p2)), bty = "n") 

# Function to simulate inbred haplotypes with LD via a simple Markov chain of 
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

# Simulate genotype matrix for inbred pools

X1 <- inbred_pool(npool1, p1, rho_switch=0.93)
X2 <- inbred_pool(npool2, p2, rho_switch=0.93)

rownames(X1) <- paste0("p1_",1:npool1)
rownames(X2) <- paste0("p2_",1:npool2)
colnames(X1) <- colnames(X2) <- paste0("S_",1:nSNP)

# Check complementarity between pools
cat("Observed allele frequency in pool 1:", mean(colMeans(X1)/2))
cat("Observed allele frequency in pool 2:", mean(colMeans(X2)/2))

plot(colMeans(X1), apply(X1, 2, sd), col="purple", xlim=c(0,2), ylim=c(0,2))
points( colMeans(X2),apply(X2, 2, sd), col="black")

hist(colMeans(X1)/2, xlim = c(0, 1), col="purple"); hist(colMeans(X2)/2, add = T, col = "grey20")

# LD proxy among sites (correlation of 0/2 genotypes across all inbreds)
Xall <- rbind(X1, X2)
LD <- cor(Xall)
# Heatmap of LD 
heatmap(LD, Rowv = NA, Colv = NA, scale = "none")
#LD profile of a single locus
plot(LD[,999], pch=20, xlab="Sites", ylab="LD");abline(h = 0, col = "red", lty = 2)  

## --------------------------
## 2) Choose QTL + effects
## --------------------------

# Sample randomly the QTL from the sites
qtl_idx <- sort(sample(1:nSNP, nQTL))

# Add the directional (positive) additive and dominance effects
a <- rep(0, nSNP)
d <- rep(0, nSNP)

################################
#####CHECK WITH DT##############
## Directional additive effects:
##   - define allele coded as "2" (reference allele count=2) as favorable at QTL
##   - set additive effects mostly positive
a_q <- abs(rnorm(nQTL, mean=0.8, sd=0.4))  # mostly positive in the QTL
a[qtl_idx] <- a_q

## Directional dominance:
##   - dominance deviations mostly positive and correlated with additive magnitude
##   - mild noise, clipped at >= 0
#d_q <- 0.6*a_q + rnorm(nQTL, mean=0.2, sd=0.15) #mostly partial dominance
d_q <- 1*a_q + rnorm(nQTL, mean=0, sd=0.35) #partial/complete dominance and overdominance
d_q <- pmax(d_q, 0)             # enforce directional (positive) dominance
d_q <- pmin(d_q, 1.8)           # cap
d[qtl_idx] <- d_q

# Check that the the QTL are ~50% partial/complete dominant and ~50% overdominant
table(ifelse(d == 0, "no dominance", ifelse(d > a, "over", "partial/complete")))
plot(p1[qtl_idx], 1-p2[qtl_idx]);abline(a=0, b=1); legend("topleft", legend = sprintf("cor = %.3f", cor(p1[qtl_idx] , 1-p2[qtl_idx])), bty = "n") 
plot(a_q , d_q);abline(a=0, b=1); legend("topleft", legend = sprintf("cor = %.3f", cor(a_q, d_q)), bty = "n") 

################################################################################

## ------------------------------------------------------------
## 3) Full factorial cross and true hybrid value H_centre(i,j)
## ------------------------------------------------------------

# Create the additive matrix for the hybrids. In the rows there are the additive 
# values for each cross p1_i x p2_j, where pool 1 cycles slowly and pool 2 cycles quickly.
# Thus, it starts with parent 1 of pool 1 (p1_1) crossed to each parent of pool 2,
# the it moves to parent 2 of pool 1 (p1_2) crossed to each parent of pool 2
X12 <- (X1[rep(1:npool1, each = npool2), ] + X2[rep(1:npool2, times = npool1), ])/2

# Create a matrix coding 0 for homozygous genotypes and 1 for heterozygous genotypes
W12 <- 1 - abs(X12 - 1)

# Create the hybrid matrix adding additive and dominance effects. For the way X12 
# is organised, we will have pool 1 in the columns and pool 2 in the rows
hybrid_mat <- matrix(X12 %*% a + W12 %*% d, ncol = npool1)
rownames(hybrid_mat) <- rownames(X2)  # Pool 2 in rows
colnames(hybrid_mat) <- rownames(X1)  # Pool 1 in columns

# Centre hybrid matrix
H_centre <- hybrid_mat - mean(hybrid_mat)

## ------------------------------------------------------------
## 4) Variance-component decomposition from the hybrid table
##    Model: Hij = mu + GCA_i + GCA_j + SCA_ij
##    (truth, no residual)
## ------------------------------------------------------------

# Calculate the general combining ability for parents of pool 1 (gca1) and
# parents of pool 2 (gca2)
gca1 <- colMeans(H_centre)
gca2 <- rowMeans(H_centre)

# Calculate the specific combining ability by subtracting the GCAs from the hybrid
# value 
sca <- sweep(sweep(H_centre, 1, gca2, "-"), 2, gca1, "-")

# Calculate the true variance components for combining abilities and the SCA to
# total genetic variance ratio
pop_var <- function(x) { mean((x - mean(x))^2) }
var_gca1 <- pop_var(gca1)
var_gca2 <- pop_var(gca2)
var_sca <- pop_var(as.vector(sca))

cat("Variance ratio SCA to total genetic variance: ", ratio<-var_sca / (var_gca1 + var_gca2+var_sca))

barplot(c(var_gca1, var_gca2, var_sca), 
        names.arg = c("GCA 1", "GCA 2", "SCA"), 
        main = paste("Variance ratio SCA to total genetic:", round(ratio,2)), 
        ylab = "Variance", 
        col = c("purple", "grey20", "orange"))

# Here variance of the GCA 2 is higher than GCA 1, which is expected given that  
# pool 2 has more allele frequency of the favourable allele coded as "2"

## ------------------------------------------------------------
## 5) Save global parameters
## ------------------------------------------------------------
if(!dir.exists("./data")) {
  dir.create("./data", recursive = TRUE)
}

save(npool1, npool2, nSNP, nQTL, 
     p1, p2, X1, X2, a, d, qtl_idx,
     H_centre, gca1, gca2, sca, 
     var_gca1, var_gca2, var_sca, ratio,
     file = "./data/S1_data.RData")
