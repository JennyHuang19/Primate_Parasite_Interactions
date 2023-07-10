# Analyzing the birds and plants data using the proposed method.

# The directory where the analysis is performed:
# wd_path <- '/Users/JennyH/Desktop/PSR_Model/'
wd_path <- '/hpc/group/dunsonlab/yjh3/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where you want to save MCMC results:
save_path <- 'outputs/'
# Where the functions are available:
source_path <- 'HelperScripts/'


# ------ STEP 0: Some functions. --------- #
# ------ STEP 0: Some libraries. --------- #
setwd(wd_path)


source(paste0(source_path, 'UpdOccur_function.R')) # (may 18th)
source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'MCMC_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))
source(paste0(source_path, 'MCMC_trim_new.R'))
# source(paste0(source_path, 'MCMC_function_trimResults.R'))


library(BayesLogit)
library(mvnfast)
library(devtools)
# devtools::install_github("gpapadog/BiasedNetwork", force = TRUE)
library(BiasedNetwork)
# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv.dat'))
load(paste0(data_path, 'obs_A2.dat'))
load(paste0(data_path, 'obs_W2.dat'))
load(paste0(data_path, 'obs_X2.dat'))
load(paste0(data_path, 'obs_F2.dat'))
load(paste0(data_path, 'obs_OB2.dat'))
load(paste0(data_path, 'obs_OP2.dat'))
# load(paste0(data_path, 'birds_232.dat'))
load(paste0(data_path, 'obs_OB2.2.dat'))
load(paste0(data_path, 'obs_OP2.2.dat'))

Cu <- Cu_phylo
Cv <- Cv

# Note that for the analysis in the appendix we used the Cu_tax and Cv_tax
# correlation matrices.

# Restricting to the analysis of the 232 bird species # SKIP FOR THE PRIMATE STUDY..

# Sample sizes of the two sets of species:
nB <- nrow(Cu)
nP <- nrow(Cv)
nStudies <- dim(obs_A2)[3]

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE  # Performing bias correction.


Nsims <- 500
burn <- 10000
thin <- 10 # 20 # 40

use_H <- 10
theta_inf <- 0.01
mh_n_pis <- 70  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 70
mh_n_rho <- 100

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

start_values <- NULL
sampling <- NULL


# --------------- STEP 2: MCMC. ----------------- #

# We run 2 chains. We suggest that you run the following code in parallel instead.
for (cc in 1 : 2) {  # Chain index:
  
  set.seed(cc)
  
  # Running the method:
  mcmc <- MCMC_trimResults(obs_A = obs_A2, focus = obs_F2, occur_B = obs_OB2.2, occur_P = obs_OP2.2,
               obs_X = obs_X2, obs_W = obs_W2, Cu = Cu, Cv = Cv,
               Nsims = Nsims, burn = burn, thin = thin,
               use_H = use_H, bias_cor = bias_cor,
               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
               stick_alpha = stick_alpha, prior_theta = prior_theta,
               prior_tau = prior_tau, prior_rho = prior_rho,
               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
               prior_sigmasq = prior_sigmasq, start_values = start_values,
               sampling = sampling)
  
  
  # Attaching the results:
  attach(mcmc)
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  # and the probability we use when sampling the interaction indicators.
  # Studying MCMC() will clarify the three quantities.
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, pL1s = mcmc$pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = rU, V = rV)
  
  # Combining the results we are interested in to a list and saving:
  res <- list(all_pred = all_pred, correlations = correlations)
  ### trimmed.
  # res <- list(all_pred = mcmc$pL1s_mean)
  ###
  save(res, file = paste0(save_path, 'res.may29th_', cc, '.dat'))
  
  
  rm(res)
  detach(mcmc)
  
}

