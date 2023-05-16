# The directory where the analysis is performed:
# at this point, you should download the results "res"..
wd_path <- '/Users/JennyH/Desktop/PSR_Model/'
# Where the processed data are saved:
data_path <- 'Data/'
# Where the results are saved:
result_path <- 'outputs/'


library(ggplot2)
library(RColorBrewer)
library(gplots)
library(superheat)
library(abind)
library(gridExtra)
library(grid)

setwd(wd_path)

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv.dat'))
load(paste0(data_path, 'obs_A2.dat'))
load(paste0(data_path, 'obs_W2.dat'))
load(paste0(data_path, 'obs_X2.dat'))
load(paste0(data_path, 'obs_F2.dat'))
load(paste0(data_path, 'obs_OB2.dat'))
load(paste0(data_path, 'obs_OP2.dat'))


# The combined network:
comb_A <- apply(obs_A2, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

nB <- nrow(obs_A2)
nP <- ncol(obs_A2)

# Number of MCMC chains for our method and for the alternative method:
nchains <- 2

# Number of cross validation repetitions:
repetitions <- 30

# Covariate names that are nicer for plotting:
good_namesX <- c('Body Mass', 'Gape Size', 'Large*', 'Fruit\nDependent*', 'Endangered*')
good_namesW <- c('Fruit\nDiameter', 'Fruit\nLength', 'Seed\nDiameter', 'Seed\nLength', 'Native*',
                 'Tree*', 'Black\nFruit*', 'Red\nFruit*', 'Yellow/Orange\nFruit*', 'Green\nFruit*',
                 'Lipid*', 'Endangered*')


# --------------- STEP 1: Getting the results together ----------------- #

all_res <- NULL
for (ii in 1 : nchains) {
  load(paste0(result_path, 'restest6_', ii, '.dat'))
  all_res[[ii]] <- res # loads the "res" object, saved in restest6_1.dat
}



# --------------- STEP 2: Plotting our analysis results ----------------- #

# Binding the posterior samples for the interactions from across chains.


# ----- Based on our model:

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1] # N_sim, 97x157, 3

# (apr 9) hunch: the all pred object changed dims after we used the MCMC_function_trimResults.

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nB, nP))
for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1]
}
dimnames(pred_ours)[2 : 3] <- list(bird = rownames(obs_A), plant = colnames(obs_A))

# Calculating the posterior probability of interaction by averaging across
# posterior samples:
pred_ours <- apply(pred_ours, c(2, 3), mean)

mean(pred_ours) # 0.898

# confusionmatrix (first), ROC


# ----------- PART A: PLOTTING THE HEATMAP:

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
bird_group <- bird_order_info$Frug_Family
plant_group <- as.character(plant_order_info$Plant_family)

# primate genus and parasite type. 

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
bird_size_cluster <- sapply(unique(bird_group), function(x) sum(bird_group == x))
plant_size_cluster <- sapply(unique(plant_group), function(x) sum(plant_group == x))

# Set plot_pred to pred_ours for results based on our method and to
# pred_alt for results based on the alternative method:
plot_pred <- pred_ours

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_bird_size <- 10
min_plant_size <- 20

keep_bird_groups <- names(which(bird_size_cluster >= min_bird_size))
keep_plant_groups <- names(which(plant_size_cluster >= min_plant_size))

keep_bird_index <- which(bird_group %in% keep_bird_groups)
keep_plant_index <- which(plant_group %in% keep_plant_groups)

# Plotting those with minimum size as specified:
superheat(plot_pred,
          # membership.rows = bird_group[keep_bird_index],
          # membership.cols = plant_group[keep_plant_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          force.bottom.label = TRUE,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.24, left.label.size = 0.24,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05))

# plot: set anything in A black.

# Array of posterior probabilities.
# df of host, parasite, posterior probabilities.
adj_list <- lapply(rownames(plot_pred),function(x)sapply(colnames(plot_pred),function(y)list(x,y,plot_pred[x,y])))
adj_listtmp <- matrix(unlist(adj_list),nrow=3)
bipartitelist <- t(adj_listtmp)
colnames(bipartitelist) <- c("Primate","Parasite","posterior_prob")
bipartitedf <- as.data.frame(bipartitelist,stringsAsFactors=F)
bipartitedf[,3] <- as.numeric(bipartitedf[,3])


# arrange in desc order of post. prob.
bipartitedf <- bipartitedf %>% 
  filter(posterior_prob < 1) %>% 
  filter(posterior_prob >= 0.9) %>% 
  arrange(desc(posterior_prob))

head(bipartitedf)
View(bipartitedf)
# https://rlbarter.github.io/superheat/clustering.html

# ----- MCMC Diagnostics:

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

## Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1] # N_sim, 97x157, 3

## rho val
rhoV = all_res[[1]]$correlations[,2]
rhoU = all_res[[1]]$correlations[,1]
probL = all_res[[1]]$all_pred[,,,3] 
probL11 = all_res[[1]]$all_pred[,1,1,3]
probL12 = all_res[[1]]$all_pred[,1,2,3]

## convert the output into coda format
coda_mcmc_chain = mcmc(cbind(rhoV, rhoU, probL11, probL12))
summary(coda_mcmc_chain)

## look at the menu options
codamenu()


