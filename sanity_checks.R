
# ------ Sanity Check: observed vs. predicted prevalence  --------- #

# and then compute the interaction prevalence in the data 
# and compare it to the mean of the posterior mean pi_hats. 
# The goal is to be able to say 
# “in the data, we see about 5% of the interactions that were observable were actually observed 
# and the model predicts an interaction prevalence of 7% - not too far off”; 
# we just want to verify that the model isn’t making interactions much more or less likely than in the data.

# METHOD A: mean( (primate, observed parasite, study_i) / (primate, potential parasites, study_i) )

obs_F_ct <- array(0, dim = c(nM, nP, nS))
dimnames(obs_F_ct) <- list(uni_primates, uni_parasites, uni_studies_lst)

# for each study: compute the number of observed parasite / number of potential parasites -- append it to a list.
prevalences <- array(0, dim = nS)

for (ss in 1 : nS) {
  wh_study <- uni_studies_lst[ss]
  
  # i. Find all primates studied in the study.
  these_primates <- unique(traits_and_studies$HostCorrectedName[traits_and_studies$Citation == wh_study]) # the primates in the study
  
  # prevalence: find all parasites observed in the study.
  observed_parasites <- unique(traits_and_studies$ParasiteCorrectedName[traits_and_studies$Citation == wh_study]) # the parasites observed in the study
  
  # ii. find the indices of all 'these_primates'
  wh_primates <- which(uni_primates %in% these_primates) # the indices of the primate.
  
  # iii. Find all parasite types studied in the study.
  wh_study_parasitetypesstudied <- unique(traits_and_studies$parasite_type[traits_and_studies$Citation == wh_study])
  
  # iv. For each parasite type, find all species in the database. We assume the test performed in the study has the potential to record such species.
  these_parasites <- unique(traits_and_studies$ParasiteCorrectedName[traits_and_studies$parasite_type %in% wh_study_parasitetypesstudied])
  
  # v. find the indices of all 'these_parasites'
  wh_parasites <- which(uni_parasites %in% these_parasites)
  
  obs_F_ct[wh_primates, wh_parasites, ss] <- 1 # records all of the potentially recorded parasites for the primates in the study.
  
  prevalences[ss] <- length(observed_parasites) / length(these_parasites)
  
}

# prevalences[i] = (primate, observed parasite, study_i) / (primate, potential parasites, study_i)
sum(prevalences) / length(prevalences) # average prevalence (0.0483)

# METHOD B: 
# prevalence: mean(A | F = 1) = out of the pairs where focus = 1, how many pairs are observed?
mean(obs_A2[obs_F2 == 1])  # 0.0312

# Prevalence predicted by the model:

# posterior probability of interaction by averaging across posterior samples:
# pred_ours from 5a.
mean(pred_ours) # 0.898
mean(pred_ours[comb_F == 1]) # 0.887. posterior probability of interaction that were the focus of some study.

# we need pairs that were ever studied in the same study.
# The combined network:
comb_F <- apply(obs_F2, c(1, 2), sum)
comb_F <- (comb_F > 1) * 1 # pairs that are the focus in >5 studies.
comb_F

# ------ Sanity Check: histogram of number of primates (top entries vs. entire dataset)  --------- #

# for each study: compute the number of observed parasite / number of potential parasites -- append it to a list.
num_parasites_tophosts <- array(0, dim = 50)

for (p in 1:50){
  primate <- bipartitedf$Primate[p]
  num_parasites <- traits_and_studies %>% 
    filter(HostCorrectedName == primate) %>% 
    select(ParasiteCorrectedName) %>% 
    distinct() %>% 
    count()
  num_parasites_tophosts[p] <- num_parasites$n
  
}

hist(num_parasites_tophosts, main="Number of Parasites Top Predicted Host Edges")

num_parasites_lst <- array(0, dim = 1000)
# number of parasites per observed host
for (i in 1:1000){
  primate <- traits_and_studies$HostCorrectedName[i]
  num_parasites <- traits_and_studies %>% 
    filter(HostCorrectedName == primate) %>% 
    select(ParasiteCorrectedName) %>% 
    distinct() %>% 
    count()
  num_parasites_lst[i] <- num_parasites$n
}
hist(num_parasites_lst, main="Number of Parasites All Hosts")

# ------ Sanity Check: model predictions for observed interactions  --------- #
sum(pred_ours == 1)
sum(comb_A) # comb_A: observed interactions combined across all studies. 

pred_ours[comb_A == 1] # posterior_mean for observed cells
mean(pred_ours[comb_A == 1]) # 1

pred_ours[comb_A == 0] # posterior_mean for unobserved cells
mean(pred_ours[comb_A == 0]) # 0.894
