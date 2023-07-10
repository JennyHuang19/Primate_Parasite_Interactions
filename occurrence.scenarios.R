# creating OBS occurrence scenarios.
# Primates: Habitat categories
all_X2 <- read.csv("new_all_X2.csv")
colnames(all_X2)

# Keeping track of the unique primates and parasite species in the study and their
# number. Also keeping track of the unique studies.
uni_primates <- unique(traits_and_studies$HostCorrectedName)
uni_parasites <- unique(traits_and_studies$ParasiteCorrectedName)
uni_studies_lst <- unique(traits_and_studies$Citation)
uni_locations_lst <- unique(traits_and_studies$LocationName)
nM <- length(uni_primates)
nP <- length(uni_parasites)
nS <- length(uni_studies_lst)
# nL <- length(uni_locations_lst) # 192 (so not much overlap between studies)
# nR <- length(uni_realms_lst)
# nH <- length(uni_habitats_lst)


## 1.
# Primates: Present in study: 1. Present at location: 0.85. 2. Present in realm: 0.5. Present in habitat: 0.25. Not present in habitat: 0.1.
# Parasites: Present in study: 1. Present at location: 0.75. Not present at location: 0.25.
## 2. (high)
# Primates: Present in study: 1. Present at location: 0.85. 2. Present in realm: 0.65. Present in habitat: 0.35. Not present in habitat: 0.25.
# Parasites: Present in study: 1. Present at location: 0.85. Not present at location: 0.25.
## 3. (low)
# Primates: Present in study: 1. Present at location: 0.75. 2. Present in realm: 0.5. Present in habitat: 0.20. Not present in habitat: 0.0.
# Parasites: Present in study: 1. Present at location: 0.75. Not present at location: 0.1.

# (May 12th Replace 0 with 0.5:) obs_OP2 <- array(0, dim = c(nP, nS))
obs_OB2.1 <- array(0.1, dim = c(nM, nS))
obs_OP2.1 <- array(0.25, dim = c(nP, nS))

dimnames(obs_OB2.1) <- list(uni_primates, uni_studies_lst)
dimnames(obs_OP2.1) <- list(uni_parasites, uni_studies_lst)

obs_OB2.2 <- array(0.25, dim = c(nM, nS))
obs_OP2.2 <- array(0.25, dim = c(nP, nS))

dimnames(obs_OB2.2) <- list(uni_primates, uni_studies_lst)
dimnames(obs_OP2.2) <- list(uni_parasites, uni_studies_lst)

obs_OB2.3 <- array(0.0, dim = c(nM, nS))
obs_OP2.3 <- array(0.1, dim = c(nP, nS))

dimnames(obs_OB2.3) <- list(uni_primates, uni_studies_lst)
dimnames(obs_OP2.3) <- list(uni_parasites, uni_studies_lst)

# start with lowest to highest probabilities.
# A: HABITAT.
# (may 28th: most primates belong to the forest)
# sum(all_X2$habitat1 == "forest") # 93

### June 1st 
traits_and_studies$study_habitats <- NA
# View(traits_and_studies)
# ---
for (ss in 1 : nS) {
  # this study
  wh_study <- uni_studies_lst[ss]
  ## find the habitats associated with the study.
  # list of all primates in the study.
  these_monkeys <- unique(traits_and_studies$HostCorrectedName[traits_and_studies$Citation == wh_study]) # primates that have been studied here.
  # list of all habitats of all primates in the study.
  these_habitats <- list()
  for (monkey in these_monkeys){ 
    # print(monkey)
    these_habitats <- append(these_habitats, all_X2$habitat1[all_X2$HostCorrectedName == monkey])
    # print(these_habitats)
    }
  # get the habitat of the study.
  wh_hab <- unique(unlist(these_habitats)) 
  
  ### June 1st
  traits_and_studies$study_habitats[traits_and_studies$Citation == wh_study] <- list(wh_hab)
  # ---
  ## find all the primates Present in habitat:
  these_primates <- list()
  for (hbt in wh_hab){
    these_primates <- append(these_primates, all_X2$HostCorrectedName[all_X2$habitat1 %in% wh_hab]) # sum(all_X2$habitat_category %in% wh_hab) >= 1. all_X2$habitat_category %in% wh_hab (list of T/F -> one T/F if present)
    these_primates <- append(these_primates, all_X2$HostCorrectedName[all_X2$habitat2 %in% wh_hab])
    these_primates <- append(these_primates, all_X2$HostCorrectedName[all_X2$habitat3 %in% wh_hab])
    
    these_parasites <- append(these_parasites, traits_and_studies$ParasiteCorrectedName[traits_and_studies$study_habitats %in% wh_hab])
    }
  wh_primates <- which(uni_primates %in% these_primates)

  obs_OB2.1[wh_primates, ss] <- 0.25 # indicator for primate in habitat. 
  obs_OB2.2[wh_primates, ss] <- 0.35
  obs_OB2.3[wh_primates, ss] <- 0.20
  
  
}
# A2. REALM.
### June 1st
traits_and_studies$study_realm <- NA
# ---
for (ss in 1 : nS) {
  # this study
  wh_study <- uni_studies_lst[ss]
  ## find the realms associated with the study.
  # list of all primates in the study.
  these_monkeys <- unique(traits_and_studies$HostCorrectedName[traits_and_studies$Citation == wh_study]) # primates that have been studied here.
  # list of all habitats of all primates in the study.
  ### Realms list of all studies.
  these_realms <- list()
  for (monkey in these_monkeys){ 
    these_realms <- append(these_realms, all_X2$realm[all_X2$HostCorrectedName == monkey]) 
    }
  # get the realm of the study.
  wh_realm <- unique(unlist(these_realms))
  ### June 1st - creating new study_realm column in traits_and_studies
  # traits_and_studies$study_realm[traits_and_studies$Citation == wh_study] <- wh_realm
  # ---
  ## find all the primates present in realm:
  these_primates <- all_X2$HostCorrectedName[all_X2$realm %in% list(wh_realm)]
  wh_primates <- which(uni_primates %in% these_primates)
  
  ## find all parasites
  these_parasites <- traits_and_studies$ParasiteCorrectedName[traits_and_studies$study_realm == wh_realm]
  wh_parasites <- which(uni_parasites %in% these_parasites)

  obs_OB2.1[wh_primates, ss] <- 0.5 # indicator for primate in realm 
  obs_OB2.2[wh_primates, ss] <- 0.65
  obs_OB2.3[wh_primates, ss] <- 0.5
}

# B. LOCATION.
for (ss in 1 : nS) {
  # this study
  wh_study <- uni_studies_lst[ss]
  # location of study
  wh_loc <- traits_and_studies$LocationName[traits_and_studies$Citation == wh_study]
  # Present in location:
  these_primates <- unique(traits_and_studies$HostCorrectedName[traits_and_studies$LocationName == wh_loc]) # primates that have been studied at this location.
  these_parasites <- unique(traits_and_studies$ParasiteCorrectedName[traits_and_studies$LocationName == wh_loc]) # parasites that have been studied at this location.
  wh_primates <- which(uni_primates %in% these_primates)
  wh_parasites <- which(uni_parasites %in% these_parasites)
  
  obs_OB2.1[wh_primates, ss] <- 0.85 # indicator that the primate has been studied here.
  obs_OP2.1[wh_parasites, ss] <- 0.85 # indicator that the parasite has been studied here.
  
  obs_OB2.2[wh_primates, ss] <- 0.85 # indicator that the primate has been studied here.
  obs_OP2.2[wh_parasites, ss] <- 0.85 # indicator that the parasite has been studied here.
  
  obs_OB2.3[wh_primates, ss] <- 0.75 # indicator that the primate has been studied here.
  obs_OP2.3[wh_parasites, ss] <- 0.75 # indicator that the parasite has been studied here.
}
# C. STUDY.
for (ss in 1 : nS) {
  wh_study <- uni_studies_lst[ss]
  # Present in study:
  these_primates <- unique(traits_and_studies$HostCorrectedName[traits_and_studies$Citation == wh_study]) # primates that have been studied here.
  these_parasites <- unique(traits_and_studies$ParasiteCorrectedName[traits_and_studies$Citation == wh_study]) # parasites that have been studied here.
  wh_primates <- which(uni_primates %in% these_primates)
  wh_parasites <- which(uni_parasites %in% these_parasites)
  obs_OB2.1[wh_primates, ss] <- 1 # indicator that the primate has been studied here.
  obs_OP2.1[wh_parasites, ss] <- 1 # indicator that the parasite has been studied here.
  obs_OB2.2[wh_primates, ss] <- 1 # indicator that the primate has been studied here.
  obs_OP2.2[wh_parasites, ss] <- 1 # indicator that the parasite has been studied here.
  obs_OB2.3[wh_primates, ss] <- 1 # indicator that the primate has been studied here.
  obs_OP2.3[wh_parasites, ss] <- 1 # indicator that the parasite has been studied here.
  
}

count_OP2 <- table(obs_OP2.1) # view proportions of every category.
count_OP2
# 0.25  0.85     1 
# 32946   552   571 # 571 parasites studied. 552 additional parasites at locations.
count_OP2.2 <- table(obs_OP2.2) # view proportions of every category.
count_OP2.2
# 0.25  0.85     1 
# 32946   552   571 
count_OP2.3 <- table(obs_OP2.3) # view proportions of every category.
count_OP2.3
# 0.1  0.75     1 
# 32946   552   571 
count_OB2 <- table(obs_OB2.1) 
count_OB2
# 0.1  0.25   0.5  0.85     1 
# 466 13709  5434   244   328
count_OB2.2 <- table(obs_OB2.2) 
count_OB2.2
# 0.25  0.35  0.65  0.85     1 
# 466 13709  5434   244   328 
count_OB2.3 <- table(obs_OB2.3) 
count_OB2.3
# 0   0.2   0.5  0.75     1 
# 466 13709  5434   244   328

### save files..
if (save_files) {
  save(obs_OB2.1, file = paste0(data_path, 'obs_OB2.1.dat'))
  save(obs_OP2.1, file = paste0(data_path, 'obs_OP2.1.dat'))
  
  save(obs_OB2.2, file = paste0(data_path, 'obs_OB2.2.dat'))
  save(obs_OP2.2, file = paste0(data_path, 'obs_OP2.2.dat'))
  
  save(obs_OB2.3, file = paste0(data_path, 'obs_OB2.3.dat'))
  save(obs_OP2.3, file = paste0(data_path, 'obs_OP2.3.dat'))
}
# ---

# Bar plot
# 0.25  0.85     1 
# 32946   552   571
ps <- data.frame(
  group = c("None", "Location", "Study"),
  value = c(32946, 552, 571)
)

# Create the bar plot
ggplot(ps, aes(x = group, y = value)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(aes(label=value), vjust=-0.3, size=6) +
  theme_minimal() +
  labs(x = "Occurrence", y = "Count", title = "Parasite Occurrences", size=8)

mk <- data.frame(
  group = c("None", "Habitat", "Realm", "Location", "Study"),
  value = c(466, 13709, 5434, 244, 328)
) 

# Create the bar plot
ggplot(mk, aes(x = group, y = value)) +
  geom_bar(stat = "identity", fill = "black") +
  geom_text(aes(label=value), vjust=-0.3, size=6) +
  theme_minimal() +
  labs(x = "Occurrence", y = "Count", title = "Primate Occurrences")
# ---

### OLD CODE.
# habitats_lst <- all_X2 %>% 
#   filter(habitat_category != "") %>%  # filter out the nas
#   mutate(habitat_category = strsplit(habitat_category, "_")) %>% 
#   select(habitat_category) # match when: a habitat in the study (each study has a list of habitats, reliant on the primate in the study) matches a habitat for that primate.
# habitats_lst <- all_X2 %>% 
#   mutate(habitat_category = strsplit(habitat_category, "_")) %>% 
#   select(habitat_category) 
# 
# # mutate: habitats list into new column.
# habitat <- habitats_lst$habitat_category
# 
# testtest <- all_X2 %>% 
#   mutate(habitat_new = habitat) # split the element in python on _..
