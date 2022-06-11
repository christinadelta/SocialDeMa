library(dplyr)
library(mosaic)
library(ggplot2)
library(readxl)

prepro_data = read.csv("githubstuff/rhul_stuff/SocialDeMa/economic_phase2_data.csv", 
                header = F)

prepro_beads = read.csv("githubstuff/rhul_stuff/SocialDeMa/beads_blockdata.csv", 
                        header = F)

colnames(prepro_beads) = c("subject", "block", "trialNo", "urntype", "draws", "response", 
                           "accuracy","rate", "condition", "balance")

colnames(prepro_data) = c("trialNo", "item", "samples", "price", "rank",
                           "reward", "balance")

# sampling mean (economic best choice task)
mean_sampling = mean(prepro_data$samples)

# mean sampling based on ranks
mean_rank1 = mean(prepro_data$samples[prepro_data$rank==1])
mean_rank2 = mean(prepro_data$samples[prepro_data$rank==2])
mean_rank3 = mean(prepro_data$samples[prepro_data$rank==3])

# ----------------------------------
# get sampling means of beads draws 
mean_easycond = mean(prepro_beads$draws[prepro_beads$condition==1])
mean_diffcond = mean(prepro_beads$draws[prepro_beads$condition==2])


# --------------
contracts = read_excel("githubstuff/rhul_stuff/SocialDeMa/experiments/excel_files/economic_best_choice.xls")
# make histogram with the prices
hist(contracts$price)

