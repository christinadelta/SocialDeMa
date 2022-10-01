# library
library(ggplot2)
library(readxl)
library(ggpubr)

# beads_all_avdraws = read.csv("~/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/behav/beads_all_avdraws.csv", 
#                             header=FALSE)


all_agent_draws = read_excel("Desktop/figures_beads/all_agent_draws.xls")

all_agent_acc = read_excel("Desktop/figures_beads/all_agent_acc.xls")

model_fitting_cs = read_excel("Desktop/figures_beads/model_fitting_cs.xlsx")
model_fitting_ll = read_excel("Desktop/figures_beads/model_fitting_ll.xlsx")

# grouped boxplot -- draws 
ggplot(all_agent_draws, aes(x=probability2, y=draws, fill=agent2)) + 
  geom_boxplot()

# grouped boxplot -- accuracy
ggplot(all_agent_acc, aes(x=probability2, y=acc, fill=agent2)) + 
  geom_boxplot()


# make grouped bar plots 
# group by: agent type
ggbarplot(
  all_agent_draws, x = "probability", y = "draws", 
  add = c("mean_sd", "jitter"), 
  add.params = list(shape = "agent_type"),
  fill= "agent_type", palette = c("#807F7F", "#BF504D"),
  position = position_dodge(0.8)
)

# try barplots for acc
ggbarplot(
  all_agent_acc, x = "probability", y = "acc", 
  add = c("mean_sd", "jitter"), 
  add.params = list(shape = "agent_type"),
  fill= "agent_type", palette = c("#807F7F", "#BF504D"),
  position = position_dodge(0.8)
)

# make bar plots for model_cd and ll parameters
# Create a simple bar plot
ggbarplot(model_fitting_cs, x = "cond", y = "cs",
          label = TRUE, label.pos = "out")


