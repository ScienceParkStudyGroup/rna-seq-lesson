# Data analysis for the life sciences
# https://github.com/genomicsclass/labs
library("tidyverse")

# This notebook https://github.com/genomicsclass/labs/blob/master/inference/random_variables.Rmd
dat <- read.csv("02.lectures/femaleMiceWeights.csv") 
head(dat)

ggplot(dat, aes(x = Diet, y = Bodyweight, fill=Diet)) + 
  geom_point() + 
  geom_boxplot() + 
  geom_jitter(width = 0.05) +
  scale_x_discrete("Diet", labels=c("Normal diet", "High-fat diet")) +
  theme(text = element_text(size=16))
ggsave(filename = "02.lectures/02_stats_femaleMiceWeights.png")
