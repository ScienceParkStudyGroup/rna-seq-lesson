suppressPackageStartupMessages(library(tidyverse))
set.seed(1234)

#######################
# Statistical refresher
#######################
set.seed(1234) # to have the same graph everytime you run this code
xp_normal_conditions <- tibble(
  expression = rnorm(             # randomly sample numbers from a normal distribution
    n = 1000,                     # number of drawings
    mean = 2,                     # mean of the normal distribution
    sd = 0.1),                    # standard deviation of the normal distribution
  condition = "normal"            # used later on for data frame row binding with heat stress
)
head(xp_normal_conditions)


p2 <- ggplot(xp_normal_conditions, aes(x = expression, fill = condition)) +
  ggtitle("Distribution of HSF2A expression levels in 1000 samples") +
  theme(legend.position = "none") +
  geom_density(color = "black", alpha = 0.5) 
p2

### heat stress distribution drawing #1
xp_heat_stress <- tibble(expression = rnorm(n = 1000, mean = 4, sd = 0.5),
                         condition = "heat stress")

xp = bind_rows(xp_normal_conditions, xp_heat_stress) # get a peek with head(xp) and tail(xp)

p3 <- ggplot(xp, aes(x = expression, fill = condition)) +
  ggtitle("Distribution of HSF2A expression levels in 1000 samples") +
  geom_density(color = "black", alpha = 0.5) 
p3


xp_heat_stress <- tibble(expression = rnorm(n = 1000, mean = 4, sd = 0.5),
                         condition = "heat stress")

xp = bind_rows(xp_normal_conditions, xp_heat_stress) # get a peek with head(xp) and tail(xp)

### More spread normal condition
xp_normal_conditions_more_spread <- tibble(expression = rnorm(n = 1000, mean = 2, sd = 0.5),
                                           condition = "normal")

xp2 = bind_rows(xp_normal_conditions_more_spread, xp_heat_stress) # heat stress values are kept unchanged

p4 <- ggplot(xp2, aes(x = expression, fill = condition)) +
  ggtitle("Distribution of HSF2A expression levels in 1000 samples") +
  geom_density(color = "black", alpha = 0.4) + 
  scale_x_continuous(limits = c(1,6))
p4

###### Only 3 data points/biol. replicate per condition

xp_normal_conditions_only_three_points <- tibble(expression = rnorm(n = 3, # change to 3 instead of 1000
                                                                    mean = 2,
                                                                    sd = 0.5),
                                                 condition = "normal")

xp_heat_stress_only_three_points <- tibble(expression = rnorm(n = 3, 
                                            mean = 4, 
                                            sd = 0.5),
                         condition = "heat stress")


xp3 = bind_rows(xp_normal_conditions_only_three_points, 
                xp_heat_stress_only_three_points) 

p4 <- ggplot(xp3, aes(x = expression, fill = condition)) +
  ggtitle("Distribution of HSF2A expression levels with n=3 samples") +
  geom_histogram(color = "black", alpha = 0.4) + 
  scale_x_continuous(limits = c(1,6)) 
p4
