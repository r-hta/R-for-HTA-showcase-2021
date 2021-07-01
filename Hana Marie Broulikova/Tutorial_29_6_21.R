######################### Trial-based cost-effectiveness analysis using R: A Tutorial #########################
# script authors HM Broulikova, A Jornada Ben, M El Alili, JM van Dongen, JE Bosmans
# script for an analysis based on the seemingly unrelated regression
# version June 29, 2021

#### Libraries ####
library(haven)
library(readxl)
library(tidyverse)
library(mice)
library(boot)
library(systemfit)
library(data.table)

#### Data preparation ####
# Import the data
dataset <- read_csv("dataset.csv")

# variables codebook: 
# id = participant number, trt = treatment arm (0 for control), age = age, gen = gender (0 for man), smok = smoking (0 if not), com = comorbidity (0 if not), uT0 = baseline utility, cT0 = baseline costs, uT1/cT1-uT4/cT4 values for the respective follow-up. Lenght of the follow-up 1 year, measurements each three months

# prevent scientific notation
options(scipen = 100)

#### Multiple imputation ####
# To be also analyzed before prediction matrix is developed: share of missing observations, test variables related to missingness: Kruskall-Wallis, logistic regression

# check for missing values
p_missing <- unlist(lapply(dataset, function(x) sum(is.na(x))))
sort(p_missing[p_missing >= 0], decreasing = TRUE)

# create two lists sorted according to the treatment arm
dataset <- group_by(dataset, trt) %>% 
  group_split() 

# prediction matrix
Prediction_matrix <- read_excel("PredM.xlsx") %>% 
  select(- ...1)

row.names(Prediction_matrix) <- c("id",
                                  "trt",
                                  "age",
                                  "gen",
                                  "educ",
                                  "smok", 
                                  "com",
                                  "uT0",
                                  "cT0",
                                  "uT1",
                                  "uT2",
                                  "uT3",
                                  "uT4",
                                  "cT1",
                                  "cT2",
                                  "cT3",
                                  "cT4"
                                  )

Prediction_matrix <- as.matrix(Prediction_matrix)

# imputation method (check/set)
mi <- lapply(dataset, mice, predictorMatrix = Prediction_matrix, m = 1, maxit = 0, print = T)
mi[[1]][["method"]]
# needed, prediction method can be changed for individual variables using the code below and then specified within mice
# Method["var_name"] <- "pmm"

mi_time <- Sys.time()
set.seed(1)
mi <- lapply(dataset, mice, predictorMatrix = Prediction_matrix,  m = 5, maxit = 20, print = T) 
#mi <- lapply(mi, mice, predictorMatrix = Prediction_matrix, method = Method, m = 5, maxit = 20, print = T)
# check for problems with imputation
mi[[1]][["loggedEvents"]]
mi[[2]][["loggedEvents"]]
Sys.time() - mi_time # calculate time to run the analusis by extracting p from the current system time  

# MI data preparation
# join data from both treatment arms (results in a list with 5 imputed datasets, each with 300 observations)
mi_1 <- mi[[1]]
mi_2 <- mi[[2]]
impds <- rbind(mi_1, mi_2)

# extracts the all imputed datasets (results in a dataset with 1500 observations)
impds <- complete(impds, action = "long", include = F) 
# if include = T also original dataset is included, thus number of datasets is m+1

# calculating outcomes
impds$QALY <- ((0.5*(impds$uT0+impds$uT1)*365/4)
               + ((0.5*(impds$uT1+impds$uT2)*365/4)
               + (0.5*(impds$uT2+impds$uT3)*365/4)
               + (0.5*(impds$uT3+impds$uT4)*365/4)))/365
  
impds$Tcosts <- impds$cT1 + impds$cT2 + impds$cT3 + impds$cT4 

#### Bootstrapping ####
# create lists to again spit by treatment arm
impds <- group_by(impds,.imp) %>% 
  group_split() 

# function to execute seemingly unrelated regression for each bootstrapped sample
fsur <- function(impds, i){
  dscc2 <- impds[i,] 
  r1 <- Tcosts ~ trt + uT0 + age + gen + educ + smok + com 
  r2 <- QALY  ~ trt + uT0 + age + gen + educ + smok + com
  fit <- with(data = impds, exp = systemfit(list(costreg = r1, effectreg = r2), "SUR", data=dscc2))
  betas <- fit$coefficients
  return(c(betas[2], betas[10])) # 
}

boot_time <- Sys.time()   
set.seed(2)
# bootstrap imputed datasets and get SUR results 
bootce <- lapply(impds, function(x) boot(data=x, statistic=fsur, R=5000)) 
# bootstrap confidence interval calculations for costs (also used to get SUR estimate t0)
bootcic <- lapply(bootce, function(x) boot.ci(boot.out = x, type = c("bca"), index = 1)) 
# bootstrap confidence interval calculations for effects (also used to get SUR estimate t0)
booteic <- lapply(bootce, function(x) boot.ci(boot.out = x, type = c("bca"), index = 2)) 
Sys.time() - boot_time

#### Analysis ####
# create the dataset to be analyzed
# this is the dataset that contains final data for cea analysis (5000 SUR results per imputed dataset)
boot_df <- lapply(bootce, function(x) as.data.frame((x[["t"]]))) 

# the observed value of the statistic; i.e. estimated difference in costs (1 number per imputed dataset); can crosschecked with sur coefficient
cost_diff <- lapply(bootcic, function(x) x[["t0"]]) 

# the observed value of the statistic; i.e. estimated difference in effects (1 number per imputed dataset)
effect_diff <- lapply(booteic, function(x) x[["t0"]]) 

# extractes the confidence intevals estimated by bca (to be used to estimate CI for costs)
bcacost <- lapply(bootcic, function(x) x[["bca"]]) 

# perform SUR to get variance
# define regressions to extract variance-covariance matrix
r1 <- Tcosts ~ trt + uT0 + age + gen + educ + smok + com 
r2 <- QALY  ~ trt + uT0 + age + gen + educ + smok + com

sur <- lapply(impds, function(x) {
  systemfit(list(costreg = r1, effectreg = r2), "SUR", data=x)})

# extract variance-covariance matrix for the effects (to be used to estimate CI for effects)
varcov <- lapply(sur, function(x) x[["coefCov"]])
var_e <- lapply(varcov, function(x) x[10,10])
var_e <- lapply(var_e, setNames, c("var_effects"))

# binding results in one dataset
boot_df <- mapply(cbind, boot_df, cost_diff, effect_diff, bcacost, var_e, SIMPLIFY=F)
boot_df <- setNames(boot_df, c("mi_1", "mi_2", "mi_3", "mi_4", "mi_5")) 
boot_df <- lapply(boot_df, setNames, c("bootcost_diff", "booteffect_diff",
                                       "obscost_diff", "obseffect_diff", 
                                       "drop_1", "drop_2", "drop_3", 
                                       "LL_cost", "UL_cost", 
                                       "var_effects"))
boot_df <- lapply(boot_df, function(x) select(x, - drop_1, - drop_2, - drop_3))
boot_df <- reduce(boot_df, bind_rows)
observed <- boot_df[!duplicated(boot_df$obscost_diff,fromLast=F),] %>% 
  select(-bootcost_diff, -booteffect_diff)

# pooling results
# cost and effect differences
effect_diff_pooled <- mean(observed$obseffect_diff)
cost_diff_pooled <- mean(observed$obscost_diff)

# lower and upper level limits for effects using Rubin's rules
Za = 1.95996
W <- mean(observed$var_effects)
B_diff <- (observed$obseffect_diff-effect_diff_pooled)^2
B_sum <- sum(B_diff)
B <- (1/(10-1))*B_sum
T_ <- W + (1+(1/10))*B
seT <- sqrt(T_)
LL_effect_pooled <- effect_diff_pooled - (Za*seT)
UL_effect_pooled <- effect_diff_pooled + (Za*seT)

# lower and upper level limits for costs using bca because of the skewness
LL_cost_pooled <- mean(boot_df$LL_cost)
UL_cost_pooled <- mean(boot_df$UL_cost)

# ICER
ICER <- cost_diff_pooled/effect_diff_pooled

# distribution of the results in the cost-effectiveness plane
boot_df <- mutate(boot_df, ICER = bootcost_diff/booteffect_diff,
                  NE = ifelse(ICER > 0 & bootcost_diff > 0, 1, 0),
                  SE = ifelse(ICER < 0 & bootcost_diff < 0, 1, 0),
                  NW = ifelse(ICER < 0 & bootcost_diff > 0, 1, 0),
                  SW = ifelse(ICER > 0 & bootcost_diff < 0, 1, 0))
NW <- sum(boot_df$NW)/nrow(boot_df)*100 
SW <- sum(boot_df$SW)/nrow(boot_df)*100 
NE <- sum(boot_df$NE)/nrow(boot_df)*100 
SE <- sum(boot_df$SE)/nrow(boot_df)*100 
# check has to be equal to 100
check <- NW + SW + NE + SE

#### Figures ####
# cost-effectiveness plane
qalys <-c(-0.1, 0.1)
data <- as.data.frame(qalys) %>% 
  mutate(costs22 = qalys*22000, 
         costs34 = qalys*34000,
         cartesian = 0, 
         cartesian_x = c(-10000,50000))

scatter <- ggplot(boot_df, aes(x=booteffect_diff, y=bootcost_diff)) +
  geom_point(color = "black", size = 1) +
  labs(x = "Incremental effects (QALYs)",
       y = "Incremental costs (€)",
       color = "WTP (€)",
       title = "2a")+
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  geom_segment(aes(x = cartesian, y = cartesian, xend = qalys, yend = 0), data = data) +
  geom_segment(aes(x = cartesian, y = cartesian, yend = cartesian_x, xend = 0), data = data) + 
  theme_classic(base_size = 20)
scatter

# the code below saves the figures
#ggsave(file = "scatter.tiff", # save name
#       plot = scatter, # object to save
#       width = 30,
#       height = 20,
#       units = "cm",
#       dpi = 300) # resolution

# CEAC
results <- data.table(boot_df)
CEAC <- data.table()
for (wtp in seq(0, 50000, 1000)) {
  propPos <- results[, mean((wtp * booteffect_diff - bootcost_diff) > 0)] # This is net benefit for everyone 
  #return(propPos)
  CEAC<- rbind(CEAC, data.table(wtp, propPos))
}

CEAC <- mutate(CEAC, Scenario = "treatment vs. control")

CEACs_plot <- ggplot(CEAC, aes(x = wtp, y = propPos, color = Scenario)) + 
  geom_line() + 
  scale_color_manual(values=c("black")) +
  labs(title="2b") +
  xlab("Willingness to pay threshold (€/QALY)") +
  ylab("Probability of cost-effectiveness") +
  scale_y_continuous(limits = c (0 , 1), breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))+
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))
CEACs_plot

#ggsave(file = "CEACs_plot.tiff", # save name
#       plot = CEACs_plot, # object to save
#       width = 25,
#       height = 15,
#       units = "cm",
#       dpi = 300) # resolution

