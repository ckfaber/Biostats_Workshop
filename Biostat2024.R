# ########################################################################
#  Biostatistics Workshop on R/RStudio at Barrow Neurological Institute
#  Topic: Common statistical Methods with R/RStudio 
#  When: May 3rd, 2024
#  Writer: Wonsuk Yoo
# ########################################################################
#
# Install packages: ######################################################
install.packages("tidyverse")
install.packages("dplyr")

# Loading the packages:
library(tidyverse)
library(dplyr)

# ### Set-up directory/folder  ###########################################
getwd()
setwd("C:/training/biostat2024/")   # designate your directory
getwd()
# Save all datasets and R-codes in this directory

# #######################################################################
# Data Manipulation using "dplyr" 
# #######################################################################
library(dplyr)
head(starwars)
dim(starwars)

# filter()
starwars %>% filter(sex =='male', height > 170)
starwars2 <- starwars %>% filter(sex =='male', height > 170)
dim(starwars2)    

# summarise()
starwars %>% summarise(height = mean(height, na.rm = TRUE))
starwars %>% group_by(sex) %>% summarise(height = mean(height, na.rm = TRUE))

# select()
starwars3 <- starwars %>% select(hair_color, skin_color, eye_color)
dim(starwars3)

# mutate()
starwars4 <- starwars %>% mutate(height_m = height/100, BMI = mass / (height_m^2)) 
head(starwars4)


# ######################################################################
# Data Importing from CSV file   
# ######################################################################
library(tidyverse)
dat <- read_csv('nhefs_del.csv') 
head(dat); dim(dat); names(dat); str(dat)

dat <- dat %>% 
       mutate(qsmk = as.factor(qsmk),
              sex = as.factor(sex), 
              exercise = as.factor(exercise))
str(dat)


# ######################################################################
# 1. Statistical Methods in Mean Difference among Independent Groups
# ######################################################################

# One-sample t-test:  ##################################################
#
# Normality Check: QQ plot (or quantile-quantile plot) 
qqnorm(dat$wt82, pch = 1, frame = FALSE)
qqline(dat$wt82, col = "steelblue", lwd = 2)
#
# Normality Check: Shapiro-Wilk test
shapiro.test(dat$wt82)

# One-sample t-test:
t.test(dat$wt82, mu=76)

#  Independent two-sample t-tests:   ###################################

# Baseline characteristics by "qsmk":
library("dplyr")  
dat %>% summarise(
        count = n(),
        mean_B = mean(wt71, na.rm = TRUE),
        sd_B = sd(wt71, na.rm = TRUE),
        mean_E = mean(wt82, na.rm = TRUE),
        sd_E = sd(wt82, na.rm = TRUE))
dat %>% group_by(qsmk) %>%
        summarise(
        count = n(),
        mean_B = mean(wt71, na.rm = TRUE),
        sd_B = sd(wt71, na.rm = TRUE),
        mean_E = mean(wt82, na.rm = TRUE),
        sd_E = sd(wt82, na.rm = TRUE))  

# Homogeneity of variance test: 
library(car)
# Levene's test
leveneTest(y=dat$wt82, group=dat$qsmk)  #, data=datdel)  
# Barlett's test
bartlett.test(wt82~qsmk, data = dat)

# two-sample independent t-tests:
t.test(wt82~qsmk, var.equal=TRUE, dat,conf.level=0.95) # when equal vars
t.test(wt82~qsmk, dat,conf.level=0.95)  # when unequal vars

#  ANOVA test of wt82 on exercise status:  #############################
#
# Descriptive stats of wt82 on exercise:
#
dat %>% group_by(exercise) %>%
        summarise(count = n(),
        mean_wt82 = mean(wt82, na.rm = TRUE),
        sd_wt82 = sd(wt82, na.rm = TRUE),
        mean_wt71 = mean(wt71, na.rm = TRUE),
        sd_wt71 = sd(wt71, na.rm = TRUE))

# Homogeneity of variance test: 
library(car)
#Levene's test 
leveneTest(y=dat$wt82, group=dat$exercise)   #"car" package

# ANOVA test
#
aov_del <- aov(wt82 ~ exercise, data = dat)
summary(aov_del)
#

# Multiple comparison tests
#
# Tukey HSD test: 
library(multcomp)
TukeyHSD(aov_del)
plot(TukeyHSD(aov_del))
#

# ######################################################################
# 2. Statistical Methods in Mean Change between baseline and endpoint
# ######################################################################

# Paired t-test on weights between 1982 and 1971: ######################
#
pairtd <- t.test(dat$wt71, dat$wt82, paired = TRUE) 
pairtd

#  Analysis of Covariance (ANCOVA): ####################################
#
# Study aim: To test the effect of change in weight between 1971
#            1982 based on the status of quit smoking after 
#            adjusting for confounding factors. 
# Potential confounding factors of age, sex, exercise, sbp

# [1] fit the ANCOVA model using "aov"
aovfit1 <- aov(wt82 ~ wt71+qsmk+sex+exercise, data = dat)
summary(aovfit1)

# [2] fit using "lm"
lmmod1 <- lm(wt82 ~ wt71+qsmk+sex+exercise, data = dat)
summary(lmmod1)


# ######################################################################
#  Non-parametric Statistical Methods in t-tests and ANOVA
# ######################################################################

# Wilcoxon rank sum test: #####
#
# This is a non-parametric version of independent t-tests.
# Independent t-test:
t.test(wt82~qsmk, var.equal=TRUE, dat, conf.level=0.95)
# Unpaired two-sample Wilcoxon test:---
wilcox.test(dat$wt71, dat$wt82)

#  Wilcoxon signed rank test (one sample):  
#
# This is a non-parametric version of paired t-tests.
# Paired t-tests:
t.test(dat$wt71, dat$wt82, paired = TRUE) 
# (Paired sample) Wilcoxon signed rank test:-----
wilcx = wilcox.test(dat$wt71, dat$wt82, paired = TRUE, exact=FALSE)
print(wilcx)

# Kruskal-Wallis test:  
# This is a non-parametric version of oneway ANOVA tests.


# ######################################################################
# 3. Statistical Methods in Proportion Differences from Frequency Table
# ######################################################################

# Study aim:---------------------
# Let's assume that we are interested in analyzing the association 
# between two categorical variables: "qsmk" and "exercise"

# Construction of Tables for analysis:  ################################
table(dat$exercise)
table(dat$exercise, dat$qsmk)

# Chi-square tests: (from data.frame):  ################################
#
qsmkexer <- table(dat$exercise, dat$qsmk)
qsmkexer
chisq.test(qsmkexer)

# Chi-square tests: Measure of Association:  
#
library(grid)
library(vcd)
assocstats(qsmkexer)

# Fisher Exact Tests:  ##################################################
#
fisher.test(qsmkexer)  # conf.int = TRUE)

# Logistic Regression Analysis:  ########################################
#
# The logistic regression is a kind of linear regression with 
# binary outcomes. Is this possible? 
# The scatter plot of sbp on qsmk is 

# fit the full model:-----
fit1 <- glm(qsmk ~ age+income+sbp+exercise, 
            data=dat, family="binomial")
summary(fit1)

# fit the reduced model:-----
fit2 <- glm(qsmk ~ sex+age+exercise+wt71, 
            data=dat, family="binomial")
summary(fit2)

################ End of Code  ###########################################






