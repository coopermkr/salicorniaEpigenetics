#################################
# Analysis code for Rubi et al. (in press) "Patterns of 
#   genetic and epigenetic diversity across a range 
#   expansion in the white-footed mouse (Peromyscus 
#   leucopus)". Integrative Organismal Biology.

##################
path <- "~"

in_epi <- "EDbyLocus.csv"
in_gen <- "populations.sumstats.forimport.txt"
################## 

library(lme4)
library(lmerTest)
library(multcomp)
library(plyr)

#############
# read in tables

df_epi <- read.csv(paste(path,in_epi,sep="/"), head=T)
df_gen <- read.table(paste(path,in_gen,sep="/"), sep="\t", head=T)

####################################################
# Genetic
####################################################

# Ho
Ho <- with(df_gen,lm(Obs.Het ~ Pop.ID))

Ho_KW <-kruskal.test(Obs.Het ~ Pop.ID, data=df_gen)
# Kruskal-Wallis rank sum test
# 
# data:  Obs.Het by Pop.ID
# Kruskal-Wallis chi-squared = 3434.2, df = 3, p-value <
#   2.2e-16

Ho_aov <-aov(Ho)
TukeyHSD(Ho_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Ho)
# 
# $Pop.ID
# diff         lwr           upr     p adj
# ChipPln-CheyPln -0.011569886 -0.01922692 -0.0039128549 0.0005999
# MenPln-CheyPln  -0.026548341 -0.03421203 -0.0188846558 0.0000000
# SchPln-CheyPln  -0.034981277 -0.04268279 -0.0272797667 0.0000000
# MenPln-ChipPln  -0.014978455 -0.02264214 -0.0073147697 0.0000031
# SchPln-ChipPln  -0.023411390 -0.03111290 -0.0157098806 0.0000000
# SchPln-MenPln   -0.008432936 -0.01614106 -0.0007248106 0.0254820


# pi
pi <- with(df_gen,lm(Pi ~ Pop.ID))

pi_KW <-kruskal.test(Pi ~ Pop.ID, data=df_gen)
# Kruskal-Wallis rank sum test
# 
# data:  Pi by Pop.ID
# Kruskal-Wallis chi-squared = 3455.8, df = 3, p-value <
#   2.2e-16


pi_aov <-aov(pi)
TukeyHSD(pi_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pi)
# 
# $Pop.ID
# diff         lwr          upr    p adj
# ChipPln-CheyPln -0.02649477 -0.03357345 -0.019416088 0.00e+00
# MenPln-CheyPln  -0.03851740 -0.04560223 -0.031432574 0.00e+00
# SchPln-CheyPln  -0.05048351 -0.05760331 -0.043363709 0.00e+00
# MenPln-ChipPln  -0.01202264 -0.01910747 -0.004937807 7.69e-05
# SchPln-ChipPln  -0.02398874 -0.03110854 -0.016868941 0.00e+00
# SchPln-MenPln   -0.01196610 -0.01909202 -0.004840189 9.44e-05

# see means
ddply(df_gen, .(Pop.ID), summarize,
      Ho_mean = mean(Obs.Het),
      Ho_SD = sd(Obs.Het),
      pi_mean = mean(Pi),
      pi_SD = sd(Pi))
# Pop.ID    Ho_mean     Ho_SD    pi_mean     pi_SD
# 1 CheyPln 0.12899500 0.1363564 0.14551212 0.1364210
# 2 ChipPln 0.11742511 0.1793513 0.11901735 0.1570214
# 3  MenPln 0.10244666 0.2005293 0.10699471 0.1831815
# 4  SchPln 0.09401372 0.2301936 0.09502861 0.2150159

####################################################
# Epigenetic
####################################################

# uh
uh <- with(df_epi,lm(uh ~ Pop))

uh_KW <-kruskal.test(uh ~ Pop, data=df_epi)
# Kruskal-Wallis rank sum test
# 
# data:  uh by Pop
# Kruskal-Wallis chi-squared = 7.7688, df = 3, p-value = 0.05104


uh_aov <-aov(uh)
TukeyHSD(uh_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = uh)
# 
# $Pop
# diff         lwr          upr     p adj
# 2Sch-1Men -0.002972141 -0.04024502  0.034300742 0.9969540
# 3Chi-1Men -0.012646628 -0.04991951  0.024626255 0.8192536
# 4Che-1Men -0.058111437 -0.09538432 -0.020838555 0.0003662
# 3Chi-2Sch -0.009674487 -0.04694737  0.027598395 0.9094336
# 4Che-2Sch -0.055139296 -0.09241218 -0.017866414 0.0008416
# 4Che-3Chi -0.045464809 -0.08273769 -0.008191927 0.0093862

# SD
sd <- with(df_epi,lm(sdev ~ Pop))
SD_KW <-kruskal.test(sdev ~ Pop, data=df_epi)
# Kruskal-Wallis rank sum test
# 
# data:  sdev by Pop
# Kruskal-Wallis chi-squared = 36.543, df = 3, p-value = 5.749e-08

sd_aov <-aov(sd)
TukeyHSD(sd_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = sd)
# 
# $Pop
# diff         lwr          upr     p adj
# 2Sch-1Men -0.009063518 -0.02768338  0.009556344 0.5941580
# 3Chi-1Men -0.002177216 -0.02079708  0.016442645 0.9905673
# 4Che-1Men -0.042807536 -0.06142740 -0.024187675 0.0000000
# 3Chi-2Sch  0.006886301 -0.01173356  0.025506163 0.7773657
# 4Che-2Sch -0.033744018 -0.05236388 -0.015124157 0.0000198
# 4Che-3Chi -0.040630319 -0.05925018 -0.022010458 0.0000001

df_epi$region <- ifelse(df_epi$Pop=="4Che","LP","UP")
# see means (group together all UP)
ddply(df_epi, .(region), summarize,
      uh_mean = mean(uh),
      uh_SD = sd(uh),
      sd_mean = mean(sdev),
      sd_SD = sd(sdev))
# region   uh_mean     uh_SD   sd_mean      sd_SD
# 1     LP 0.2072038 0.2000322 0.1686077 0.08740696
# 2     UP 0.2601090 0.2867041 0.2076683 0.14595192
