source('analyze EBMH.R')
source('fun to plot EBMH.R')


length(unique(antidep$Study_No)) # number of studies
sum(antidep$No_randomised) # number of participants
nrow(antidep) # number of dose arms

# the maximum probability and the  dose with the maximum
plotdata2 = plotdata.fun(drma = rcs_pooled1,
                         data = antidep,
                         knots=knots)
plotdata2[which.max(plotdata2$prob),]

#  1-study - linear 
round(lin_1study$coefficients,4)
round(confint(lin_1study), 4)
exp(round(lin_1study$coefficients,4)*10)-1 # for 10-fold increase in dsoe, OR raise by 

# 1-study - RCS
rcs_1study$coefficients
round(confint(rcs_1study), 4)

# multi-study - RCS - 2stage
length(unique(antidep_2stage$Study_No))

round(rcs_pooled2$coefficients,4)

lb1 <- round(rcs_pooled2$coefficients[1]-1.96*sqrt(rcs_pooled2$vcov[1,1]),4)
ub1 <- round(rcs_pooled2$coefficients[1]+1.96*sqrt(rcs_pooled2$vcov[1,1]),4)

lb2 <- round(rcs_pooled2$coefficients[2]-1.96*sqrt(rcs_pooled2$vcov[2,2]),4)
ub2 <- round(rcs_pooled2$coefficients[2]+1.96*sqrt(rcs_pooled2$vcov[2,2]),4)


# multi-study - RCS - 1stage
round(rcs_pooled1$coefficients,4)
round(confint(rcs_pooled1),4)


# placebo effect
antidep_p <- antidep[antidep$Drug=='placebo',]
r <- sum(antidep_p$Responders)
n <- sum(antidep_p$No_randomised)
r/n

# p-value - 2stage
summary(rcs_pooled2)


# VPC at 20
df <- antidep[!is.na(antidep$selogOR),]
df$vpc <- vpc(rcs_pooled1)
min(df$vpc[df$hayasaka_ddd==20])*100
max(df$vpc[df$hayasaka_ddd==20])*100
max(df$vpc)*100


