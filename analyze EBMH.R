# Run dose-effect meta-analysis: 
# 1. one-study analysis (with linear and RCS)
# 2. multi-studies analysis: 1stage and 2stage (with RCS)


library(rms) # rcs()
library(dosresmeta) # dosresmeta()
library(dplyr)
library(meta) # metaprop


# ---------- load data an prepare ----------

# load and exclude single arm studies
mydata <-  read.csv('DOSEmainanalysis.csv')
antidep=mydata[mydata$exc==F,]
# add OR
source('fun to analyze EBMH.R') # has createORreference.fun()
antidep$studyid <- as.numeric(as.factor(antidep$Study_No))
antidep <- antidep%>%arrange(Study_No,hayasaka_ddd)
antidep$nonResponders <- antidep$No_randomised- antidep$Responders
logORmat <- sapply(unique(antidep$studyid),function(i) createORreference.fun(antidep$Responders[antidep$studyid==i],antidep$No_randomised[antidep$studyid==i]),simplify = FALSE)
logORmat <- do.call(rbind,logORmat)
antidep$logOR <- c(logORmat[,1])
antidep$selogOR <- c(logORmat[,2])

# knots
knots= quantile(antidep$hayasaka_ddd[antidep$hayasaka_ddd!=0],c(0.10,0.50,0.90))

# ---------- 1.one-study analysis ----------
study_87 <- antidep[antidep$Study_No=='87',] 

# linear
lin_1study <- dosresmeta(formula=logOR~hayasaka_ddd, 
                         id=Study_No, 
                         type=type,
                         cases=Responders,
                         n=No_randomised,
                         se=selogOR,
                         data=study_87,
                         method = 'reml')

summary(lin_1study)

# RCS
rcs_1study <- dosresmeta(formula=logOR~rcs(hayasaka_ddd,knots), 
                         id=Study_No, 
                         type=type,
                         cases=Responders,
                         n=No_randomised,
                         se=selogOR,
                         data=study_87,
                         method = 'reml')
summary(rcs_1study)
waldtest(b=coef(rcs_1study)[2], 
               Sigma=vcov(rcs_1study)[2,2],
               Terms=1) # wald test for spline coefficient
# ---------- 1.multi-study analysis ----------
# 1-stage
rcs_pooled1 <- dosresmeta(formula=logOR~rcs(hayasaka_ddd,knots), 
                          proc="1stage",
                          id=Study_No, 
                          type=type,
                          cases=Responders,
                          n=No_randomised,
                          se=selogOR,
                          data=antidep,
                          method = 'reml')
summary(rcs_pooled1)

# 2-stage
# include studies with at least 3 arms
studies_2arm <- unique(antidep$Study_No)[table(antidep$Study_No)<3]
antidep_2stage <- antidep[!antidep$Study_No%in%studies_2arm,]

rcs_pooled2 <- dosresmeta(formula=logOR~rcs(hayasaka_ddd,knots), 
                          proc="2stage",
                          id=Study_No, 
                          type=type,
                          cases=Responders,
                          n=No_randomised,
                          se=selogOR,
                          data=antidep_2stage,
                          method = 'reml')
summary(rcs_pooled2)
# placebo effect - meta-analysis
antidep_p <- antidep[antidep$Drug=='placebo',]
meta_pl<-metaprop(event=Responders, 
                       n=No_randomised, 
                       data=antidep_p[!(is.na(antidep_p$Responders)|is.na(antidep_p$No_randomised)),], 
                       studlab=Study_No,
                       comb.fixed = FALSE)
meta_pl# proportion
pl_eff <- exp(meta_pl$TE.random)/(1+exp(meta_pl$TE.random))
# r <- sum(antidep_p$Responders)
# n <- sum(antidep_p$No_randomised)
# r/n

