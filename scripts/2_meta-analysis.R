#######################################################
#
#     EU meadow birds meta-analysis
#
#######################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

#=================================  SET DIRECTORY STRUCTURE  ================================

# LOCAL
if(.Platform$OS =='windows') {
  cluster <- FALSE
  Mac <- FALSE
}

# HPCBTO
if(.Platform$OS=='unix' & Sys.getenv('USER')=='samf') {
  cluster <- TRUE
  Mac <- FALSE
  Wales <- FALSE
}

# Mac
if(.Platform$OS=='unix' & Sys.getenv('USER')=='samantha') {
  cluster <- FALSE
  Mac <- TRUE
  Wales <- FALSE
}

#### SET DIRECTORY PATHS
# # Wales HPC cluster
# if (cluster) parentwd <- c("/home/samantha.franks/")

if (cluster) parentwd <- c("/users1/samf") # BTO cluster
if (!cluster) {
  if (!Mac) parentwd <- c("C:/Users/samf/Documents/Git/eu_meadow_birds")
  if (Mac) parentwd <- c("/Volumes/SAM250GB/BTO PC Documents/Git/eu_meadow_birds")
}

scriptswd <- paste(parentwd, "scripts", sep="/")
datawd <- paste(parentwd, "data", sep="/")
outputwd <- paste(parentwd, "output", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


#=================================  LOAD DATA  ===============================

dat0 <- readRDS(paste(workspacewd, "meadow birds analysis dataset.rds", sep="/"))

mgmtvars <- c("AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# subset dataset for analysis to desired columns only
dat <- subset(dat0, select=c("reference","country","study.length","habitat","species","overall.metric","sample.size","analysis2","success",mgmtvars))


#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}

#=================================  ANALYSIS  ===============================

# Analysis evaluates whether one or several management interventions has a significantly positive outcome (i.e. is successful) compared to the baseline reference situation (which in most cases is no management)
# reference level for all interventions is 'none', so management (either applied or reduced, depending on the intervention) is evaluated against the reference
# are studies which don't apply a particular management intervention more likely to have a positive outcome for the metric being measured than studies that don't apply that intervention? 

# 0a) success of individual management types overall
# 0b) success of individual management types by species
# 1a) success of AES and nature reserves
# 1b) success of specific management interventions
# 2a) success of AES and nature reserves by species
# 2b) success of specific management interventions by species
# 3a) success of AES and nature reserves by different metrics (abundance vs reproductive success)
# 3b) success of specific management interventions by different metrics (abundance vs reproductive success)
# 4a) success of AES and nature reserves by reproductive success metrics (nest vs chick)
# 4b) success of specific management interventions by reproductive success metrics (nest vs chick)
# 5a) success of AES and nature reserves by abundance metrics (abundance vs abundance change)
# 5b) success of specific management interventions by abundance metrics (abundance vs abundance change)
# 6a) success of AES and nature reserves by habitat and species
# 6b) success of specific management interventions by habitat and species

### Nuisance explanatory variables and random effects
# ‘nuisance explanatory variables: study duration (continuous), sample size (categorical: small, medium, large), multivariate/univariate (categorical)
# Random effects: country, study id - habitat will be included later on

#------------------------------ 0a) success of individual management types -------------------------

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$species, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# subset out categories with few observations
mdat <- subset(dat, mowing!="applied" & fertpest!="applied" & predator.control!="reduced")
mdat <- droplevels(mdat)

m.ind <- list()

for (i in 1:length(mgmtvars)) {
  
  modform <- as.formula(paste("success ~ ", mgmtvars[i], " + study.length + analysis2 + sample.size + (1|reference) + (1|country)", sep=""))
  m.ind[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
}

warningmessages <- lapply(m.ind, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages

### Output model results ###

setwd(outputwd)
sink(paste("abundance_0a.txt", sep=" "))
cat("\n########==========  0a) success of individual management types ==========########\n", sep="\n")
print(lapply(m.ind, summary))
print(warningmessages)
sink()

#------------------------------ 0b) success of individual management types by species -------------------------

# subset dataset to remove dunlin and ruff, not enough data for these species for any analysis that includes species as a covariate
mdat <- subset(mdat, species!="dunlin" & species!="ruff")
mdat <- droplevels(mdat)

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(mdat$species, mdat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# problems with all 0's or all 1's in a single category level - complete separation
# try Ben Bolker's fix using the blme package
# Your intuition is exactly right. This phenomenon is called complete separation. You can find quite a lot (now that you know its name) Googling around ... It is fairly thoroughly discussed here in a general context, and here in the context of GLMMs. The standard solution to this problem is to add a small term that pushes the parameters back toward zero -- in frequentist contexts this is called a penalized or bias-corrected method. The standard algorithm is due to Firth (1993, "Bias reduction of maximum likelihood estimates" Biometrika 80, 27-38), and is implemented in the logistf package on CRAN. In Bayesian contexts this is framed as adding a weak prior to the fixed-effect parameters.
# 
# To my knowledge Firth's algorithm hasn't been extended to GLMMs, but you can use the Bayesian trick by using the blme package, which puts a thin Bayesian layer over the top of the lme4 package. Here's an example from the above-linked GLMM discussion:
# cmod_blme_L2 <- bglmer(predation~ttt+(1|block),data=newdat,
# family=binomial,
# fixef.prior = normal(cov = diag(9,4)))
# The first two lines in this example are exactly the same as we would use in the standard glmer model; the last specifies that the prior for the fixed effects is a multivariate normal distribution with a diagonal variance-covariance matrix. The matrix is 4x4 (because we have 4 fixed-effect parameters in this example), and the prior variance of each parameter is 9 (corresponding to a standard deviation of 3, which is pretty weak -- that means +/- 2SD is (-6,6), which is a very large range on the logit scale).
# 
# The very large standard errors of the parameters in your example are an example of a phenomenon closely related to complete separation (it occurs whenever we get extreme parameter values in a logistic model) called the Hauck-Donner effect.

m.ind.sp <- list()

for (i in 1:length(mgmtvars)) {
  
#   # remove additional combinations/categories of interactions/species for which there aren't enough data
#   if (mgmtvars[i]=="mowing") {
#     mdat <- subset(mdat, mowing!="applied")
#     mdat <- droplevels(mdat)
#   }
#   
#   if (mgmtvars[i]=="grazing") {
#     mdat <- subset(mdat, mowing!="applied")
#     mdat <- droplevels(mdat)
#   }
  
  
  modform <- as.formula(paste("success ~ ", mgmtvars[i], "*species + study.length + analysis2 + sample.size + (1|reference) + (1|country)", sep=""))
  m.ind.sp[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
}

warningmessages <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages

mdat <- subset(dat, species!="dunlin" & species!="ruff" & mowing!="applied")
mdat <- droplevels(mdat)

m.test <- glmer(success ~ AE.level*species + study.length + analysis2 + sample.size + (1|country), data=mdat, family=binomial)
summary(m.test)
slot(m.test, "optinfo")$conv$lme4$messages

tt <- getME(m.test,"theta")
ll <- getME(m.test,"lower")
min(tt[ll==0])

print(Sys.time())

setwd(outputwd)
sink(paste("abundance_0b.txt", sep=" "))
cat("\n########==========  0b) success of individual management types by species ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))
sink()


#------------------------------ 1a) overall success of AES and nature reserves -------------------------

m <- list()

print(Sys.time())

m[[1]] <- glmer(success ~ AE.level + reserve.desig + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())

#------------------------------ 1b) success of specific management interventions -------------------------

print(Sys.time())

m[[2]] <- glmer(success ~ mowing + grazing + fertpest + nest.protect + predator.control + water + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())

#------------------------------ 2a) success of AES and nature reserves by species -------------------------

# subset dataset to remove dunlin and ruff, not enough data for these species
mdat <- subset(dat, species!="dunlin" & species!="ruff")
mdat <- droplevels(mdat)


print(Sys.time())

m[[3]] <- glmer(success ~ AE.level*species + reserve.desig*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())


#------------------------------ 2b) success of specific management interventions by species -------------------------

# subset dataset to remove dunlin and ruff, not enough data for these species
mdat <- subset(dat, species!="dunlin" & species!="ruff")
mdat <- droplevels(mdat)

print(Sys.time())

m[[4]] <- glmer(success ~ mowing*species + grazing*species + fertpest*species + nest.protect*species + predator.control*species + water*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################




