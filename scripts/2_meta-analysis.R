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

dat0 <- readRDS(paste(workspacewd, "meadow birds analysis dataset_full.rds", sep="/"))

mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# subset dataset for analysis to desired columns only
dat <- subset(dat0, select=c("reference","lit.type","country","study.length","habitat","species","overall.metric","metric","sample.size","analysis2","success",mgmtvars))



#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}

#=================================  ANALYSIS  ===============================

# Analysis evaluates whether one or several management interventions has a significantly positive outcome (i.e. is successful) compared to the baseline reference situation (which in most cases is no management)
# reference level for all interventions is 'none', so management (either applied or reduced, depending on the intervention) is evaluated against the reference

# 0a) success of individual management types overall
# 0b) success of individual management types by species
# 0c) success of individual management types by metric
# 0d) success of individual management types by habitat

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
# 'nuisance explanatory variables: study duration (continuous), sample size (categorical: small, medium, large), multivariate/univariate (categorical)
# Random effects: study id - habitat will be included later on - have taken out country as a random effect since it causes problems in the analysis


#------------------------------ Test effect of nuisance variables on success for the full dataset -------------------------

m.nui1 <- glmer(success ~ study.length + sample.size + analysis2 + lit.type + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.nui1)

# only significant effect of a nuisance variable is literature type
# primary literature study is more likely to be unsuccessful than successful
# include lit.type in subsequent analyses

#------------------------------ 0a) success of individual management types -------------------------

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$success, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind <- list()
m.ind.blme <- list()
usedat <- list() # data subset used to run a model

for (i in 1:length(mgmtvars)) {
  
  print(mgmtvars[i])
  mdat <- dat[dat[,mgmtvars[i]]!="none",]
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  #   if (mgmtvars[i]=="mowing") {
  #     mdat <- subset(mdat, mowing!="applied")
  #   }
  #   if (mgmtvars[i]=="fertpest") {
  #     mdat <- subset(mdat, fertpest!="applied")
  #   }
  if (mgmtvars[i]=="predator.control") {
    mdat <- subset(mdat, predator.control!="reduced")
  }
  #   if (mgmtvars[i]=="water") {
  #     mdat <- subset(mdat, water!="reduced")
  #   }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$success))
  
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], " + lit.type + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ 1 + lit.type + (1|reference)")
  }
  
  # run a normal glmer model
  m.ind[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  # vars with warnings using glmer: mowing, fertpest, predator.control
  
  # use bglmer since there are some cases of singularity produced by 0/1 not having any observations for some of the categorical variables
  # calculate the dimensions of the covariance matrix for bglmer (diagonal matrix number of fixed effects) based on the number of levels present of all factor variables in the model - 1; if lit.type is included as a variable, then add the number of levels of lit.type to the calculations for the dimensions of the covariance matrix
  # vcov.dim <- length(levels(mdat[,mgmtvars[i]])) + length(levels(mdat$lit.type)) - 1
  
  vcov.dim <- nrow(vcov(m.ind[[i]]))
  m.ind.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))
  
  
}

names(m.ind) <- mgmtvars
names(m.ind.blme) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme

# ### lme4 convergence troubleshooting
# # check singularity
# tt <- getME(m.ind.blme[[5]],"theta")
# ll <- getME(m.ind.blme[[5]],"lower")
# min(tt[ll==0])
# 
# # check gradient calculations
# derivs1 <- m.ind.blme[[5]]@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

### Output model results ###

setwd(outputwd)
sink(paste("model output_0a.txt", sep=" "))

cat("\n########==========  0a) success of individual management types - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.blme, summary))

cat("\n########==========  0a) success of individual management types - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind, file=paste(workspacewd, "models_0a_lme4.rds", sep="/"))
saveRDS(m.ind.blme, file=paste(workspacewd, "models_0a_blme.rds", sep="/"))

### Save dataset for 0a models
saveRDS(usedat, file=paste(workspacewd, "model dataset_0a.rds", sep="/"))



#------------------------------ 0b) success of individual management types by species -------------------------

### ANALYTICAL ISSUES ###

# as a result of a dataset with lots of categories and sometimes few values in those categories, there have been multiple issues with convergence, huge parameter estimats and standard errors caused by complete separation or singularity
# often these problems can be solved by using the blme package, which uses a Bayesian framework superimposed over lme4 methods to add a weak prior to the fixed-effect parameters
# sometimes further reducing the complexity of the models is needed if even bglmer is throwing convergence warnings

# see http://stats.stackexchange.com/questions/132677/binomial-glmm-with-a-categorical-variable-with-full-successes and
# http://rpubs.com/bbolker/glmmchapter for how to solve complete separation issues
# problems with all 0's or all 1's in a single category level - complete separation
# try Ben Bolker's fix using the blme package
# this phenomenon is called complete separation. You can find quite a lot (now that you know its name) Googling around ... It is fairly thoroughly discussed here in a general context, and here in the context of GLMMs. The standard solution to this problem is to add a small term that pushes the parameters back toward zero -- in frequentist contexts this is called a penalized or bias-corrected method. The standard algorithm is due to Firth (1993, "Bias reduction of maximum likelihood estimates" Biometrika 80, 27-38), and is implemented in the logistf package on CRAN. In Bayesian contexts this is framed as adding a weak prior to the fixed-effect parameters.
# To my knowledge Firth's algorithm hasn't been extended to GLMMs, but you can use the Bayesian trick by using the blme package, which puts a thin Bayesian layer over the top of the lme4 package. Here's an example from the above-linked GLMM discussion:
# cmod_blme_L2 <- bglmer(predation~ttt+(1|block),data=newdat,
# family=binomial,
# fixef.prior = normal(cov = diag(9,4)))
# The first two lines in this example are exactly the same as we would use in the standard glmer model; the last specifies that the prior for the fixed effects is a multivariate normal distribution with a diagonal variance-covariance matrix. The matrix is 4x4 (because we have 4 fixed-effect parameters in this example), and the prior variance of each parameter is 9 (corresponding to a standard deviation of 3, which is pretty weak -- that means +/- 2SD is (-6,6), which is a very large range on the logit scale).
# The very large standard errors of the parameters in your example are an example of a phenomenon closely related to complete separation (it occurs whenever we get extreme parameter values in a logistic model) called the Hauck-Donner effect.



### ANALYSIS ###

# subset dataset to remove dunlin and ruff, not enough data for these species for any analysis that includes species as a covariate
# sppdat <- subset(dat, species!="dunlin" & species!="ruff")
sppdat <- dat
sppdat <- droplevels(sppdat)

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(sppdat$species, sppdat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind.sp <- list()
m.ind.sp.blme <- list()
usedat <- list() # data subset used to run a model

for (i in 1:length(mgmtvars)) {
  
  mgmtvars[i]
  mdat <- sppdat[sppdat[,mgmtvars[i]]!="none",]
  table(mdat[,mgmtvars[i]], mdat$species, mdat$success)
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  # dunlin and ruff will need to be removed from most categories due to lack of observations, but there may be some with sufficient, or ok if combined
  if (mgmtvars[i]=="AE") {
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
  }
  
  if (mgmtvars[i]=="AE.level") {
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
  }
  
  if (mgmtvars[i]=="reserve.desig") {
    mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff", "dunlin/ruff", as.character(mdat$species)) # combine dunlin + ruff since only 1 obs of each
  }
  
  if (mgmtvars[i]=="mowing") {
    mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff" | mdat$species=="curlew", "dunlin/ruff/curlew", as.character(mdat$species)) # combine dunlin + ruff + curlew
    # mdat <- subset(mdat, species!="dunlin" & species!="ruff" & species!="curlew")
    mdat <- subset(mdat, mowing!="applied")
  }
  
  if (mgmtvars[i]=="grazing") {
    # mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff", "dunlin/ruff", as.character(mdat$species)) # combine dunlin + ruff since only 1 obs of each
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
    # mdat <- subset(mdat, grazing!="reduced")
    # mdat <- subset(mdat, species!="snipe" & species!="curlew")
  }
  
  if (mgmtvars[i]=="fertpest") {
    mdat <- subset(mdat, fertpest!="applied") # remove applied level as causes non-convergence even for bglmer model
    mdat <- subset(mdat, species!="snipe" & species!="curlew" & species!="dunlin" & species!="ruff")
  }
  
  if (mgmtvars[i]=="nest.protect") {
    mdat <- subset(mdat, species!="curlew" & species!="ruff") # 1 obs each; exclude rather than combine as interpretation for curlew/ruff combined doesn't make sense (too dissimilar in ecology)
  }
  
  if (mgmtvars[i]=="predator.control") {
    mdat <- subset(mdat, predator.control!="reduced")
    # mdat <- subset(mdat, species!="dunlin" & species!="ruff")
  }
  
  if (mgmtvars[i]=="water") {
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
    
    mdat <- subset(mdat, water!="reduced")
    # model runs ok without 2 of curlew, oystercatcher and snipe, but model won't converge and has a very high max|grad| value when more than 1 of these species is included. Since they all have no successes for this management intervention and similar levels of failure, then combine together for the water analysis
    # mdat <- subset(mdat, species!="curlew" & species!="oystercatcher")
    mdat$species <- ifelse(mdat$species=="curlew" | mdat$species=="snipe", "curlew/snipe", as.character(mdat$species))
  }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$species, mdat$success))
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*species + lit.type + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ species + lit.type + (1|reference)")
  }
  
  # run a normal glmer model
  m.ind.sp[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  # if ANY checkzeros==0, then there are categories which are missing any observations whatsoever, so will have problems with complete separation and/or convergence
  # use bglmer since there are some cases of singularity produced by 0/1 not having any observations for some of the categorical variables
  # calculate the dimensions of the covariance matrix for bglmer, based on the dimensions of the covariance matrix from the regular glmer model
  
  # if (any(checkzeros==0)) {
  vcov.dim <- nrow(vcov(m.ind.sp[[i]]))
  m.ind.sp.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))
  # }
  
}

names(m.ind.sp) <- mgmtvars
names(m.ind.sp.blme) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.sp.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme


setwd(outputwd)
sink(paste("model output_0b.txt", sep=" "))
cat("\n########==========  0b) success of individual management types by species - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.sp.blme, summary))

cat("\n########==========  0b) success of individual management types by species - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind.sp.blme, file=paste(workspacewd, "models_0b_blme.rds", sep="/"))
saveRDS(m.ind.sp, file=paste(workspacewd, "models_0b_lme4.rds", sep="/"))

### Save dataset for 0b models
saveRDS(usedat, file=paste(workspacewd, "model dataset_0b.rds", sep="/"))

#------------------------------ 0c) success of individual management types by metric -------------------------

### ANALYSIS ###

# subset dataset to remove recruitment and survival metrics, too few observations
metricdat <- subset(dat, metric!="recruitment" & metric!="survival")
metricdat <- droplevels(metricdat)

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(metricdat$metric, metricdat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind.sp <- list()
m.ind.sp.blme <- list()
usedat <- list() # data subset used to run a model

for (i in 1:length(mgmtvars)) {
  
  mgmtvars[i]
  mdat <- metricdat[metricdat[,mgmtvars[i]]!="none",]
  table(mdat[,mgmtvars[i]], mdat$metric, mdat$success)
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  #   if (mgmtvars[i]=="reserve.desig") {
  #     mdat <- subset(mdat, metric!="productivity")
  #   }
  
  #   if (mgmtvars[i]=="mowing") {
  #     mdat <- subset(mdat, mowing!="applied")
  #   }
  
  #   if (mgmtvars[i]=="grazing") {
  #     mdat <- subset(mdat, metric!="occupancy")
  #   }
  
  if (mgmtvars[i]=="fertpest") {
    # mdat <- subset(mdat, fertpest!="applied")
    mdat <- subset(mdat, metric!="occupancy")
  }
  
  if (mgmtvars[i]=="predator.control") {
    mdat <- subset(mdat, predator.control!="reduced")
  }
  
  if (mgmtvars[i]=="water") {
    # mdat <- subset(mdat, water!="reduced")
    mdat <- subset(mdat, metric!="occupancy")
    # model runs ok without 2 of curlew, oystercatcher and snipe, but model won't converge and has a very high max|grad| value when more than 1 of these species is included. Since they all have no successes for this management intervention and similar levels of failure, then combine together for the water analysis
    # mdat <- subset(mdat, species!="curlew" & species!="oystercatcher")
    # mdat$species <- ifelse(mdat$species=="oystercatcher" | mdat$species=="curlew" | mdat$species=="snipe", "OC/CU/SN", as.character(mdat$species))
  }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$metric, mdat$success))
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*metric + lit.type + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ metric + lit.type + (1|reference)")
  }
  
  # run a normal glmer model
  m.ind.sp[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  # if ANY checkzeros==0, then there are categories which are missing any observations whatsoever, so will have problems with complete separation and/or convergence
  # use bglmer since there are some cases of singularity produced by 0/1 not having any observations for some of the categorical variables
  # calculate the dimensions of the covariance matrix for bglmer, based on the dimensions of the covariance matrix from the regular glmer model
  
  # if (any(checkzeros==0)) {
  vcov.dim <- nrow(vcov(m.ind.sp[[i]]))
  m.ind.sp.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))
  # }
  
}

names(m.ind.sp) <- mgmtvars
names(m.ind.sp.blme) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.sp.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme


setwd(outputwd)
sink(paste("model output_0c.txt", sep=" "))
cat("\n########==========  0b) success of individual management types by metric - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.sp.blme, summary))

cat("\n########==========  0b) success of individual management types by metric - lme4 models (also good) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (also good) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind.sp.blme, file=paste(workspacewd, "models_0c_blme.rds", sep="/"))
saveRDS(m.ind.sp, file=paste(workspacewd, "models_0c_lme4.rds", sep="/"))

### Save dataset for 0b models
saveRDS(usedat, file=paste(workspacewd, "model dataset_0c.rds", sep="/"))

#######################################################
#######################################################
#######################################################        End of analysis for 18 March report delivery
#######################################################
#######################################################


#------------------------------ 0d) success of individual management types by habitat -------------------------

### NOTES on Analysis 0d ###

# the effect of habitat type on intervention success won't be very applicable to certain habitat types (e.g. no mowing done in arable habitats, pastoral only), so look at only a select few management types used across habitats to gauge intervention effectiveness in different habitats
# management types: AE, AE.level, reserve.desig, nest.protect, water
# reserve.desig model run without literature type because not enough observations in the different literature types to produce convergence (model with lit.type was rank deficient)

mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")


# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$habitat, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

mgmtvars <- c("AE","AE.level","reserve.desig","nest.protect","water")

### ANALYSIS ###

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$habitat, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind.sp <- list()
m.ind.sp.blme <- list()
usedat <- list() # data subset used to run a model

for (i in 1:length(mgmtvars)) {
  
  
  mgmtvars[i]
  mdat <- dat[dat[,mgmtvars[i]]!="none",]
  table(mdat[,mgmtvars[i]], mdat$habitat, mdat$success)
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  if (mgmtvars[i]=="AE") {
    mdat <- subset(mdat, habitat!="unenclosed")
  }
  if (mgmtvars[i]=="AE.level") {
    mdat <- subset(mdat, habitat!="unenclosed")
  }
  
#   # for the following categories, subset further because there aren't enough observations of either 0,1 or both
#   if (mgmtvars[i]=="reserve.desig") {
#     mdat <- subset(mdat, !(grepl("arable", mdat$habitat)) & habitat!="mixed pastoral/unenclosed")
#     mdat <- subset(mdat, habitat=="pastoral")
#   }
  

  if (mgmtvars[i]=="water") {
    mdat <- subset(mdat, water!="reduced")
    # mdat <- subset(mdat, metric!="occupancy")
    # model runs ok without 2 of curlew, oystercatcher and snipe, but model won't converge and has a very high max|grad| value when more than 1 of these species is included. Since they all have no successes for this management intervention and similar levels of failure, then combine together for the water analysis
    # mdat <- subset(mdat, species!="curlew" & species!="oystercatcher")
    # mdat$species <- ifelse(mdat$species=="oystercatcher" | mdat$species=="curlew" | mdat$species=="snipe", "OC/CU/SN", as.character(mdat$species))
  }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$habitat, mdat$success))
  

  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*habitat + lit.type + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ habitat + lit.type + (1|reference)")
  }
  
  # if model produces rank deficiency because of missing observations from certain habitat types for certain literature types, use below formula - only the case for nature reserves management type
  if (mgmtvars[i] %in% c("reserve.desig")) {
    modform <- as.formula("success ~ habitat + (1|reference)")
  }
  
  # run a normal glmer model
  m.ind.sp[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  # if ANY checkzeros==0, then there are categories which are missing any observations whatsoever, so will have problems with complete separation and/or convergence
  # use bglmer since there are some cases of singularity produced by 0/1 not having any observations for some of the categorical variables
  # calculate the dimensions of the covariance matrix for bglmer, based on the dimensions of the covariance matrix from the regular glmer model
  
  # if (any(checkzeros==0)) {
  vcov.dim <- nrow(vcov(m.ind.sp[[i]]))
  m.ind.sp.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))
  # }
  

  
}

names(m.ind.sp) <- mgmtvars
names(m.ind.sp.blme) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.sp.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme


setwd(outputwd)
sink(paste("model output_0d.txt", sep=" "))
cat("\n########==========  0d) success of individual management types (subset only) by habitat - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.sp.blme, summary))

cat("\n########==========  0d) success of individual management types (subset only) by habitat - lme4 models (also good) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (also good) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind.sp.blme, file=paste(workspacewd, "models_0d_blme.rds", sep="/"))
saveRDS(m.ind.sp, file=paste(workspacewd, "models_0d_lme4.rds", sep="/"))

### Save dataset for 0b models
saveRDS(usedat, file=paste(workspacewd, "model dataset_0d.rds", sep="/"))

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################

# 
# #------------------------------ 1a) overall success of AES and nature reserves -------------------------
# 
# m1 <- list()
# 
# print(Sys.time())
# 
# m1[[1]] <- glmer(success ~ AE.level + reserve.desig + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# print(Sys.time())
# 
# #------------------------------ 1b) success of specific management interventions -------------------------
# 
# # subset out management categories with few observations
# mdat <- subset(dat, mowing!="applied" & fertpest!="applied" & predator.control!="reduced")
# mdat <- droplevels(mdat)
# 
# print(Sys.time())
# 
# m1[[2]] <- glmer(success ~ mowing + grazing + fertpest + nest.protect + predator.control + water + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# print(Sys.time())
# 
# warningmessages <- lapply(m1, function(x) slot(x, "optinfo")$conv$lme4$messages)
# warningmessages
# 
# ### Output model results ###
# 
# setwd(outputwd)
# sink(paste("model output_1ab.txt", sep=" "))
# cat("\n########==========  1) success of AES/natures reserves and specific management interventions (species pooled) ==========########\n", sep="\n")
# print(lapply(m1, summary))
# print(warningmessages)
# sink()
# 
# ### Save individual interventions models
# saveRDS(m1, file=paste(workspacewd, "models_1ab.rds", sep="/"))
# 
# 
# 
# #------------------------------ 2a) success of AES and nature reserves by species -------------------------
# 
# # subset dataset to remove dunlin and ruff, not enough data for these species
# mdat <- subset(dat, species!="dunlin" & species!="ruff")
# mdat <- droplevels(mdat)
# 
# 
# print(Sys.time())
# 
# m[[3]] <- glmer(success ~ AE.level*species + reserve.desig*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# print(Sys.time())
# 
# #######################################################
# #######################################################
# #######################################################
# #######################################################
# #######################################################
# 
# 
# #------------------------------ 2b) success of specific management interventions by species -------------------------
# 
# # subset dataset to remove dunlin and ruff, not enough data for these species
# mdat <- subset(dat, species!="dunlin" & species!="ruff")
# mdat <- droplevels(mdat)
# 
# print(Sys.time())
# 
# m[[4]] <- glmer(success ~ mowing*species + grazing*species + fertpest*species + nest.protect*species + predator.control*species + water*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# 
# print(Sys.time())
# 
# 
# #######################################################
# #######################################################
# #######################################################
# #######################################################
# #######################################################
# 
# 
# 
# 
