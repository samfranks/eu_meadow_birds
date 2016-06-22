#################################################################
#
#     Step 2: EU meadow birds meta-analysis - models evaluating success
#
#################################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr")

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
dat <- subset(dat0, select=c("reference","lit.type","score","country","study.length","habitat","species","overall.metric","metric","sample.size","analysis2","success","failure","outcome",mgmtvars))

### Combine habitat categories ###

# combine habitat categories even further into 3 different categories
# arable: includes arable and mixed arable/pastoral
# pastoral: includes pastoral and mixed pastoral/unenclosed
# unenclosed
dat$newhabitat <- ifelse((dat$habitat %in% c("arable","mixed arable/pastoral")), "arable", ifelse((dat$habitat %in% c("pastoral","mixed pastoral/unenclosed")), "pastoral", "unenclosed"))

### Identify studies where specific intervention is evaluated ###

# identify studies where a specific intervention is evaluated (could be in combination with others, and also with higher-level interventions)
# if any specific measure is used, then dat$spec.int.used=1, otherwise if all specific measures are 'none', then spec.int.used=0
dat$spec.int.used <- with(dat, ifelse(mowing=="none" & grazing=="none" & fertpest=="none" & nest.protect=="none" & predator.control=="none" & water=="none", 0, 1))


### Identify studies where high level interventions are evaluated ###

# identify studies where a high level intervention is evaluated (could be in combination with others, and also with specific interventions)
# if any higher level measure is used, then dat$high.int.used=1, otherwise if all high level measures are 'none', then high.int.used=0
dat$high.int.used <- with(dat, ifelse(AE.level=="none" & reserve.desig=="none", 0, 1))


### Identify abundance change studies ###

dat$new.metric <- with(dat, ifelse(dat$overall.metric=="abundance change", "abundance change", metric))

#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}

#================================== Test effect of nuisance variables on success for the full dataset ===========================

m.nui1 <- glmer(success ~ study.length + sample.size + analysis2 + lit.type*score + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.nui1)

setwd(outputwd)
sink(paste("model output_nuisance variables.txt", sep=" "))

cat("\n########==========  Nuisance variables - lme4 models ==========########\n", sep="\n")
print(summary(m.nui1))

sink()

# significance is given by Wald Z tests (default for summary.glmer())
# only significant effect of a nuisance variable is literature type (when 'score' is not included as a variable)
# primary literature study is more likely to be unsuccessful than successful
# include lit.type in subsequent analyses
# controlling for interaction between lit.type and score removes significant effect of lit.type, suggesting that variance in the data explained by literature type could be accounted for by the quality of the analysis

# proportion of scores for each literature type
#             good   medium     poor
# grey    0.224299 0.485981 0.289720
# primary 0.405350 0.465021 0.129630

#=================================  ANALYSIS  1 ===============================

# Analysis evaluates whether one or several management interventions has a significantly positive outcome (i.e. is successful) compared to the baseline reference situation (which in most cases is no management)
# reference level for all interventions is 'none', so management (either applied or reduced, depending on the intervention) is evaluated against the reference

# 0a) success of individual management types overall
# 0b) success of individual management types by species
# 0c) success of individual management types by metric
# 0d) success of individual management types by habitat

### Nuisance explanatory variables and random effects
# 'nuisance explanatory variables: study duration (continuous), sample size (categorical: small, medium, large), multivariate/univariate (categorical)
# Random effects: study id - habitat will be included later on - have taken out country as a random effect since it causes problems in the analysis




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
usedat <- list() # data subset used to run a model

for (i in 1:length(mgmtvars)) {
  
  print(mgmtvars[i])
  mdat <- dat[dat[,mgmtvars[i]]!="none",]
  mdat <- subset(mdat, species!="ruff") # subset out ruff because there are too few studies
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  #   if (mgmtvars[i]=="mowing") {
  #     mdat <- subset(mdat, mowing!="applied")
  #   }
  #   if (mgmtvars[i]=="fertpest") {
  #     mdat <- subset(mdat, fertpest!="applied")
  #   }
  #   if (mgmtvars[i]=="predator.control") {
  #     mdat <- subset(mdat, predator.control!="reduced")
  #   }
  #   if (mgmtvars[i]=="water") {
  #     mdat <- subset(mdat, water!="reduced")
  #   }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$success))
  
  table(mdat[,mgmtvars[i]], mdat$species,mdat$success)
  
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], " + (1|reference) + (1|species)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ 1 + (1|reference) + (1|species)")
  }
  
  # run a normal glmer model
  m.ind[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  
}

names(m.ind) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

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

cat("\n########==========  0a) success of individual management types - lme4 models ==========########\n", sep="\n")
print(lapply(m.ind, summary))

cat("\n########==========  Warning messages lme4 models ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind, file=paste(workspacewd, "models_0a_lme4.rds", sep="/"))

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

sppdat <- dat

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
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
    
    # mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff", "dunlin/ruff", as.character(mdat$species)) # combine dunlin + ruff since only 1 obs of each
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
    mdat <- subset(mdat, species=="black-tailed godwit" | species=="lapwing") # model convergence issues when oyc and redshank are included
    # mdat <- subset(mdat, species!="snipe" & species!="curlew" & species!="dunlin" & species!="ruff")
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
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*species + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ species + (1|reference)")
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
  # mdat <- subset(mdat, species!="ruff" & species!="dunlin")
  
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
    mdat <- subset(mdat, fertpest!="applied")
    # mdat <- subset(mdat, metric!="occupancy")
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
    #     modform <- as.formula(paste("success ~ ", mgmtvars[i], "*metric + (1|reference) + (1|species)", sep=""))
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*metric + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    #     modform <- as.formula("success ~ metric + (1|reference) + (1|species)")
    modform <- as.formula("success ~ metric + (1|reference)")
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

cat("\n########==========  0b) success of individual management types by metric - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
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

# mgmtvars <- c("AE","AE.level","reserve.desig","nest.protect","water")


### ANALYSIS ###

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$newhabitat, dat[,mgmtvars[i]])
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
  table(mdat[,mgmtvars[i]], mdat$newhabitat, mdat$success)
  
  if (mgmtvars[i]=="AE.level") {
    mdat <- subset(mdat, newhabitat!="unenclosed")
  }
  
  if (mgmtvars[i]=="mowing") {
    mdat <- subset(mdat, newhabitat!="unenclosed")
  }
  
  if (mgmtvars[i]=="fertpest") {
    mdat <- subset(mdat, newhabitat!="unenclosed")
    mdat <- subset(mdat, fertpest!="applied")
  }
  
  if (mgmtvars[i]=="predator.control") {
    mdat <- subset(mdat, predator.control!="reduced")
  }
  
  if (mgmtvars[i]=="water") {
    mdat <- subset(mdat, water!="reduced")
  }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$newhabitat, mdat$success))
  
  
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*newhabitat + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("success ~ newhabitat + (1|reference)")
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
  summary(m.ind.sp.blme[[i]])
  
  
  
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

cat("\n########==========  0d) success of individual management types (subset only) by habitat - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
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





#=================================  ANALYSIS  2 ===============================


#-----------------   HIGH-LEVEL INTERVENTION SUCCESS  ------------------

mdat <- subset(dat, high.int.used==1)
mdat <- droplevels(mdat)

newlevels <- data.frame(unique(mdat[,c("AE.level","reserve.desig")]), AE.reserve=c("AE.basic-no reserve","AE.higher-no reserve","no AE-reserve","AE.basic-reserve","AE.higher-reserve"))

mdat <- merge(mdat, newlevels, by=c("AE.level","reserve.desig"))

mdat$AE.reserve <- relevel(mdat$AE.reserve, ref="no AE-reserve")
print(levels(mdat$AE.reserve))

m.high <- lme(success ~ AE.reserve + species, random = ~1|reference, data=mdat)

summary(m.high)

setwd(outputwd)
sink(paste("model output_2a.txt", sep=" "))
cat("\n########==========  Success of higher-level interventions combined ==========########\n", sep="\n")
print(summary(m.high))

sink()

### Save individual interventions models
saveRDS(m.high, file=paste(workspacewd, "models_2a.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2a.rds", sep="/"))



#-----------------   SPECIFIC INTERVENTION SUCCESS  ------------------

# when every specific intervention type = none, we won't evaluate studies which don't explicitly test the effect of a specific intervention (i.e. remove these records from being included in the model evaluating success of different interventions)
# For all other interventions, they may or may not be being used in combination with the focal intervention of the study, but often, this is not stated in the paper or is not tested explicitly in the study - ASSUMPTION is that these 'other interventions' which may be used are staying constant between the treatments of the focal intervention and the control (i.e. all else assumed equal between control + treatment)

# identify studies where a specific intervention is evaluated (could be in combination with others, and also with higher-level interventions)
# if any specific measure is used, then dat$spec.int.used=1, otherwise if all specific measures are 'none', then spec.int.used=0

# table(dat$spec.int.used)
# 0   1 
# 242 351

# subset data to only use records where specific interventions were tested
mdat <- subset(dat, spec.int.used==1)
mdat <- droplevels(mdat)
m.spec <- lme(success ~ mowing + grazing + fertpest + nest.protect + predator.control + water + species, random = ~1|reference, data=mdat)

summary(m.spec)

setwd(outputwd)
sink(paste("model output_2b.txt", sep=" "))

cat("\n########==========  Success of specific interventions combined ==========########\n", sep="\n")
print(summary(m.spec))

sink()

### Save individual interventions models
saveRDS(m.spec, file=paste(workspacewd, "models_2b.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2b.rds", sep="/"))



#-----------------   HIGH-LEVEL INTERVENTION SUCCESS by metric  ------------------

mdat <- subset(dat, new.metric!="survival" & new.metric!="recruitment")
mdat <- droplevels(mdat)

m.high.metric <- lme(success ~ AE.level*new.metric + reserve.desig*new.metric + species, random = ~1|reference, data=mdat)

summary(m.high.metric)

setwd(outputwd)
sink(paste("model output_2c.txt", sep=" "))
cat("\n########==========  Success of higher-level interventions combined, by metric ==========########\n", sep="\n")
print(summary(m.high.metric))

sink()

### Save individual interventions models
saveRDS(m.high.metric, file=paste(workspacewd, "models_2c.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2c.rds", sep="/"))



#-----------------   HIGH-LEVEL INTERVENTION SUCCESS by habitat  ------------------

out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$newhabitat, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# $AE
# 
# none applied
# arable       55      46
# pastoral    205     208
# unenclosed   71       8
# 
# $AE.level
# 
# none basic higher
# arable       55    11     35
# pastoral    205   132     76
# unenclosed   71     8      0
# 
# $reserve.desig
# 
# none applied
# arable      101       0
# pastoral    274     139
# unenclosed   59      20

# m.high.hab <- lme(success ~ AE*newhabitat + reserve.desig*newhabitat + species, random = ~1|reference, data=dat)

##### model doesn't run because of singularities (probably due to low-no unenclosed habitats using AES and no arable habitats in reserves) ####

#-----------------   SPECIFIC INTERVENTION SUCCESS by metric  ------------------

# remove categories with too few records, which make the model collapse
mdat <- subset(dat, new.metric!="survival" & new.metric!="recruitment" & new.metric!="occupancy" & new.metric!="abundance change")
mdat <- subset(mdat, grazing!="reduced" & mowing!="reduced" & fertpest!="applied" & predator.control!="reduced" & water!="reduced")

# subset data to only use records where specific interventions were tested
mdat <- subset(mdat, spec.int.used==1)
mdat <- droplevels(mdat)


# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(mdat$new.metric, mdat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# m.spec.metric <- lme(success ~ new.metric*mowing + new.metric*grazing + new.metric*fertpest + new.metric*nest.protect + new.metric*predator.control + new.metric*water + species, random = ~1|reference, data=mdat)

##### model doesn't run because of singularities  ####



#-----------------   SPECIFIC INTERVENTION SUCCESS by habitat  ------------------

##### model doesn't run because of singularities ####


#=================================  PLOT OUTPUTS - ANALYSIS  2 ===============================

setwd(outputwd)


### Read individual interventions models
m.high <- readRDS(file=paste(workspacewd, "models_2a.rds", sep="/"))
m.spec <- readRDS(file=paste(workspacewd, "models_2b.rds", sep="/"))
m.high.metric <- readRDS(file=paste(workspacewd, "models_2c.rds", sep="/"))

### Read dataset
dat.high <- readRDS(file=paste(workspacewd, "model dataset_2a.rds", sep="/"))
dat.spec <- readRDS(file=paste(workspacewd, "model dataset_2b.rds", sep="/"))
dat.high.metric <- readRDS(file=paste(workspacewd, "model dataset_2c.rds", sep="/"))


#-------------  HIGHER LEVEL INTERVENTIONS  -------------#

###----  Produce plotting dataset predictions ----###

plotmod <- m.high # model to plot results
origdat <- dat.high # original dataset

# unique.mgmtvars <- unique(origdat[,c("AE.level","reserve.desig")]) # unique combination of mgmtvars appearing in the original dataset
# newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

unique.mgmtvars <- unique(origdat$AE.reserve)

newdat <- data.frame(AE.reserve=rep(unique.mgmtvars, times=length(levels(origdat$species))), species=rep(levels(origdat$species), each=length(unique.mgmtvars)))

newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))

Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar)
newdat$lwr <- newdat$pred - (1.96*newdat$SE)
newdat$upr <- newdat$pred + (1.96*newdat$SE)
newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)

# fits <- newdat[,c("AE.level","reserve.desig","species","pred","SE","lwr","upr")]
# fits <- unique(fits)

fits <- newdat[,c("AE.reserve","species","pred","SE","lwr","upr")]
fits <- unique(fits)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.reserve=fits$AE.reserve), mean)

plotdat <- sum.fits

plotdat <- merge(plotdat, unique(origdat[,c("AE.level","reserve.desig","AE.reserve")]), by="AE.reserve")
plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE.level),]

###---- Output plot ----###

setwd(outputwd)
png("2a_high level combination intervention success.png", res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,6.5,2,3))

x <- c(1:nrow(plotdat))

plotdat$pch <- c(1,2,15,16,17)

plot(plotdat$pred~x, pch=plotdat$pch, cex=1.5, ylim=c(-0.2,1.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=4.5)

dev.off()

# 
# setwd(outputwd)
# png("2a_high level combination intervention success.png", res=300, height=12, width=15, units="in", pointsize=20)
# 
# par(mar=c(6,6,2,2))
# 
# x <- c(1:nrow(plotdat))
# 
# plotdat$pch <- c(0,1,2,15,16,17)
# 
# plot(plotdat$pred~x, pch=plotdat$pch, cex=1.5, ylim=c(0,1), xaxt="n", xlab="", ylab="", las=1, bty="n")
# arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
# arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
# abline(h=0.05, lty=3, lwd=2)
# # axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
# axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
# title(xlab="Intervention combination", cex.lab=1.5, font=2, line=5)
# 
# dev.off()

#-------------  SPECIFIC INTERVENTIONS  -------------#

###----  Produce plotting dataset predictions ----###

plotmod <- m.spec # model to plot results
origdat <- dat.spec # original dataset

unique.mgmtvars <- unique(origdat[,mgmtvars[4:9]]) # unique combination of mgmtvars appearing in the original dataset

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))

Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar)
newdat$lwr <- newdat$pred - (1.96*newdat$SE)
newdat$upr <- newdat$pred + (1.96*newdat$SE)
newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)

fits <- newdat
fits <- unique(fits)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(mowing=fits$mowing, grazing=fits$grazing, fertpest=fits$fertpest, nest.protect=fits$nest.protect, predator.control=fits$predator.control, water=fits$water), mean)

plotdat <- sum.fits
# plotdat <- plotdat[order(plotdat$mowing, plotdat$grazing, plotdat$fertpest, plotdat$nest.protect, plotdat$predator.control, plotdat$water),]

plotdat <- plotdat[order(plotdat$pred),]




###---- Output plot ----###

setwd(outputwd)
png("2b_specific combination intervention success.png", res=300, height=15, width=30, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(1,8,1,2))

x <- c(1:nrow(plotdat))-0.5

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(-0.3,1), xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(1,27))
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
# axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=0, xpd=TRUE)
title(ylab="Predicted probability of success \n (significant positive impact) ", cex.lab=1.5, font=2, line=3)

# dev.off()


# par(new=T)

y <- length(mgmtvars[4:9]):1

x <- c(1:nrow(plotdat))
new.x <- rep(x, each=max(y))
new.y <- rep(y, times=max(x))
tab <- data.frame(x=new.x,y=new.y)
tab <- tab[order(tab$y, decreasing=TRUE),]
labs <- plotdat[,1:6]

labs.long <- gather(labs, intervention, level, mowing:water)

# unicode up arrow: &#x2191;
# unicode down arrow: &#x2193;

tab.filled <- data.frame(tab, labs.long)
tab.filled[4] <- apply(tab.filled[4], 2, function(x) {
  gsub("none", "", x)
})
# tab.filled[4] <- apply(tab.filled[4], 2, function(x) {
#   gsub("applied", 2191, x)
# })
# tab.filled[4] <- apply(tab.filled[4], 2, function(x) {
#   gsub("reduced", 2193, x)
# })

par(mar=c(1,8,1,2))
plot(tab, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(1.2,6.8), xlim=c(1,27))
# points(x=tab$x-0.5, y=tab$y+0.5, pch=tab.filled$level, cex=1.5)
abline(v=c(0,x,29), h=c(y,7), lty=1)
axis(2, y+0.5, labels=c("mowing","grazing","fertiliser/\npesticides", "nest protection","predator control", "water"), las=1, cex.axis=1, font=2,tick=FALSE)
text(tab$x-0.5, tab$y+0.5, labels=ifelse(tab.filled$level=="applied", "\U2191", ifelse(tab.filled$level=="reduced", "\U2193", tab.filled$level)), cex=2)

dev.off()



#-------------  HIGHER LEVEL INTERVENTIONS - METRIC  -------------#

###----  Produce plotting dataset predictions ----###

plotmod <- m.high.metric # model to plot results
origdat <- dat.high.metric # original dataset

unique.mgmtvars <- unique(origdat[,c("AE.level","reserve.desig")]) # unique combination of mgmtvars appearing in the original dataset

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

newdat <- data.frame(newdat[rep(seq_len(nrow(newdat)), times=length(levels(as.factor(origdat$new.metric)))),], new.metric=rep(levels(as.factor(origdat$new.metric)), each=nrow(newdat)))

newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))

Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar)
newdat$lwr <- newdat$pred - (1.96*newdat$SE)
newdat$upr <- newdat$pred + (1.96*newdat$SE)
newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)

fits <- newdat[,c("AE.level","reserve.desig","species","new.metric","pred","SE","lwr","upr")]
fits <- unique(fits)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.level=fits$AE.level, reserve.desig=fits$reserve.desig, new.metric=fits$new.metric), mean)

plotdat <- sum.fits

# set up colours and plotting points for plot
cols <- rev(grey(seq(from=0,to=1,length.out = 6)))
colpch <- data.frame(cols, AE.level=rep(levels(plotdat$AE.level), times=2), reserve.desig=rep(levels(plotdat$reserve.desig), each=3))
pchs <- data.frame(pch=c(21,22,23,24), new.metric=levels(plotdat$new.metric))

plotdat <- merge(plotdat, colpch, by=c("AE.level","reserve.desig"))
plotdat <- merge(plotdat, pchs, by=c("new.metric"))

plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE.level, plotdat$new.metric),]


###---- Output plot ----###

setwd(outputwd)
png("2c_high level combination intervention success_metric.png", res=300, height=10, width=25, units="in", pointsize=20)

par(mar=c(7,6,2,2))

x <- c(1:nrow(plotdat))

plot(plotdat$pred~x, pch=plotdat$pch, bg=as.character(plotdat$cols), cex=1.5, ylim=c(-0.2,1.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
abline(v=12.5, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
axis(1, seq(from=2.5, by=4, length.out=6), labels=rep("", 6), tick=TRUE)
text(seq(from=2.5, by=4, length.out=6), par("usr")[3]*1.15, srt = 0, pos=1, xpd = TRUE, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), cex=1)
text(c(4,16), par("usr")[3]*2, srt = 0, pos=4, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1.2)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=5.5)

legend(min(x)-0.5,1.2, legend=unique(plotdat$new.metric), pch=plotdat$pch, pt.cex=1.2, bty="n", xpd=TRUE)

dev.off()




