###################################################################################################################
#
#     Step 4: EU meadow birds analysis - INDIVIDUAL INTERVENTION MODELS (overall, by metric, by species)
#
###################################################################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016


# =================================  SET LOGIC STATEMENTS  ====================

bias <- FALSE # setting this logic control variable to true means introducing the confounding covariate fixed effect 'biased.metric' into the intervention models
# we have decided not to do this, so this control variable is set permanently to 'FALSE'

# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)


# =================================  LOAD FUNCTIONS =================================



# =================================  SET DIRECTORY STRUCTURE  ================================

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
outputwd <- paste(parentwd, "output/revision Dec 2016", sep="/")
workspacewd <- paste(parentwd, "workspaces/revision Dec 2016", sep="/")

options(digits=6)


# =================================  LOAD DATA  ===============================

source(paste(scriptswd, "source_model data preparation.R", sep="/"))





# =================================  ANALYSIS  1 - INDIVIDUAL INTERVENTIONS OVERALL ===============================


### NOTES
# Analysis evaluates the overall success rate (probability that an intervention has a significant positive outcome) of each intervention individually

### Confounding covariates
# Literature type was the only significant confounding covariate in the nuisance variable analysis, but inclusion as a fixed effect causes issues with model convergence
# instead, include study ID as a random effect to account for some of this variation


# ---------------- Models ------------------

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
  
  if (!bias) {
    # create different formulas to use depending on whether management variable is 1 or 2 levels
    if (length(levels(mdat[,mgmtvars[i]])) > 1) {
      modform <- as.formula(paste("success ~ ", mgmtvars[i], " + (1|reference) + (1|species)", sep=""))
    }
    
    if (length(levels(mdat[,mgmtvars[i]])) < 2) {
      modform <- as.formula("success ~ 1 + (1|reference) + (1|species)")
    }
  }
  
  # if (bias) {
  #   # create different formulas to use depending on whether management variable is 1 or 2 levels
  #   if (length(levels(mdat[,mgmtvars[i]])) > 1) {
  #     modform <- as.formula(paste("success ~ ", mgmtvars[i], " + biased.metric + (1|reference) + (1|species)", sep=""))
  #   }
  #   
  #   if (length(levels(mdat[,mgmtvars[i]])) < 2) {
  #     modform <- as.formula("success ~ 1 + biased.metric + (1|reference) + (1|species)")
  #   }
  # }
  # 
  # # run a normal glmer model
  # m.ind[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  
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



# -------------- Output model results ----------------------

if (!bias) {
  setwd(outputwd)
  sink(paste("model output_analysis 1.txt", sep=" "))
  
  cat("\n########==========  Analysis 1) success of individual management types - lme4 models ==========########\n", sep="\n")
  print(lapply(m.ind, summary))
  
  cat("\n########==========  Warning messages lme4 models ==========########\n", sep="\n")
  print(warningmessages.lme4)
  sink()
  
  ### Save individual interventions models
  saveRDS(m.ind, file=paste(workspacewd, "models_analysis 1_lme4.rds", sep="/"))
  
  ### Save dataset for 0a models
  saveRDS(usedat, file=paste(workspacewd, "model dataset_analysis 1.rds", sep="/"))
}


# if (bias) {
#   setwd(outputwd)
#   sink(paste("model output_0a_bias.txt", sep=" "))
#   
#   cat("\n########==========  0a) success of individual management types - lme4 models ==========########\n", sep="\n")
#   print(lapply(m.ind, summary))
#   
#   cat("\n########==========  Warning messages lme4 models ==========########\n", sep="\n")
#   print(warningmessages.lme4)
#   sink()
#   
#   ### Save individual interventions models
#   saveRDS(m.ind, file=paste(workspacewd, "models_0a_lme4_bias.rds", sep="/"))
#   
#   ### Save dataset for 0a models
#   saveRDS(usedat, file=paste(workspacewd, "model dataset_0a_bias.rds", sep="/"))
# }



# ===========================  ANALYSIS  2 - INTERVENTIONS BY SPECIES and METRIC-SPECIFIC RESPONSE  ===========================

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




# ------------------- MODELS:  2a) success of individual management types by species  ----------------

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
    # mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff" | mdat$species=="curlew", "dunlin/ruff/curlew", as.character(mdat$species)) # combine dunlin + ruff + curlew
    mdat <- subset(mdat, species!="dunlin" & species!="ruff" & species!="curlew")
    mdat <- subset(mdat, mowing!="applied")
  }
  
  if (mgmtvars[i]=="grazing") {
    # mdat$species <- ifelse(mdat$species=="dunlin" | mdat$species=="ruff", "dunlin/ruff", as.character(mdat$species)) # combine dunlin + ruff since only 1 obs of each
    mdat <- subset(mdat, species!="dunlin" & species!="ruff" & species!="curlew")
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
    # mdat <- subset(mdat, predator.control!="reduced")
    mdat <- subset(mdat, species!="snipe" & species!="oystercatcher" & species!="redshank")
  }
  
  if (mgmtvars[i]=="water") {
    mdat <- subset(mdat, species!="dunlin" & species!="ruff")
    mdat <- subset(mdat, water!="reduced")
    # model runs ok without 2 of curlew, oystercatcher and snipe, but model won't converge and has a very high max|grad| value when more than 1 of these species is included. Since they all have no successes for this management intervention and similar levels of failure, then combine together for the water analysis
    mdat <- subset(mdat, species!="curlew" & species!="snipe")
    # mdat$species <- ifelse(mdat$species=="curlew" | mdat$species=="snipe", "curlew/snipe", as.character(mdat$species))
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
  
  # if ANY checkzeros==0, then there are categories which are missing  observations, so will have problems with complete separation and/or convergence
  # use bglmer since there are some cases of singularity produced by 0/1 not having any observations for some of the categorical variables
  # calculate the dimensions of the covariance matrix for bglmer, based on the dimensions of the covariance matrix from the regular glmer model
  
  vcov.dim <- nrow(vcov(m.ind.sp[[i]]))
  m.ind.sp.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))

  
}


names(m.ind.sp) <- mgmtvars
names(m.ind.sp.blme) <- mgmtvars
names(usedat) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.sp.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme



# -------------- Output model results ----------------------

setwd(outputwd)
sink(paste("model output_analysis 2a.txt", sep=" "))

cat("\n########==========  Analysis 2a) success of individual management types by species - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.sp.blme, summary))

cat("\n########==========  Analysis 2a) success of individual management types by species - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind.sp.blme, file=paste(workspacewd, "models_analysis 2a_blme.rds", sep="/"))
saveRDS(m.ind.sp, file=paste(workspacewd, "models_analysis 2a_lme4.rds", sep="/"))

### Save dataset for 0b models
saveRDS(usedat, file=paste(workspacewd, "model dataset_analysis 2a.rds", sep="/"))






# ------------------- MODELS:  2b) success of individual management types by metric  ----------------


# subset dataset to remove recruitment and survival metrics, too few observations
metricdat <- subset(dat, new.metric!="recruitment" & new.metric!="survival")
metricdat <- droplevels(metricdat)

# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(metricdat$new.metric, metricdat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind.sp <- list()
m.ind.sp.blme <- list()
usedat <- list() # data subset used to run a model
checkzeros <- list()

for (i in 1:length(mgmtvars)) {
  
  mgmtvars[i]
  mdat <- metricdat[metricdat[,mgmtvars[i]]!="none",]
  table(mdat[,mgmtvars[i]], mdat$new.metric, mdat$success)
  # mdat <- subset(mdat, species!="ruff" & species!="dunlin")
  
  # for the following categories, subset further because there aren't enough observations of either 0,1 or both
  #   if (mgmtvars[i]=="reserve.desig") {
  #     mdat <- subset(mdat, new.metric!="productivity")
  #   }
  
  #   if (mgmtvars[i]=="mowing") {
  #     mdat <- subset(mdat, mowing!="applied")
  #   }
  
  if (mgmtvars[i]=="grazing") {
    mdat <- subset(mdat, new.metric!="abundance/occupancy change")
  }
  
  if (mgmtvars[i]=="fertpest") {
    mdat <- subset(mdat, fertpest!="applied")
  }

  if (mgmtvars[i]=="water") {
    mdat <- subset(mdat, new.metric!="abundance/occupancy change")
  }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros[[i]] <- table(mdat[,mgmtvars[i]], mdat$new.metric, mdat$success))
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    #     modform <- as.formula(paste("success ~ ", mgmtvars[i], "*new.metric + (1|reference) + (1|species)", sep=""))
    modform <- as.formula(paste("success ~ ", mgmtvars[i], "*new.metric + (1|reference)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    #     modform <- as.formula("success ~ new.metric + (1|reference) + (1|species)")
    modform <- as.formula("success ~ new.metric + (1|reference)")
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
names(checkzeros) <- mgmtvars

warningmessages.lme4 <- lapply(m.ind.sp, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.lme4

warningmessages.blme <- lapply(m.ind.sp.blme, function(x) slot(x, "optinfo")$conv$lme4$messages)
warningmessages.blme


# -------------- Output model results ----------------------


setwd(outputwd)
sink(paste("model output_analysis 2b.txt", sep=" "))

cat("\n########==========  Analysis 2b) success of individual management types by metric - BLME models (good) ==========########\n", sep="\n")
print(lapply(m.ind.sp.blme, summary))

cat("\n########==========  Analysis 2b) success of individual management types by metric - lme4 models (convergence issues) ==========########\n", sep="\n")
print(lapply(m.ind.sp, summary))

cat("\n########==========  Warning messages BLME models (good) ==========########\n", sep="\n")
print(warningmessages.blme)

cat("\n########==========  Warning messages lme4 models (convergence issues) ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind.sp.blme, file=paste(workspacewd, "models_analysis 2b_blme.rds", sep="/"))
saveRDS(m.ind.sp, file=paste(workspacewd, "models_analysis 2b_lme4.rds", sep="/"))

### Save dataset for 0b models
saveRDS(usedat, file=paste(workspacewd, "model dataset_analysis 2b.rds", sep="/"))



##########################         END OF ANALYSIS FOR J Applied Ecology submission/revision     ##################




