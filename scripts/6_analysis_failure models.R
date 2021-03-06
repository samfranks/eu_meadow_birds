#################################################################
#
#     Step 6: EU meadow birds analysis - FAILURE MODELS
#
#################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016


# =================================  SET LOGIC STATEMENTS  ====================


# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)


# =================================  LOAD FUNCTIONS =================================

### Ben Bolker's function for calculating CIs on predictions from a merMod object and plotting the results from his RPubs GLMM worked examples
# http://rpubs.com/bbolker/glmmchapter
# by specifying re.form=NA we're saying that we want the population-level prediction, i.e. setting the random effects to zero and getting a prediction for an average (or unknown) group
# Computing confidence intervals on the predicted values is relatively easy if we're willing to completely ignore the random effects, and the uncertainty of the random effects
# this easy method produces similar width CIs to using the bootMer function in lme4, perhaps slightly wider CIs in some cases

# can change to alpha=0.16, approximately equal to 84% CIs
easyPredCI <- function(model,newdata,alpha=alphalevel) {
  
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
                    newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link (logistic) function: could also use plogis()
  linkinv <- model@resp$family$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(lwr=pred0-crit*pred.se,
                upr=pred0+crit*pred.se))
  
}

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



#=================================  ANALYSIS  1b - failures of management interventions ===============================


# identify which categories have low numbers
out <- list()
for(i in 1:length(mgmtvars)) {
  out[[i]] <- table(dat$failure, dat[,mgmtvars[i]])
}
names(out) <- mgmtvars
out

# set up list to output models and model datasets to
m.ind <- list()
usedat <- list() # data subset used to run a model

# remove predator control from list of interventions because it causes unsolvable convergence issues and huge SEs
mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","water","nest.protect","predator.control")
mgmtvars <- mgmtvars[-which(mgmtvars %in% "predator.control")]  # even if 'reduced' predator control is converted to applied, model outputs are still way off with really large variance on the random effects and huge error on the parameter estimate. Continue to leave this intervention out.

# ---------------- Models ------------------

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
  #     if (mgmtvars[i]=="predator.control") {
  #       mdat <- subset(mdat, predator.control!="reduced")
  #     }
  #   if (mgmtvars[i]=="water") {
  #     mdat <- subset(mdat, water!="reduced")
  #   }
  
  mdat <- droplevels(mdat)
  usedat[[i]] <- mdat
  
  (checkzeros <- table(mdat[,mgmtvars[i]], mdat$failure))
  
  # create different formulas to use depending on whether management variable is 1 or 2 levels
  if (length(levels(mdat[,mgmtvars[i]])) > 1) {
    modform <- as.formula(paste("failure ~ ", mgmtvars[i], " + (1|reference) + (1|species)", sep=""))
  }
  
  if (length(levels(mdat[,mgmtvars[i]])) < 2) {
    modform <- as.formula("failure ~ 1 + (1|reference) + (1|species)")
  }
  
  # run a normal glmer model
  m.ind[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
  # run a bglmer model  
  #   vcov.dim <- nrow(vcov(m.ind[[i]]))
  #   m.ind.blme[[i]] <- bglmer(modform, data=mdat, family=binomial, fixef.prior = normal(cov = diag(9,vcov.dim)), control=glmerControl(optimizer="bobyqa"))
  
  
}

names(m.ind) <- mgmtvars
# names(m.ind.blme) <- mgmtvars
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

setwd(outputwd)
sink(paste("model output_analysis 1b_failure.txt", sep=" "))

cat("\n########==========  Analysis 1b) failure of individual management types - lme4 models ==========########\n", sep="\n")
print(lapply(m.ind, summary))

cat("\n########==========  Warning messages lme4 models ==========########\n", sep="\n")
print(warningmessages.lme4)
sink()

### Save individual interventions models
saveRDS(m.ind, file=paste(workspacewd, "models_analysis 1b_lme4.rds", sep="/"))

### Save dataset for 0a models
saveRDS(usedat, file=paste(workspacewd, "model dataset_analysis 1b.rds", sep="/"))



#=================================  OUTPUT PARAMETER TABLES  ===============================

setwd(outputwd)

moddat <- readRDS(paste(workspacewd, "model dataset_analysis 1b.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_analysis 1b_lme4.rds", sep="/"))
alphalevel <- 0.05
successlevel <- 0.05


#------------------------------ 1b) failure of individual management types -------------------------


# -------    Output model parameter table   -----------

# Output model coefficient tables for each management type, and convert parameter table to a dataframe instead of a matrix
coeftab <- lapply(mod, function(x) summary(x)$coefficients)
coeftab <- lapply(coeftab, function(x) {
  out <- as.data.frame(x)
  out$parameter <- rownames(out)
  rownames(out) <- 1:nrow(out)
  return(out) })
coeftab2 <- do.call(rbind, coeftab)
coeftab3 <- data.frame(coeftab2, mgmtvar=rep(names(coeftab),lapply(coeftab,nrow)))
coeftab3$mgmtvar <- as.character(coeftab3$mgmtvar)
rownames(coeftab3) <- c(1:nrow(coeftab3))
partable <- coeftab3
partable <- partable[,c(6,5,1,2,4)] # omit z value column
names(partable) <- c("Management intervention","Parameter level","Estimate","SE","p-value")

n <- unlist(lapply(moddat, nrow))
n <- data.frame(management=names(n), n)
names(n) <- c("Management intervention","n")

partable <- merge(partable, n, sort=FALSE)

# Write the parameter table
write.csv(format(partable, scientific=FALSE, digits=2),  "parameter table_analysis 1b_failure_lme4.csv", row.names=FALSE)




