#############################################################################
#
#     Step 5: EU meadow birds analysis - COMBINATION INTERVENTION MODELS
#
#############################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016


# =================================  SET LOGIC STATEMENTS  ====================


# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme","dplyr")

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
easyPredCI <- function(model,newdata,alpha=0.05) {
  
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



# =================================  ANALYSIS  3 - INTERVENTIONS IN COMBINATION ===============================


#-----------------   3a) POLICY INTERVENTION SUCCESS  ------------------


#### ---- Method 1 (USED IN SUBMISSION): Run models using only those records which evaluate either AES and/or nature reserves ---- ####


# this tells us the probability of success IF THESE INTERVENTIONS WERE ATTEMPTED (subtly different than below Method 2, which tells us the probability of success if they were used vs not used)
# re-structured the data to combine AE and reserves into a single variable because otherwise lack of the none/none category produces rank deficiency in the model and it drops a coefficient

# > table(mdat$AE,mdat$reserve.desig)
# 
#         none applied
# none       0      85
# applied  190      72

mdat <- subset(dat, high.int.used==1)
mdat <- subset(mdat, species!="ruff" & species!="dunlin")
mdat <- droplevels(mdat)

newlevels <- data.frame(AE=c("applied","none","applied"), reserve.desig=c("none","applied","applied"), AE.reserve=c("AE-no reserve","no AE-reserve","AE-reserve"))

mdat <- merge(mdat, newlevels, by=c("AE","reserve.desig"))

mdat$AE.reserve <- relevel(mdat$AE.reserve, ref="no AE-reserve")
table(mdat$AE.reserve)
# no AE-reserve AE-no reserve    AE-reserve 
#            85           190            72

newlevels2 <- data.frame(AE.level=c("basic","higher","none","basic","higher"), reserve.desig=c("none","none","applied","applied","applied"), AE.reserve2=c("AE.basic-no reserve","AE.higher-no reserve","no AE-reserve","AE.basic-reserve","AE.higher-reserve"))

mdat <- merge(mdat, newlevels2, by=c("AE.level","reserve.desig"))



# --------  MODEL ---------

m.high <- glmer(success ~ AE.reserve + species + (1|reference), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.high)


# -------------- Output model results ----------------------

library(multcomp) # contains the ghlt function

setwd(outputwd)
sink(paste("model output_analysis 3a_method 1.txt", sep=" "))

cat("\n########==========  3a) Success of policy interventions combined ==========########\n", sep="\n")
print(summary(m.high))

cat("\n###----  Tukey contrasts  ---###\n", sep="\n")
print(summary(glht(m.high, linfct = mcp(AE.reserve="Tukey"))))

sink()

### Save individual interventions models
saveRDS(m.high, file=paste(workspacewd, "models_analysis 3a_method 1.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_analysis 3a_method 1.rds", sep="/"))





#-----------------   SPECIFIC INTERVENTION SUCCESS  ------------------

# when every specific intervention type = none, we won't evaluate studies which don't explicitly test the effect of a specific intervention (i.e. remove these records from being included in the model evaluating success of different interventions)
# this means that we are not evaluating the success of used vs not used, but rather if used, which interventions are most successful (while controlling for the use of other interventions), and also what combination of interventions are most successful?

# need summary outputs showing the overall effect of each intervention (while controlling for the others)

# For all other interventions, they may or may not be being used in combination with the focal intervention of the study, but often, this is not stated in the paper or is not tested explicitly in the study - ASSUMPTION is that these 'other interventions' which may be used are staying constant between the treatments of the focal intervention and the control (i.e. all else assumed equal between control + treatment)

# identify studies where a specific intervention is evaluated (could be in combination with others, and also with higher-level interventions)
# if any specific measure is used, then dat$spec.int.used=1, otherwise if all specific measures are 'none', then spec.int.used=0

# table(dat$spec.int.used)
# 0   1 
# 242 351

# subset data to only use records where specific interventions were tested
mdat <- subset(dat, spec.int.used==1)
mdat <- subset(mdat, species!="ruff")
mdat <- droplevels(mdat)

# add a variable totalling the number of interventions tested at once
interventions.only <- mdat[,c("mowing","grazing","fertpest","nest.protect","predator.control","water")]
num.interventions <- ifelse(interventions.only=="none",0,1)
num.interventions.sum <- apply(num.interventions, 1, sum)

mdat <- data.frame(mdat, num.interventions.sum)

# model which includes number of interventions tested at once is rank deficient, so omit this variable
m.spec <- glmer(success ~ mowing + grazing + fertpest + water + nest.protect + predator.control + species + (1|reference), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.spec)

# ------ LRT single-term deletion of combination management interventions from global model --------------

vars <- c("mowing","grazing","fertpest","water","nest.protect","predator.control")

m <- list()
for (i in 1:length(vars)) {
  m[[i]] <- drop1(m.spec, scope = as.formula(paste("~", vars[i], sep="")), test="Chisq")
}

setwd(outputwd)
sink(paste("model output_analysis 3b.txt", sep=" "))

cat("\n########==========  3b) Success of specific interventions combined ==========########\n", sep="\n")
print(summary(m.spec))


cat("\n###---  Likelihood Ratio Tests ---###\n", sep="\n")
print(m)

cat("\n###----  Significance of each intervention (in the presence of others)  ---###\n", sep="\n")
print(drop1(m.spec, scope = ~mowing, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~grazing, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~fertpest, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~water, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~nest.protect, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~predator.control, test="Chisq"))


sink()

### Save individual interventions models
saveRDS(m.spec, file=paste(workspacewd, "models_analysis 3b.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_analysis 3b.rds", sep="/"))


# ---------- OUTPUT TABLE of combination management interventions LRT results (dplyr/piping way) -------

lapply(m, function(x) {data.frame(df=x[,"Df"][2], LRT=x[,"LRT"][2], pval=x[,"Pr(Chi)"][2]) %>% return
}) %>%
  do.call(rbind, .) %>%
  data.frame(vars) %>%
  write.csv(paste(outputwd, "model output_combination intervention LRT table.csv", sep="/"), row.names=FALSE)



# --------  Produce plotting dataset predictions for interventions applied in combination ---------

# Read model
plotmod <- readRDS(file=paste(workspacewd, "models_analysis 3b.rds", sep="/"))
# Read model dataset
origdat <- readRDS(file=paste(workspacewd, "model dataset_analysis 3b.rds", sep="/")) # original dataset

unique.mgmtvars <- unique(origdat[,mgmtvars[4:9]]) # unique combination of mgmtvars appearing in the original dataset

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

pred <- predict(plotmod, newdat, type="response", re.form=NA)
pred.CI <- easyPredCI(plotmod, newdat)
fits <- data.frame(newdat, pred, pred.CI)

# produce mean population level prediction for all intervention combinations across species
sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(mowing=fits$mowing, grazing=fits$grazing, fertpest=fits$fertpest, water=fits$water, nest.protect=fits$nest.protect, predator.control=fits$predator.control), mean)
plotdat <- sum.fits
plotdat <- plotdat[order(plotdat$mowing, plotdat$grazing, plotdat$fertpest, plotdat$water, plotdat$nest.protect, plotdat$predator.control),]

# ---- Success of individual combinations of interventions -----

# count up number of interventions applied at once and add sum to plotdat
interventions.only <- plotdat[,1:6]
num.interventions <- ifelse(interventions.only=="none",0,1)
num.interventions.sum <- apply(num.interventions, 1, sum)

plotdat <- data.frame(plotdat, num.interventions.sum)

# add variable showing whether success is 'significant' (i.e. lwr does not overlap 0.05)
plotdat$sig <- ifelse(plotdat$lwr > 0.05, "Y", "N")

# order by 'significance' and magnitude of success
plotdat <- plotdat[order(plotdat$sig, plotdat$pred),]

# Save dataframe file for predicted success of all combinations of management interventions applied
saveRDS(plotdat, file=paste(workspacewd, "analysis 3b_all combinations plotting dataset.rds", sep="/"))

# ---- Averaged success of interventions, accounting for application in combiantion -----

# produce overall effect of intervention level, averaged across predicted values of other intervention levels, and also averaged across all species
avg.fits <- list()
for (i in 1:6) {
  temp <- aggregate(sum.fits[,c("pred","lwr","upr")], by=list(level=sum.fits[,i]), mean)
  temp <- data.frame(temp, intervention=names(sum.fits[i]))
  avg.fits[[i]] <- temp[-which(temp[,1]=="none"),]
}
plotdat <- do.call(rbind, avg.fits)

saveRDS(plotdat, file=paste(workspacewd, "analysis 3b_average intervention plotting dataset.rds", sep="/"))
write.csv(plotdat, "analysis 3b_average intervention probabilities and CIs.csv", row.names=FALSE)


