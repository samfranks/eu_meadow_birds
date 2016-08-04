#######################################################################
#
#     Step 5: EU meadow birds meta-analysis - output results for J Appl Ecol submission
#
#######################################################################

# Samantha Franks
# 4 Aug 2016

set.seed(2)

#=================================  SET LOGIC STATEMENTS  ====================

# default to plot when all are FALSE is results from overall analysis (0a)
species <- FALSE # plot the species-specific model results (0b)
metric <- FALSE # plot the metric-specific model results (0c)
habitat <- FALSE # plot the habitat-specific model results (0d)


alphalevel <- 0.05
successlevel <- 0.05

#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)


#=================================  LOAD FUNCTIONS =================================

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
outputwd <- paste(parentwd, "output/submission", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)

#==============================  SET OUTPUT DIRECTORY  ===================================

setwd(outputwd)

#================================    MANAGEMENT VARIABLE NAMES    ===============================


mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")



#==============================  FIGURE 1  ===================================

# Fig 1a = AES, site protection (Analysis 0a)
# Fig 1b = basic vs higher AES (Analysis 0a)
# Fig 1c = combined effects of AES and site protection (Analysis 2a)


moddat <- readRDS(paste(workspacewd, "model dataset_0a.rds", sep="/")) #[c("AE","reserve.desig")]
mod <- readRDS(paste(workspacewd, "models_0a_lme4.rds", sep="/")) #[c("AE","reserve.desig")]

plotdat <- list()
n <- list()

for (i in 1:length(mod)) {
  
  # dataset to predict over is the same as the original dataset
  pred <- predict(mod[[i]], type="response", re.form=NA)
  pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
  n[i] <- nrow(moddat[[i]])
  
  fits <- data.frame(pred,pred.CI,lit.type=moddat[[i]][,"lit.type"],mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]))
  unique.fits <- unique(fits)
  
  plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar), mean)
  
}

names(plotdat) <- mgmtvars #c("AE","reserve.desig")

fig1a <- do.call(rbind, plotdat[c("AE","reserve.desig")])
fig1b <- do.call(rbind, plotdat[c("AE.level")])


###-------- Output plot --------###


###---- Fig 1a: Policy level interventions plot ----

png("Fig1a_AES_site protection.png", res=300, height=12, width=20, units="in", pointsize=24)

par(mfrow=c(1,2))

par(mar=c(7,6,2,2))

plotfinal <- fig1a

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2)

text(0.6, 1, "a)")



###---- Fig 1b: AES level interventions plot ----

par(mar=c(7,3,2,2))

plotfinal <- fig1b

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES","higher-level AES"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1, line=3)
abline(h=successlevel, lty=3, lwd=2)

text(0.6, 1, "b)")

mtext("Intervention", side=1, outer=TRUE, line=-2, cex=1.2)

dev.off()







# load model data & models: use blme models for all 0a-d analyses

if (species & !habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0b.rds", sep="/"))
if (!species & !habitat & metric) moddat <- readRDS(paste(workspacewd, "model dataset_0c.rds", sep="/"))
if (!species & habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0d.rds", sep="/"))
if (bias & !species & !habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0a_bias.rds", sep="/"))

if (bias) moddat <- moddat[-8] # remove predator control model (huge variance on the random effects, and SEs on fixed effects are very large)

# load models
if (!bias & !species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0a_lme4.rds", sep="/"))
if (species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0b_blme.rds", sep="/"))
if (!species & !habitat & metric) mod <- readRDS(paste(workspacewd, "models_0c_blme.rds", sep="/"))
if (!species & habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0d_blme.rds", sep="/"))
if (bias & !species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0a_lme4_bias.rds", sep="/"))

if (bias) mod <- mod[-8] # remove predator control model (huge variance)

#=================================  OUTPUT PARAMETER TABLES  ===============================

setwd(outputwd)