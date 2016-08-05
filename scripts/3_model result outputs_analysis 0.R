#######################################################################
#
#     Step 3: EU meadow birds meta-analysis - output Step 2 model results
#
#######################################################################

# Samantha Franks
# 18 March 2016

set.seed(2)

# =================================  SET LOGIC STATEMENTS  ====================

# default to plot when all are FALSE is results from overall analysis (0a)
species <- FALSE # plot the species-specific model results (0b)
metric <- FALSE # plot the metric-specific model results (0c)
habitat <- TRUE # plot the habitat-specific model results (0d)

bias <- FALSE

alphalevel <- 0.05
successlevel <- 0.05

# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme")

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
outputwd <- paste(parentwd, "output", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


# =================================  LOAD DATA & MODELS  ===============================


mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")
if (bias) mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","water")

# load model data & models: use blme models for all 0a-d analyses

if (!bias & !species & !habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0a.rds", sep="/"))
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

# =================================  OUTPUT PARAMETER TABLES  ===============================

setwd(outputwd)

# ------------------------------ 0a) success of individual management types -------------------------

if (!bias & !species & !metric & !habitat) {
  
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
  
  partable <- merge(partable, n)
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "0a_overall parameter table.csv", row.names=FALSE)
  
}

if (bias & !species & !habitat & !metric) {
  
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
  
  partable <- merge(partable, n)
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "0a_overall parameter table_bias.csv", row.names=FALSE)
  
  
}

# ------------------------------ 0b) success of individual management types by species -------------------------


if (species & !metric & !habitat) {
  
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
  
  partable <- merge(partable, n)
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "0b_species-specific parameter table.csv", row.names=FALSE)
  
}

# ------------------------------ 0c) success of individual management types by metric -------------------------


if (!species & metric & !habitat) {
  
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
  
  partable <- merge(partable, n)
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "0c_metric-specific parameter table.csv", row.names=FALSE)
  
}

# ------------------------------ 0d) success of individual management types by habitat -------------------------


if (!species & !metric & habitat) {
  
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
  
  partable <- merge(partable, n)
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "0d_habitat-specific parameter table.csv", row.names=FALSE)
  
}


# =================================  PLOT MODEL OUTPUTS  ===============================

# set the wd to output to
setwd(outputwd)

# ------------------------------ 0a) success of individual management types -------------------------

if (!bias & !species & !metric & !habitat) {
  
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
  
  names(plotdat) <- mgmtvars
  
  plotfinal.policy <- do.call(rbind, plotdat[c("AE","reserve.desig")])
  plotfinal.AE <- do.call(rbind, plotdat[c("AE.level")])
  plotfinal.mgmt <- do.call(rbind, plotdat[4:9])
  
  ### -------- Output plot --------###
  
  ### ---- Policy level interventions plot ----
  png("0aA_overall model results.png", res=300, height=12, width=20, units="in", pointsize=20)
  
  par(mfrow=c(1,2))
  
  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.policy
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
  # axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
  text(x, par("usr")[3]-0.06, srt = 0, pos=1, xpd = TRUE, labels=c("AES","nature reserve"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  # dev.off()
  
  
  ### ---- AES level interventions plot ----
  # png("0aB_overall model results.png", res=300, height=12, width=12, units="in", pointsize=20)
  
  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.AE
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
  # axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
  text(x, par("usr")[3]-0.06, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  dev.off()
  
  ### ---- Specific management interventions plot ----
  
  png("0aB_overall model results.png", res=300, height=12, width=30, units="in", pointsize=20)
  
  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.mgmt
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
  axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
  text(x, par("usr")[3]-0.06, srt = 0, pos=1, xpd = TRUE, labels=c("mowing applied", "mowing reduced", "grazing applied", "grazing reduced", "fertiliser/pesticides \n applied","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water \n applied", "water \n reduced"))
  # text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE, labels=c("AES","basic-level \n AES","higher-level \n AES","nature reserve/ \n designation", "mowing applied", "mowing reduced", "grazing applied", "grazing reduced", "fertiliser/pesticides \n applied","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water \n applied", "water \n reduced"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  dev.off()
  
  ### -------- Output table of predicted probabilities +/- CIs --------###
  write.csv(do.call(rbind,plotdat), "0a_overall probabilities and CIs.csv", row.names=FALSE)
  
}


if (bias & !species & !metric & !habitat) {
  
  plotdat <- list()
  n <- list()
  
  for (i in 1:length(mod)) {
    
    # dataset to predict over is the same as the original dataset
    pred <- predict(mod[[i]], type="response", re.form=NA)
    pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
    
    fits <- data.frame(pred, pred.CI, biased.metric=moddat[[i]]$biased.metric, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]), mgmt.type=i)
    unique.fits <- unique(fits)
    
    plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar, mgmt.type=unique.fits$mgmt.type, biased.metric=unique.fits$biased.metric), mean)
    
    
  }
  
  plotfinal <- do.call(rbind, plotdat)
  
  n <- length(levels(plotfinal$biased.metric))
  
  pch <- data.frame(biased.metric=levels(plotfinal$biased.metric), pch=c(21,22), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  plotfinal <- merge(plotfinal,pch, by="biased.metric")
  
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$biased.metric),]
  
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  
  ### -------- Output plot --------###
  
  if (alphalevel==0.05) {
    png("0a_overall model results_bias.png", res=300, height=12, width=30, units="in", pointsize=20)
  } else {png("0a_overall model results_84CIs_bias.png", res=300, height=12, width=30, units="in", pointsize=20)}
  
  par(mar=c(7,6,3,2))
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \napplied", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","more water \napplied", "water \nreduced"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Management intervention evaluated", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  
  legend(min(plotfinal$rowid)-1,1.1, legend=c("change/productivity metrics","abundance/occupancy metrics"), pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  ### -------- Output table of predicted probabilities +/- CIs --------###
  if (alphalevel==0.05) {
    write.csv(plotfinal[,c("biased.metric","pred","lwr","upr","mgmtvar")], "0a_overall probabilities and CIs_bias.csv", row.names=FALSE)
  } else {write.csv(plotfinal[,c("biased.metric","pred","lwr","upr","mgmtvar")], "0a_overall probabilities and 84CIs_bias.csv", row.names=FALSE)}
  
  
}

# ------------------------------ 0b) success of individual management types by species -------------------------

if (species & !metric & !habitat) {
  
  plotdat <- list()
  
  for (i in 1:length(mod)) {
    
    # dataset to predict over is the same as the original dataset
    pred <- predict(mod[[i]], type="response", re.form=NA)
    pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
    
    fits <- data.frame(pred, pred.CI, species=moddat[[i]]$species, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]), mgmt.type=i)
    unique.fits <- unique(fits)
    
    plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar, mgmt.type=unique.fits$mgmt.type, species=unique.fits$species), mean)
    
  }
  
  names(plotdat) <- mgmtvars
  
  plotfinal.policy <- do.call(rbind, plotdat[c("AE","reserve.desig")])
  plotfinal.AE <- do.call(rbind, plotdat[c("AE.level")])
  plotfinal.mgmt <- do.call(rbind, plotdat[4:9])
  
  ### -------- Output plot --------###
  
  ### ---- Policy level interventions plot ----
  
  png("0bA_species-specific model results.png", res=300, height=20, width=20, units="in", pointsize=20)
  
  par(mfrow=c(2,1))
  
  par(mar=c(3,6,4,2))
  
  plotfinal <- plotfinal.policy
  maxspecies <- levels(do.call(rbind, plotdat)$species)
  n <- length(maxspecies)
  
  set.seed(2)
  pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), replace=TRUE, n))
  pch
  plotfinal <- merge(plotfinal,pch, by="species")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","nature reserve"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  # title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.2, legend=levels(plotfinal$species), pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  # dev.off()
  
  ### ---- AES level interventions plot ----
  
  # png("0bB_species-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)

  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.AE
  maxspecies <- levels(do.call(rbind, plotdat)$species)
  n <- length(maxspecies)
  
  set.seed(2)
  pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), replace=TRUE, n))
  pch
  plotfinal <- merge(plotfinal,pch, by="species")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  # text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  
  # legend(min(plotfinal$rowid)-1,1.2, legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  ### ---- Specific management interventions plot ----
  png("0bB_species-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)
  
  par(mar=c(7,6,4,2))
  
  plotfinal <- plotfinal.mgmt
  maxspecies <- levels(do.call(rbind, plotdat)$species)
  n <- length(maxspecies)
  
  set.seed(2)
  pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), replace=TRUE, n))
  pch
  plotfinal <- merge(plotfinal,pch, by="species")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  # text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing reduced", "grazing applied", "grazing reduced","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water applied"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.2, legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  ### -------- Output table of predicted probabilities +/- CIs for all 0b --------###
  write.csv(do.call(rbind,plotdat), "0b_species-specific probabilities and CIs.csv", row.names=FALSE)

}

# ------------------------------ 0c) success of individual management types by metric -------------------------

if (!species & metric & !habitat) {
  
  plotdat <- list()
  
  for (i in 1:length(mod)) {
    
    # dataset to predict over is the same as the original dataset
    pred <- predict(mod[[i]], type="response", re.form=NA)
    pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
    
    fits <- data.frame(pred, pred.CI, metric=moddat[[i]]$new.metric, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]), mgmt.type=i)
    unique.fits <- unique(fits)
    
    plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar, mgmt.type=unique.fits$mgmt.type, metric=unique.fits$metric), mean)
    
  }
  
  names(plotdat) <- mgmtvars
  
  plotfinal.policy <- do.call(rbind, plotdat[c("AE","reserve.desig")])
  plotfinal.AE <- do.call(rbind, plotdat[c("AE.level")])
  plotfinal.mgmt <- do.call(rbind, plotdat[4:9])
  
  ### -------- Output plot --------###
  
  ### ---- Policy level interventions plot ----
  
  png("0cA_metric-specific model results.png", res=300, height=20, width=18, units="in", pointsize=20)
  
  par(mfrow=c(2,1))
  
  par(mar=c(3,6,4,2))
  
  plotfinal <- plotfinal.policy
  maxmetric <- levels(do.call(rbind, plotdat)$metric)
  n <- length(maxmetric)
  
  set.seed(2)
  pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="metric")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid)+0.2), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","nature reserve"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  # title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.1, legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  # dev.off()
  
  
  ### ---- AES level interventions plot ----
  

  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.AE
  maxmetric <- levels(do.call(rbind, plotdat)$metric)
  n <- length(maxmetric)
  
  set.seed(2)
  pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="metric")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid)+0.2), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  # text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  

  dev.off()
  
  
  
  ### ---- Specific management interventions plot ----
  
  png("0cB_metric-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)
  
  par(mar=c(7,6,4,2))
  
  plotfinal <- plotfinal.mgmt
  maxmetric <- levels(do.call(rbind, plotdat)$metric)
  n <- length(maxmetric)
  
  set.seed(2)
  pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="metric")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing applied","mowing reduced", "grazing applied", "grazing reduced","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water applied", "water reduced"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.2, legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  ### -------- Output table of predicted probabilities +/- CIs --------###
  write.csv(do.call(rbind, plotdat), "0c_metric-specific probabilities and CIs.csv", row.names=FALSE)

  
}

# ------------------------------ 0d) success of individual management types by habitat -------------------------

if (!species & !metric & habitat) {
  
  plotdat <- list()
  
  for (i in 1:length(mod)) {
    
    # dataset to predict over is the same as the original dataset
    pred <- predict(mod[[i]], type="response", re.form=NA)
    pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
    
    fits <- data.frame(pred, pred.CI, habitat=moddat[[i]]$newhabitat, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]), mgmt.type=i)
    unique.fits <- unique(fits)
    
    plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar, mgmt.type=unique.fits$mgmt.type, habitat=unique.fits$habitat), mean)
    
    
  }
  
  names(plotdat) <- mgmtvars
  
  plotfinal.policy <- do.call(rbind, plotdat[c("AE","reserve.desig")])
  plotfinal.AE <- do.call(rbind, plotdat[c("AE.level")])
  plotfinal.mgmt <- do.call(rbind, plotdat[4:9])
  
  ### -------- Output plot --------###
  
  ### ---- Policy level interventions plot ----
  
  png("0dA_habitat-specific model results.png", res=300, height=20, width=16, units="in", pointsize=20)
  
  par(mfrow=c(2,1))
  
  par(mar=c(3,6,4,2))
  
  plotfinal <- plotfinal.policy
  maxhabitat <- levels(do.call(rbind, plotdat)$habitat)
  n <- length(maxhabitat)
  
  set.seed(2)
  pch <- data.frame(habitat=maxhabitat, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="habitat")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$habitat),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid))+0.2, pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","nature reserve"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  # title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=successlevel, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.1, legend=pch$habitat, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  # dev.off()
  
  
  ### ---- AES level interventions plot ----
  
  
  par(mar=c(7,6,2,2))
  
  plotfinal <- plotfinal.AE
  maxhabitat <- levels(do.call(rbind, plotdat)$habitat)
  n <- length(maxhabitat)
  
  set.seed(2)
  pch <- data.frame(habitat=maxhabitat, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="habitat")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$habitat),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid))+0.2, pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  # axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  # text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  
  
  dev.off()
  
  
  
  ### ---- Specific management interventions plot ----
  
  png("0dB_habitat-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)
  
  par(mar=c(7,6,4,2))
  
  plotfinal <- plotfinal.mgmt
  maxhabitat <- levels(do.call(rbind, plotdat)$habitat)
  n <- length(maxhabitat)
  
  set.seed(2)
  pch <- data.frame(habitat=maxhabitat, pch=rep(c(21,22,23),length.out=n), col=sample(grey(seq(from=0,to=1,length.out = n)), n))
  pch
  plotfinal <- merge(plotfinal,pch, by="habitat")
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$habitat),]
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  x <- c(min(plotfinal$rowid):max(plotfinal$rowid))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing applied","mowing reduced", "grazing applied", "grazing reduced","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water applied"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.05, lty=3, lwd=2)
  
  legend(min(plotfinal$rowid),1.2, legend=pch$habitat, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  
  ### -------- Output table of predicted probabilities +/- CIs --------###
  write.csv(do.call(rbind, plotdat), "0d_habitat-specific probabilities and CIs.csv", row.names=FALSE)

  
}


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################