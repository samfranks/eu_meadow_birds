#######################################################################
#
#     EU meadow birds meta-analysis - OUTPUTTING MODEL RESULTS (graphs and tables)
#
#######################################################################

# Samantha Franks
# 18 March 2016

set.seed(2)

#=================================  SET LOGIC STATEMENTS  ====================

species <- TRUE # plot the species-specific model results (0b)
metric <- FALSE # plot the metric-specific model results (0c)
alphalevel <- 0.05

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
outputwd <- paste(parentwd, "output", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


#=================================  LOAD DATA & MODELS  ===============================

# load master dataset
dat0 <- readRDS(paste(workspacewd, "meadow birds analysis dataset.rds", sep="/"))
mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# subset master dataset to desired columns only
dat <- subset(dat0, select=c("reference","country","study.length","habitat","species","overall.metric","metric","sample.size","analysis2","success",mgmtvars))

# load model data & models
if (!metric) {
  if (!species) moddat <- readRDS(paste(workspacewd, "model dataset_0a.rds", sep="/"))
  if (species) moddat <- readRDS(paste(workspacewd, "model dataset_0b.rds", sep="/"))
  
  # load models
  if (!species) mod <- readRDS(paste(workspacewd, "models_0a.rds", sep="/"))
  if (species) mod <- readRDS(paste(workspacewd, "models_0b_blme.rds", sep="/"))
}

if (metric) {
  
  # load models & model data
  moddat <- readRDS(paste(workspacewd, "model dataset_0c.rds", sep="/"))
  mod <- readRDS(paste(workspacewd, "models_0c_blme.rds", sep="/"))
}

#=================================  LOAD DATA & MODELS  ===============================

if (!metric) {
  if (!species) {
    
    # Output model coefficient tables for each management type, and convert parameter table to a dataframe instead of a matrix
    coeftab <- lapply(mod, function(x) summary(x)$coefficients)
    coeftab2 <- do.call(rbind, coeftab)
    coeftab3 <- coeftab2
    rownames(coeftab3) <- c(1:nrow(coeftab3))
    coeftab3 <- as.data.frame(coeftab3)
    
    # create a variable with the names of the different management interventions and their levels (if present)
    mgmtvarlevels <- list()
    for (i in 1:length(mgmtvars)) {
      mgmtvarlevels[[i]] <- paste(mgmtvars[i], levels(moddat[[i]][,mgmtvars[i]]))
    }
    mgmtvarlevels <- unlist(mgmtvarlevels)
    
    # calculate sample sizes of the datasets
    n <- lapply(moddat, nrow)
    n <- unlist(n)
    
    # bind the two together to create a named parameter table with coefficients and their SEs against the relevant intervention
    partable <- cbind(mgmtvarlevels,coeftab3)
    partable <- partable[,c(1,2,3,5)] # omit z value column
    names(partable) <- c("Management intervention evaluated","Estimate","SE","p-value")
    
    # Write the parameter table
    write.csv(format(partable, scientific=FALSE, digits=2),  "overall parameter table.csv", row.names=FALSE)
    
  }
  
  if (species) {
    
    # Output model coefficient tables for each management type, and convert parameter table to a dataframe instead of a matrix
    coeftab <- lapply(mod, function(x) summary(x)$coefficients)
    coeftab2 <- do.call(rbind, coeftab)
    coeftab3 <- coeftab2
    rownames(coeftab3) <- c(1:nrow(coeftab3))
    coeftab3 <- as.data.frame(coeftab3)
    
    # bind the two together to create a named parameter table with coefficients and their SEs against the relevant intervention
    partable <- cbind(parameters=rownames(coeftab2),coeftab3)
    partable <- partable[,c(1,2,3,5)] # omit z value column
    names(partable) <- c("Management intervention evaluated","Estimate","SE","p-value")
    
    # Write the parameter table
    write.csv(format(partable, scientific=FALSE, digits=2),  "species-specific parameter table.csv", row.names=FALSE)
    
  }
}

if (metric) {
  
  # Output model coefficient tables for each management type, and convert parameter table to a dataframe instead of a matrix
  coeftab <- lapply(mod, function(x) summary(x)$coefficients)
  coeftab2 <- do.call(rbind, coeftab)
  coeftab3 <- coeftab2
  rownames(coeftab3) <- c(1:nrow(coeftab3))
  coeftab3 <- as.data.frame(coeftab3)
  
  # bind the two together to create a named parameter table with coefficients and their SEs against the relevant intervention
  partable <- cbind(parameters=rownames(coeftab2),coeftab3)
  partable <- partable[,c(1,2,3,5)] # omit z value column
  names(partable) <- c("Management intervention evaluated","Estimate","SE","p-value")
  
  # Write the parameter table
  write.csv(format(partable, scientific=FALSE, digits=2),  "metric-specific parameter table.csv", row.names=FALSE)
  
}


#=================================  PLOT MODEL OUTPUTS  ===============================

# set the wd to output to
setwd(outputwd)

if (!metric) {
  
  if (!species) {
    
    
    plotdat <- list()
    
    for (i in 1:length(mod)) {
      
      # dataset to predict over is the same as the original dataset
      pred <- predict(mod[[i]], type="response", re.form=NA)
      pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
      
      plotdat[[i]] <- unique(data.frame(pred, pred.CI, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]])))
      
    }
    
    plotfinal <- do.call(rbind, plotdat)
    
    ###-------- Output plot --------###
    png("overall model results.png", res=300, height=12, width=20, units="in", pointsize=20)
    par(mar=c(7,6,2,2))
    
    x <- c(1:nrow(plotfinal))
    
    plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n")
    axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
    text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE, labels=c("AES","basic-level \n AES","higher-level \n AES","nature reserve/ \n designation", "mowing reduced", "grazing applied", "grazing reduced", "fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","more water \n applied"))
    arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
    arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
    title(xlab="Management intervention evaluated", cex.lab=1.5, font=2, line=5)
    title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
    abline(h=0.5, lty=3, lwd=2)
    
    dev.off()
    
    ###-------- Output table of predicted probabilities +/- CIs --------###
    write.csv(plotfinal[,c("pred","lwr","upr","mgmtvar")], "overall probabilities and CIs.csv", row.names=FALSE)
    
  }
  
  if (species) {
    
    plotdat <- list()
    
    for (i in 1:length(mod)) {
      
      # dataset to predict over is the same as the original dataset
      pred <- predict(mod[[i]], type="response", re.form=NA)
      pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
      
      plotdat[[i]] <- unique(data.frame(pred, pred.CI, species=moddat[[i]]$species, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]])))
      
    }
    
    plotfinal <- do.call(rbind, plotdat)
    
    pch <- data.frame(species=levels(plotfinal$species), pch=c(21,22,24,25,21,22,24), col=sample(grey(seq(from=0,to=1,length.out = 7)), 7))
    plotfinal <- merge(plotfinal,pch, by="species")
    
    plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
    
    plotfinal$rowid <- 1:nrow(plotfinal)
    
    xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
    xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
    xloc.divide <- xloc.divide[-length(xloc.divide)]
    
    ###-------- Output plot --------###
    
    if (alphalevel==0.05) {
      png("species-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)
    } else {png("species-specific model results_84CIs.png", res=300, height=12, width=30, units="in", pointsize=20)}
    
    par(mar=c(7,6,3,2))
    
    x <- c(1:nrow(plotfinal))
    
    plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(2,50), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
    axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
    abline(v=xloc.divide, lty=3, lwd=1.5)
    text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
    arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
    arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
    title(xlab="Management intervention evaluated", cex.lab=1.5, font=2, line=5)
    title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
    abline(h=0.5, lty=3, lwd=2)
    
    
    legend(0.5,1.1, legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
    
    dev.off()
    
    ###-------- Output table of predicted probabilities +/- CIs --------###
    if (alphalevel==0.05) {
      write.csv(plotfinal[,c("species","pred","lwr","upr","mgmtvar")], "species-specific probabilities and CIs.csv", row.names=FALSE)
    } else {write.csv(plotfinal[,c("species","pred","lwr","upr","mgmtvar")], "species-specific probabilities and 84CIs.csv", row.names=FALSE)}
    
    
    
  }
  
}

if (metric) {
  
  plotdat <- list()
  
  for (i in 1:length(mod)) {
    
    # dataset to predict over is the same as the original dataset
    pred <- predict(mod[[i]], type="response", re.form=NA)
    pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
    
    plotdat[[i]] <- unique(data.frame(pred, pred.CI, metric=moddat[[i]]$metric, mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]])))
    
  }
  
  plotfinal <- do.call(rbind, plotdat)
  
  pch <- data.frame(metric=levels(plotfinal$metric), pch=c(21,22,24), col=sample(grey(seq(from=0,to=1,length.out = 3)), 3))
  plotfinal <- merge(plotfinal,pch, by="metric")
  
  plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
  
  plotfinal$rowid <- 1:nrow(plotfinal)
  
  xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
  xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
  xloc.divide <- xloc.divide[-length(xloc.divide)]
  
  ###-------- Output plot --------###
  
  if (alphalevel==0.05) {
    png("metric-specific model results.png", res=300, height=12, width=30, units="in", pointsize=20)
  } else {png("metric-specific model results_84CIs.png", res=300, height=12, width=30, units="in", pointsize=20)}
  
  par(mar=c(7,6,3,2))
  
  x <- c(1:nrow(plotfinal))
  
  plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(1.5,26), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n")
  axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
  abline(v=xloc.divide, lty=3, lwd=1.5)
  text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
  arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
  arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
  title(xlab="Management intervention evaluated", cex.lab=1.5, font=2, line=5)
  title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
  abline(h=0.5, lty=3, lwd=2)
  
  
  legend(0.5,1.1, legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)
  
  dev.off()
  
  ###-------- Output table of predicted probabilities +/- CIs --------###
  if (alphalevel==0.05) {
    write.csv(plotfinal[,c("metric","pred","lwr","upr","mgmtvar")], "metric-specific probabilities and CIs.csv", row.names=FALSE)
  } else {write.csv(plotfinal[,c("metric","pred","lwr","upr","mgmtvar")], "metric-specific probabilities and 84CIs.csv", row.names=FALSE)}
  
  
  
}


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################