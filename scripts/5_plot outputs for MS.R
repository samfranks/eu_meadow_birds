#######################################################################
#
#     Step 5: EU meadow birds meta-analysis - output results for J Appl Ecol submission
#
#######################################################################

# Samantha Franks
# 4 Aug 2016

set.seed(2)

# =================================  SET LOGIC STATEMENTS  ====================

# default to plot when all are FALSE is results from overall analysis (0a)
species <- FALSE # plot the species-specific model results (0b)
metric <- FALSE # plot the metric-specific model results (0c)
habitat <- FALSE # plot the habitat-specific model results (0d)


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
outputwd <- paste(parentwd, "output/submission", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)

# ==============================  SET OUTPUT DIRECTORY  ===================================

setwd(outputwd)

# ================================    MANAGEMENT VARIABLE NAMES    ===============================


mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")



# ==============================  FIGURE 1  ===================================

# Fig 1a = AES, site protection (Analysis 0a)
# Fig 1b = basic vs higher AES (Analysis 0a)
# Fig 1c = combined effects of AES and site protection (Analysis 2a)


png("Fig1_AES_site protection.png", res=300, height=15, width=17, units="in", pointsize=24)

par(mfrow=c(2,2), oma=c(3,5,1,1))


# ==== Fig 1a and 1b ====


# -------    Load data and model   -----------


moddat <- readRDS(paste(workspacewd, "model dataset_0a.rds", sep="/")) #[c("AE","reserve.desig")]
mod <- readRDS(paste(workspacewd, "models_0a_lme4.rds", sep="/")) #[c("AE","reserve.desig")]


# -------   Produce plotting dataset predictions   ---------

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




### ---- Fig 1a: Policy level interventions plot ----


par(mar=c(3,2,2,2))

plotfinal <- fig1a

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("a)", side=3, adj=0)


### ---- Fig 1b: AES level interventions plot ----


par(mar=c(3,2,2,2))

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

mtext("b)", side=3, adj=0)

mtext("Intervention", side=1, outer=TRUE, line=1.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)


# ==== Fig 1c ====


# -------    Load data and model   -----------

plotmod <- readRDS(file=paste(workspacewd, "models_2a_method 1.rds", sep="/")) # model to plot results
origdat <- readRDS(file=paste(workspacewd, "model dataset_2a_method 1.rds", sep="/")) # original dataset

unique.mgmtvars <- unique(origdat$AE.reserve)

newdat <- data.frame(AE.reserve=rep(unique.mgmtvars, times=length(levels(origdat$species))), species=rep(levels(origdat$species), each=length(unique.mgmtvars)))

pred <- predict(plotmod, newdat, type="response", re.form=NA)
pred.CI <- easyPredCI(plotmod, newdat)
fits <- data.frame(newdat, pred, pred.CI)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(AE.reserve=fits$AE.reserve), mean)

plotdat <- sum.fits
plotdat <- merge(plotdat, unique(origdat[,c("AE","reserve.desig","AE.reserve")]), by="AE.reserve")
plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE),]



### ---- Fig 1c: AES*site protection in combination plot (Method 1 creating a new combined AE-reserve variable) ----

par(mar=c(3,2,2,2))

x <- c(1:nrow(plotdat))

plotdat$pch <- c(16,16,16)
# plotdat$pch <- c(1,2,15,16,17)

plot(plotdat$pred~x, pch=plotdat$pch, cex=2, ylim=c(0,1), xlim=c(0.8,3.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
# axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES only","site protection \nonly", "AES + \nsite protection"), cex=1)
# text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
# title(xlab="Intervention combination", cex.lab=1.2, line=4.5)

mtext("c)", side=3, adj=0)

dev.off()



# ==============================  FIGURE 2  ===================================

# Fig 2a = AES, site protection by species (Analysis 0b)
# Fig 2b = AES level by species (Analysis 0b)
# Fig 2c = AES, site protection by metric (Analysis 0c)
# Fig 2d = AES level by metric (Analysis 0c)



png("Fig2_AES_site protection_species_metric.png", res=300, height=15, width=20, units="in", pointsize=20)

par(mfcol=c(2,2), oma=c(3,5,3,1))


# ==== Fig 2a and 2b - SPECIES ====


# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_0b.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_0b_blme.rds", sep="/"))


# -------   Produce plotting dataset predictions   ---------

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

fig2a <- do.call(rbind, plotdat[c("AE","reserve.desig")])
fig2b <- do.call(rbind, plotdat[c("AE.level")])



### ---- Fig 2a: Policy level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2a
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
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("a)", side=3, adj=0)

legend("topleft", legend=levels(plotfinal$species), pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)

mtext("Intervention", side=1, outer=TRUE, line=1.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)


### ---- Fig 2b: AES level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2b
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
abline(h=0.05, lty=3, lwd=2)

mtext("b)", side=3, adj=0)




# ==== Fig 2c and 2d - METRIC ====


# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_0c.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_0c_blme.rds", sep="/"))


# -------   Produce plotting dataset predictions   ---------

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

fig2c <- do.call(rbind, plotdat[c("AE","reserve.desig")])
fig2d <- do.call(rbind, plotdat[c("AE.level")])



### ---- Fig 2c: Policy level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2c
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
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2)

legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)

mtext("c)", side=3, adj=0)


### ---- Fig 2d: AES level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2d
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
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2)

mtext("d)", side=3, adj=0)


dev.off()



# ==============================  FIGURE 3  - SPECIFIC MANAGEMENT INTERVENTIONS ===================================

# Fig 3a = specific measures overall (Analysis 0a)
# Fig 3b = specific measures by species (Analysis 0b)
# Fig 3c = specific measures by metric (Analysis 0c)



png("Fig3_specifc management_overall_species_metric.png", res=300, height=50, width=50, units="in", pointsize=60)

par(mfrow=c(3,1), oma=c(3,5,1,1))


# ==== Fig 3a - OVERALL ====


# -------    Load data and model   -----------


moddat <- readRDS(paste(workspacewd, "model dataset_0a.rds", sep="/")) #[c("AE","reserve.desig")]
mod <- readRDS(paste(workspacewd, "models_0a_lme4.rds", sep="/")) #[c("AE","reserve.desig")]


# -------   Produce plotting dataset predictions   ---------

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

fig3a <- do.call(rbind, plotdat[4:9])


### ---- Fig 3a: Specific management interventions plot ----


par(mar=c(4,3,2,2))

plotfinal <- fig3a

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \napplied", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/\npesticides \napplied","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied", "water \nreduced"))
# text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE, labels=c("AES","basic-level \n AES","higher-level \n AES","nature reserve/ \n designation", "mowing applied", "mowing reduced", "grazing applied", "grazing reduced", "fertiliser/pesticides \n applied","fertiliser/pesticides \n reduced","nest protection \n applied","predator control \n applied","water \n applied", "water \n reduced"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("a)", side=3, adj=0)

mtext("Intervention", side=1, outer=TRUE, line=1.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)



# ==== Fig 3b - SPECIES ====

# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_0b.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_0b_blme.rds", sep="/"))


# -------   Produce plotting dataset predictions   ---------

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

fig3b <- do.call(rbind, plotdat[4:9])


### ---- Fig 3b: Specific management interventions plot - SPECIES ----


par(mar=c(4,3,2,2))

plotfinal <- fig3b
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
text(xloc.mgmtvars, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \nreduced", "grazing \napplied", "grazing \nreduced","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2)

mtext("b)", side=3, adj=0)

legend("topleft", legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)




# ==== Fig 3b - METRIC ====

# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_0c.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_0c_blme.rds", sep="/"))


# -------   Produce plotting dataset predictions   ---------

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

fig3c <- do.call(rbind, plotdat[4:9])


### ---- Fig 3c: Specific management interventions plot - METRIC ----


par(mar=c(4,3,2,2))

plotfinal <- fig3c
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
# axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
abline(v=xloc.divide, lty=3, lwd=1.5)
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \napplied","mowing \nreduced", "grazing \napplied", "grazing \nreduced","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied", "water \nreduced"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2)

mtext("c)", side=3, adj=0)

legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE)



dev.off()




# 
# 
# # load model data & models: use blme models for all 0a-d analyses
# 
# if (species & !habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0b.rds", sep="/"))
# if (!species & !habitat & metric) moddat <- readRDS(paste(workspacewd, "model dataset_0c.rds", sep="/"))
# if (!species & habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0d.rds", sep="/"))
# if (bias & !species & !habitat & !metric) moddat <- readRDS(paste(workspacewd, "model dataset_0a_bias.rds", sep="/"))
# 
# if (bias) moddat <- moddat[-8] # remove predator control model (huge variance on the random effects, and SEs on fixed effects are very large)
# 
# # load models
# if (!bias & !species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0a_lme4.rds", sep="/"))
# if (species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0b_blme.rds", sep="/"))
# if (!species & !habitat & metric) mod <- readRDS(paste(workspacewd, "models_0c_blme.rds", sep="/"))
# if (!species & habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0d_blme.rds", sep="/"))
# if (bias & !species & !habitat & !metric) mod <- readRDS(paste(workspacewd, "models_0a_lme4_bias.rds", sep="/"))
# 
# if (bias) mod <- mod[-8] # remove predator control model (huge variance)
# 
# #=================================  OUTPUT PARAMETER TABLES  ===============================
# 
# setwd(outputwd)