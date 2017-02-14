#######################################################################################
#
#     Step 7: EU meadow birds analysis - output results for J Appl Ecol revision
#
#######################################################################################

# Samantha Franks
# 4 Aug 2016
# 22 Dec 2016

set.seed(2)

# =================================  SET LOGIC STATEMENTS  ====================

alphalevel <- 0.05
successlevel <- 0.05

# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","tidyr")

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

# ==============================  SET OUTPUT DIRECTORY  ===================================

setwd(outputwd)

# ================================    MANAGEMENT VARIABLE NAMES    ===============================


mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","water","nest.protect","predator.control")


# ==============================  FIGURE 1 - overall success of interventions ===================================

png("Fig1_overall success of interventions.png", res=300, height=15, width=18, units="in", pointsize=22)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(oma=c(3,5,1,1))


# -------    Load data and model   -----------


moddat <- readRDS(paste(workspacewd, "model dataset_analysis 1.rds", sep="/")) 
mod <- readRDS(paste(workspacewd, "models_analysis 1_lme4.rds", sep="/"))


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
fig1c <- do.call(rbind, plotdat[c("mowing","grazing","fertpest","water","nest.protect","predator.control")])


# -------- Plot policy (AES, site protection) interventions ---------

par(mar=c(3,3,3,4))

plotfinal <- fig1a

x <- c(2:3)

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(1.5,3.5))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site\nprotection"))
# text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("a)", side=3, adj=0, line=1)



# -------- Fig 1b: plot AES levels ---------

par(mar=c(3,3,3,4))

plotfinal <- fig1b

x <- c(2:3)

plot(plotfinal$pred~x, xlim=c(1.5,3.5), ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n")
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("basic\nAES","higher\nAES"))
# text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("b)", side=3, adj=0, line=1)



# -------- Fig 1c: Plot management interventions ---------

par(mar=c(5,3,3,2))

plotfinal <- fig1c

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x), max(x)))
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing\napplied", "mowing\nreduced", "grazing\napplied", "grazing\nreduced", "agro-\nchemicals\napplied","agro-\nchemicals\nreduced","water\napplied","water\nreduced","nest\nprotection\napplied", "predator\ncontrol\napplied"))
# text(x, par("usr")[3]-0.06, srt = 30, pos=1, xpd = TRUE
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2)

mtext("c)", side=3, adj=0, line=1)

mtext("Intervention", side=1, outer=TRUE, line=1.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)

dev.off()



# ==============================  FIGURE 2  - species-specific success of interventions ===================================

# Fig 2a = AES, site protection by species (Analysis 2a)
# Fig 2b = AES level by species (Analysis 2a)
# Fig 2c = management by species (Analysis 2a)


png("Fig2_intervention success_species.png", res=300, height=15, width=22, units="in", pointsize=22)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(oma=c(3,5,3,1))


# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_analysis 2a.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_analysis 2a_blme.rds", sep="/"))


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
fig2c <- do.call(rbind, plotdat[c("mowing","grazing","fertpest","water","nest.protect","predator.control")])



### ---- Fig 2a: Policy level interventions plot - SPECIES ----


par(mar=c(3,3,2,2))

plotfinal <- fig2a
maxspecies <- levels(do.call(rbind, plotdat)$species)
n <- length(maxspecies)

set.seed(2)
pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(grey(seq(from=0.2,to=1,length.out = n)), replace=TRUE, n))
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

mtext("a)", side=3, adj=0, line=1)

# legend("topleft", legend=levels(fig2a$species), pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.05))

legend("topleft", legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.05,-0.1))


### ---- Fig 2b: AES level interventions plot - SPECIES ----


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

mtext("b)", side=3, adj=0, line=1)




### ---- Fig 2c: Specific management interventions plot - SPECIES ----


par(mar=c(5,3,3,2))

plotfinal <- fig2c
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
text(xloc.mgmtvars, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing\nreduced", "grazing\napplied", "grazing\nreduced","agro-\nchemicals\nreduced","water\napplied","nest\nprotection\napplied","predator \ncontrol\napplied"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2)

mtext("c)", side=3, adj=0, line=1)

mtext("Intervention", side=1, outer=TRUE, line=1, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)


dev.off()




# ==============================  FIGURE 3  - metric-specific success of interventions ===================================


# Fig 3a = AES, site protection by metric (Analysis 2b)
# Fig 3b = AES level by metric (Analysis 2b)
# Fig 3c = management by metric (Analysis 2b)


png("Fig3_intervention success_metric.png", res=300, height=15, width=22, units="in", pointsize=22)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(oma=c(3,5,3,1))


# -------    Load data and model   -----------

moddat <- readRDS(paste(workspacewd, "model dataset_analysis 2b.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_analysis 2b_blme.rds", sep="/"))



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

fig3a <- do.call(rbind, plotdat[c("AE","reserve.desig")])
fig3b <- do.call(rbind, plotdat[c("AE.level")])
fig3c <- do.call(rbind, plotdat[c("mowing","grazing","fertpest","water","nest.protect","predator.control")])



### ---- Fig 3a: Policy level interventions plot - METRIC ----


par(mar=c(3,3,2,2))

plotfinal <- fig3a
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

# legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.05))

legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.04,-0.1))

mtext("a)", side=3, adj=0, line=1)


### ---- Fig 3b: AES level interventions plot - METRIC ----


par(mar=c(3,3,2,2))

plotfinal <- fig3b
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

mtext("b)", side=3, adj=0, line=1)



### ---- Fig 3c: Specific management interventions plot - METRIC ----

par(mar=c(5,3,3,2))

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
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing\napplied","mowing\nreduced", "grazing\napplied", "grazing\nreduced","agro-\nchemicals\nreduced","water\napplied", "water\nreduced","nest\nprotection\napplied","predator \ncontrol\napplied"))

arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="grey30")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="grey30")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2)

mtext("c)", side=3, adj=0, line=1)

mtext("Intervention", side=1, outer=TRUE, line=1, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)



dev.off()



# ==========================  FIGURE 4  - success of AES*site protection combined  ============================



# -------    Load data and model   -----------

plotmod <- readRDS(file=paste(workspacewd, "models_analysis 3a_method 1.rds", sep="/")) # model to plot results
origdat <- readRDS(file=paste(workspacewd, "model dataset_analysis 3a_method 1.rds", sep="/")) # original dataset

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



# ----- Fig 4 ----------

png("Fig4_intervention success_AES x site protection.png", res=300, height=12, width=14, units="in", pointsize=20)

par(oma=c(3,5,0,1))
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
text(x, par("usr")[3]-0.01, srt = 0, pos=1, xpd = TRUE, labels=c("AES only","site protection \nonly", "AES +\nsite protection"), cex=1)
# text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
# title(xlab="Intervention combination", cex.lab=1.2, line=4.5)

mtext("Intervention", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)

dev.off()



# =========================   SUPPLEMENTARY INFORMATION FIGURES  ========================



# ===============  Fig S6a - Average invervention combination success rates  ==================


setwd(outputwd)

png("FigS6a_specific combination intervention success_average.png", res=300, height=15, width=25, units="in", pointsize=20)

plotdat <- readRDS(file=paste(workspacewd, "analysis 3b_average intervention plotting dataset.rds", sep="/"))

par(mar=c(3,2,2,2))
par(oma=c(3,4,1,1))

x <- c(1:nrow(plotdat)) # gives an x axis to plot against

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(0,1), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05) # add error bars, lwr and upr to each prediction
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3], srt = 0, pos=1, xpd = TRUE, labels=c("mowing\napplied","mowing\nreduced","grazing\napplied","grazing\nreduced","agro-\nchemicals\napplied","agro-\nchemicals\nreduced", "water\napplied", "water\nreduced","nest\nprotection\napplied","predator\ncontrol\napplied"), cex=1)

mtext("Intervention (controlling for other interventions used)", side=1, outer=TRUE, line=0.5, cex=1.2)
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5)

mtext("a)", side=3, adj=0, line=0.5, cex=1.2)


dev.off()




# ===============  Fig S6b - Ranked invervention combination success rates  ==================


setwd(outputwd)

png("FigS6b_specific combination intervention success.png", res=300, height=15, width=30, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(1,8,3,2))

plotdat <- readRDS(file=paste(workspacewd, "analysis 3b_all combinations plotting dataset.rds", sep="/"))

x <- c(1:nrow(plotdat))-0.5 # gives an x axis to plot against

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(0,1), xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(1,26))
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05) # add error bars, lwr and upr to each prediction
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2) # add a 'significance' line (what is the threshold for 'success'?)
abline(v=max(which(plotdat$sig=="N")), lty=3, lwd=2) # add a line dividing 'successful' vs 'unsuccessful' intervention combos
title(xlab="Intervention combination", cex.lab=1.2, font=2, line=0, xpd=TRUE)
title(ylab="Predicted probability of success\n(significant positive impact) ", cex.lab=1.2, font=2, line=3)

mtext("b)", side=3, adj=0, line=1, cex=1.2)


### Create a 'table' of intervention combinations to display below the plot showing the predicted success ###

y <- length(mgmtvars[4:9]):1 # how many interventions there are (will be labels down the y-axis of the table starting at the top and working down)

# create the table of coordinates for the table (centres of the grid cells for the table), which is 28 across (the number of different intervention combinations) x 6 down (the number of types of interventions)
x <- c(1:nrow(plotdat))
new.x <- rep(x, each=max(y))
new.y <- rep(y, times=max(x))
tab <- data.frame(x=new.x,y=new.y)
tab <- tab[order(tab$y, decreasing=TRUE),]
labs <- plotdat[,1:6]

labs.long <- gather(labs, intervention, level, mowing:predator.control) # convert the interventions from wide to long format using tidyr

# merge the table of x,y coordinates with the intervention labels
# replace all levels of 'none' with a blank
tab.filled <- data.frame(tab, labs.long)
tab.filled[4] <- apply(tab.filled[4], 2, function(x) {
  gsub("none", "", x)
})

# create a dataframe with x,y coordinates for the locations of the total number of interventions used
total.interventions <- data.frame(x, y=rep(0, times=length(x)), intervention=rep("", length(x)), sum=plotdat$num.interventions.sum)

# set up the plot of the 'table'
par(mar=c(1,8,1,2)) # give a big left margin to write the intervention names; margin for above plot needs to be exactly the same so everything lines up
plot(tab, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0.2,6.8), xlim=c(1,26)) # plot a completely blank plot, with ylims and xlims only (will draw the table using abline; use exactly the same xlims as the actual data plot above so that everything lines up
abline(v=c(min(x)-1,x,max(x)), h=c(7,y,0), lty=1) # draw the lines of the table (really the 'inner lines'), but adding 1 extra to each of top, bottom and sides to complete the outline of the table (otherwise will be inner lines only)
abline(v=max(which(plotdat$sig=="N")), lty=1, lwd=4) # thick line showing division between 'successful' and 'unsuccessful' combinations of interventions
axis(2, y+0.5, labels=c("mowing","grazing","agrochemicals",  "water", "nest protection","predator control"), las=1, cex.axis=1, font=2,tick=FALSE) # draw the labels for the rows, from top to bottom
text(tab$x-0.5, tab$y+0.5, labels=ifelse(tab.filled$level=="applied", "\U2191", ifelse(tab.filled$level=="reduced", "\U2193", tab.filled$level)), cex=2) # fill in the values of the grid cells, but if an intervention was applied then use a unicode 'up' arrow, and if it was reduced than use a down arrow
text(tab$x-0.5, 0.5, labels=total.interventions$sum, cex=1) # add the total number of interventions to the bottom 'row' of the table


dev.off()




# ===============  FIGURE S7 - Failed intervention probability ==================


setwd(outputwd)

moddat <- readRDS(paste(workspacewd, "model dataset_analysis 1b.rds", sep="/"))
mod <- readRDS(paste(workspacewd, "models_analysis 1b_lme4.rds", sep="/"))

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

plotfinal <- do.call(rbind, plotdat)


###-------- Output plot --------###

png("FigS7_overall model results_failures.png", res=300, height=12, width=28, units="in", pointsize=20)
par(mar=c(7,6,2,2))

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n")
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.06, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic\nAES","higher\nAES","site\nprotection", "mowing\napplied", "mowing\nreduced", "grazing\napplied", "grazing\nreduced", "agro-\nchemicals\napplied","agro-\nchemicals\nreduced","water\napplied", "water\nreduced","nest\nprotection\napplied"))
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05)
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05)
title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
title(ylab="Predicted probability of failure\n(significant negative impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2)
abline(v=4.5, lty=3, lwd=2)

dev.off()




