#######################################################
#
#     Step 1: EU meadow birds meta-analysis - DATA SUMMARY
#
#######################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================
# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","lme4","tidyr")

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
outputwd <- paste(parentwd, "output/presentation", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


#=================================  LOAD DATA  ===============================

source(paste(scriptswd, "2_meta-analysis_data preparation.R", sep="/"))



#=================================  SUMMARY STATISTICS  ===============================

# total number of studies in dataset
num.studies <- length(unique(dat$reference))

#---------- Proportion of studies by country ----------

# summarize proportions of studies from different countries
# create based on a unique dataset of reference id and country
countrysum <- table(unique(dat[,c("reference","country")])$country)
countrysum.prop <- countrysum/num.studies


png(paste(outputwd, "summary_proportion of studies by country species.png", sep="/"), res=300, height=18, width=15, units="in", pointsize=20, bg="transparent")

# png(paste(outputwd, "summary_proportion of studies by country.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(6,5,2,1))
# x <- barplot(countrysum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.6), xaxt="n")
# text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c(names(countrysum.prop)[-which(names(countrysum.prop) %in% "United Kingdom")], "United \n Kingdom"))
# title(xlab="Country", font=2, cex.lab=1.2, line=4.5)
# title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
# text(x, countrysum.prop+0.02, countrysum) # sample sizes for each country

# dev.off()

x <- barplot(countrysum.prop, space=0.1, las=1, col="white", xaxt="n", ylim=c(0,0.6), col.axis="white", fg="white")
text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c(names(countrysum.prop)[-which(names(countrysum.prop) %in% "United Kingdom")], "United \n Kingdom"), col="white")
title(xlab="Country", font=2, cex.lab=1.2, line=4, col.lab="white")
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3, col.lab="white")
text(x, countrysum.prop+0.02, countrysum, col="white") # sample sizes for each country



#---------- Proportion of studies by species ----------

speciessum <- table(unique(dat[,c("reference","species")])$species)
speciessum.prop <- speciessum/num.studies

# png(paste(outputwd, "summary_proportion of studies by species.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(speciessum.prop, space=0.1, las=1, col="white", ylim=c(0,1), xaxt="n", col.axis="white", fg="white")
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels = c( "black-tailed \n godwit",names(speciessum.prop)[-which(names(speciessum.prop) %in% "black-tailed godwit")]), col="white")
title(xlab="Species", font=2, cex.lab=1.2, line=4.5, col.lab="white")
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3, col.lab="white")
text(x, speciessum.prop+0.02, speciessum, col="white") # sample sizes for each species

dev.off()



#---------- Proportion of studies by overall metric ----------

# create an ordered factor of metrics for this summary (particularly for productivity so it's nest, chick, nest+chick)
# don't want to change this factor to ordered for the analysis though, because it will be treated differently than a regular factor in the model

d1.metric <- dat
d1.metric$new.metric <- factor(d1.metric$new.metric, c("abundance/occupancy","abundance/occupancy change","productivity","recruitment","survival"), ordered=TRUE)

metricsum <- table(unique(d1.metric[,c("reference","new.metric")])$new.metric)
metricsum.prop <- metricsum/num.studies


#---------- Proportion of studies by specific metric (for productivity)----------

# create an ordered factor of metrics for this summary (particularly for productivity so it's nest, chick, nest+chick)

d1.metric <- dat

# pool abundance/occupancy and change metrics
# show different productivity metrics separately
d1.metric$prod.metric <- with(d1.metric, ifelse(d1.metric$overall.metric=="abundance" | d1.metric$overall.metric=="occupancy", "abundance/occupancy", ifelse(d1.metric$overall.metric=="abundance change" | d1.metric$overall.metric=="occupancy change", "abundance/occupancy change", as.character(overall.metric))))

# create ordered factor
d1.metric$prod.metric <- factor(d1.metric$prod.metric, c("abundance/occupancy","abundance/occupancy change","productivity (nest level)", "productivity (chick level)", "productivity (nest + chick)","recruitment","survival"), ordered=TRUE)

metricsumprod <- table(unique(d1.metric[,c("reference","prod.metric")])$prod.metric)
metricsumprod.prop <- metricsumprod/num.studies

png(paste(outputwd, "summary_proportion of studies by productivity metric.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20, bg="transparent")

par(mar=c(6,5,2,1))
x <- barplot(c(metricsumprod.prop[1:2], metricsum.prop["productivity"], metricsumprod.prop[3:7]), space=0.1, las=1, col="white", ylim=c(0,0.8), xaxt="n", col.axis="white", fg="white")
text(x, par("usr")[3]-0.01, srt = 0, pos=1, xpd = TRUE, labels = c("abundance/\noccupancy","abundance/\noccupancy \nchange", "productivity \n(all)", "productivity \n(nest)", "productivity \n(chick)", "productivity \n(nest+chick)","recruitment","survival"), col="white")
title(xlab="Study metric", font=2, cex.lab=1.2, line=4.5, col.lab="white")
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3, col.lab="white")
text(x, c(metricsumprod.prop[1:2]+0.02, metricsum.prop["productivity"]+0.02, metricsumprod.prop[3:7]+0.02), c(metricsumprod[1:2], metricsum["productivity"], metricsumprod[3:7]), col="white") # sample sizes for each metric

dev.off()

#---------- Proportion of studies examining the effect of an intervention type ----------

### determine how many studies evaluated the effect of the intervention (repeat for each intervention)
# create blank objects
intervensum <- numeric()
intervensum.prop <- numeric()
intervensum.level <- list()
intervensum.level.prop <- list()

# mgmtvars to evluate, minus AE.level
eval.mgmtvars <- mgmtvars[-which(mgmtvars %in% "AE.level")]

for (i in 1:length(mgmtvars)) {
  
  # number of unique cases (i.e. unique studies) where mgmtvar level != 'none'
  x <- unique(dat[,c("reference",mgmtvars[i])]) # unique references and levels of the intervention
  y <- x[x[mgmtvars[i]] != "none",] # remove the cases where intervention not evaluated
  intervensum[i] <- length(unique(y$reference)) # the number of studies that evaluated the intervention
  intervensum.prop[i] <- intervensum[i]/num.studies
  intervensum.level[[i]] <- table(y[,mgmtvars[i]])
  intervensum.level.prop[[i]] <- intervensum.level[[i]]/num.studies
}

names(intervensum) <- mgmtvars
names(intervensum.prop) <- mgmtvars
names(intervensum.level) <- mgmtvars
names(intervensum.level.prop) <- mgmtvars

# remove AE level from being plotted
intervensum.prop <- intervensum.prop[-which(names(intervensum.prop) %in% "AE.level")]
intervensum <- intervensum[-which(names(intervensum) %in% "AE.level")]

png(paste(outputwd, "summary_proportion of studies by intervention.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20, bg="transparent")

par(mar=c(6,5,2,1))
x <- barplot(intervensum.prop, space=0.1, las=1, col="white", ylim=c(0,0.6), xaxt="n", col.axis="white", fg="white")
text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c("AES","site\nprotection", "mowing","grazing","fertiliser/\npesticides","nest\nprotection","predator\ncontrol","water\nmanagement"), col="white")
title(xlab="Intervention", font=2, cex.lab=1.2, line=4.5, col.lab="white")
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3, col.lab="white")
text(x, intervensum.prop+0.02, intervensum, col="white") # sample sizes for each intervention type

dev.off()


#=================================  EFFECTIVENESS OF INTERVENTIONS  ===============================


set.seed(2)

# =================================  SET LOGIC STATEMENTS  ====================

# default to plot when all are FALSE is results from overall analysis (0a)
species <- FALSE # plot the species-specific model results (0b)
metric <- FALSE # plot the metric-specific model results (0c)
habitat <- FALSE # plot the habitat-specific model results (0d)


alphalevel <- 0.05
successlevel <- 0.05


# =================================  LOAD FUNCTIONS =================================


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


# ==============================  SET OUTPUT DIRECTORY  ===================================


setwd(outputwd)

# ================================    MANAGEMENT VARIABLE NAMES    ===============================


mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")


# ==============================  FIGURE 1  ===================================

# Fig 1a = AES, site protection (Analysis 0a)
# Fig 1b = basic vs higher AES (Analysis 0a)
# Fig 1c = combined effects of AES and site protection (Analysis 2a)


png("Fig1_AES_site protection.png", res=300, height=10, width=18, units="in", pointsize=24, bg="transparent")

par(mfrow=c(1,2), oma=c(3,5,1,1))


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

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5), col.axis="white", fg="white", col="white")
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
abline(h=successlevel, lty=3, lwd=2, col="white")



### ---- Fig 1b: AES level interventions plot ----


par(mar=c(3,2,2,2))

plotfinal <- fig1b

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5), col.axis="white", fg="white", col="white")
# axis(1, x, labels=rep("",nrow(plotfinal)), tick=TRUE)
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES","higher-level AES"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1, line=3)
abline(h=successlevel, lty=3, lwd=2, col="white")

mtext("Intervention", side=1, outer=TRUE, line=1.5, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")

dev.off()


# # ==== Fig 1c ====
# 
# 
# # -------    Load data and model   -----------
# 
# plotmod <- readRDS(file=paste(workspacewd, "models_2a_method 1.rds", sep="/")) # model to plot results
# origdat <- readRDS(file=paste(workspacewd, "model dataset_2a_method 1.rds", sep="/")) # original dataset
# 
# unique.mgmtvars <- unique(origdat$AE.reserve)
# 
# newdat <- data.frame(AE.reserve=rep(unique.mgmtvars, times=length(levels(origdat$species))), species=rep(levels(origdat$species), each=length(unique.mgmtvars)))
# 
# pred <- predict(plotmod, newdat, type="response", re.form=NA)
# pred.CI <- easyPredCI(plotmod, newdat)
# fits <- data.frame(newdat, pred, pred.CI)
# 
# # produce mean population level prediction for interventions across species
# sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(AE.reserve=fits$AE.reserve), mean)
# 
# plotdat <- sum.fits
# plotdat <- merge(plotdat, unique(origdat[,c("AE","reserve.desig","AE.reserve")]), by="AE.reserve")
# plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE),]
# 
# 
# 
# ### ---- Fig 1c: AES*site protection in combination plot (Method 1 creating a new combined AE-reserve variable) ----
# 
# par(mar=c(3,2,2,2))
# 
# x <- c(1:nrow(plotdat))
# 
# plotdat$pch <- c(16,16,16)
# # plotdat$pch <- c(1,2,15,16,17)
# 
# plot(plotdat$pred~x, pch=plotdat$pch, cex=2, ylim=c(0,1), xlim=c(0.8,3.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
# arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
# arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
# abline(h=0.05, lty=3, lwd=2)
# # axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
# # axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
# text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("AES only","site protection \nonly", "AES + \nsite protection"), cex=1)
# # text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# # text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# # text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
# # title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.2, line=3)
# # title(xlab="Intervention combination", cex.lab=1.2, line=4.5)
# 
# mtext("c)", side=3, adj=0)
# 
# dev.off()

# ==============================  FIGURE 2 - by species ===================================



png("Fig2_AES_site protection_species.png", res=300, height=10, width=18, units="in", pointsize=20, bg="transparent")

par(mfcol=c(1,2), oma=c(3,5,3,1))


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
pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(c("blue","darkmagenta","orange"), replace=TRUE, n))
pch
plotfinal <- merge(plotfinal,pch, by="species")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(min(plotfinal$rowid):max(plotfinal$rowid))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
# axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2, col="white")

legend("topleft", legend=levels(plotfinal$species), pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.05), text.col="white", col="white")

mtext("Intervention", side=1, outer=TRUE, line=0.5, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")


### ---- Fig 2b: AES level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2b
maxspecies <- levels(do.call(rbind, plotdat)$species)
n <- length(maxspecies)

set.seed(2)
pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(c("blue","darkmagenta","orange"), replace=TRUE, n))
pch
plotfinal <- merge(plotfinal,pch, by="species")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(min(plotfinal$rowid):max(plotfinal$rowid))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
# axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
# text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","basic-level \nAES","higher-level \nAES","nature reserve/ \ndesignation", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/ \npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","more water \napplied"))
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
abline(h=0.05, lty=3, lwd=2, col="white")


dev.off()



# ==============================  FIGURE 3 - by metric ===================================



png("Fig3_AES_site protection_metric.png", res=300, height=10, width=18, units="in", pointsize=20, bg="transparent")

par(mfcol=c(1,2), oma=c(3,5,3,1))



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
pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(c("blue","darkmagenta","orange"), n))
pch
plotfinal <- merge(plotfinal,pch, by="metric")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid)+0.2), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
# axis(1, xloc.mgmtvars, labels=rep("",length(xloc.mgmtvars)), tick=TRUE)
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("AES","site protection"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=successlevel, lty=3, lwd=2)

legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.05), text.col="white", col="white")

mtext("Intervention", side=1, outer=TRUE, line=0.5, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")

### ---- Fig 2d: AES level interventions plot ----


par(mar=c(3,3,2,2))

plotfinal <- fig2d
maxmetric <- levels(do.call(rbind, plotdat)$metric)
n <- length(maxmetric)

set.seed(2)
pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(c("blue","darkmagenta","orange"), n))
pch
plotfinal <- merge(plotfinal,pch, by="metric")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(min(plotfinal$rowid):max(plotfinal$rowid))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid)-0.5,max(plotfinal$rowid)+0.2), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("basic AES","higher AES"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
abline(h=0.05, lty=3, lwd=2, col="white")



dev.off()


# ==============================  FIGURE 4 - habitat management ===================================


png("Fig4_specifc management_overall.png", res=300, height=10, width=18, units="in", pointsize=20, "transparent")

par(oma=c(3,5,1,0))


# ==== OVERALL ====


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


par(mar=c(5,3,3,1))

plotfinal <- fig3a

x <- c(1:nrow(plotfinal))

plot(plotfinal$pred~x, ylim=c(0,1), pch=16, cex=2, xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(min(x)-0.5, max(x)+0.5), col.axis="white", fg="white", col="white")
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \napplied", "mowing \nreduced", "grazing \napplied", "grazing \nreduced", "fertiliser/\npesticides \napplied","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied", "water \nreduced"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
abline(h=successlevel, lty=3, lwd=2, col="white")


mtext("Intervention", side=1, outer=TRUE, line=0.5, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")

dev.off()



# ====  SPECIES ====


png("Fig5_specifc management_species.png", res=300, height=9, width=18, units="in", pointsize=20, "transparent")

par(oma=c(3,5,1,0))

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


par(mar=c(5,3,3,0))

plotfinal <- fig3b
maxspecies <- levels(do.call(rbind, plotdat)$species)
n <- length(maxspecies)

set.seed(2)
pch <- data.frame(species=maxspecies, pch=rep(c(21,22,23,24,25),length.out=n), col=sample(c("blue","darkmagenta","orange"), replace=TRUE, n))

pch
plotfinal <- merge(plotfinal,pch, by="species")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$species),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(min(plotfinal$rowid):max(plotfinal$rowid))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
text(xloc.mgmtvars, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \nreduced", "grazing \napplied", "grazing \nreduced","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
# title(xlab="Intervention", cex.lab=1.5, font=2, line=5)
# title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
abline(h=0.05, lty=3, lwd=2, col="white")


legend("topleft", legend=pch$species, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.15), text.col="white", col="white")


mtext("Intervention", side=1, outer=TRUE, line=-1, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")

dev.off()



# ====  METRIC ====

png("Fig6_specifc management_metric.png", res=300, height=9, width=18, units="in", pointsize=20, "transparent")

par(oma=c(3,5,1,0))


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


par(mar=c(5,3,3,0))

plotfinal <- fig3c
maxmetric <- levels(do.call(rbind, plotdat)$metric)
n <- length(maxmetric)

set.seed(2)
pch <- data.frame(metric=maxmetric, pch=rep(c(21,22,23),length.out=n), col=sample(c("blue","darkmagenta","orange"), n))
pch
plotfinal <- merge(plotfinal,pch, by="metric")
plotfinal <- plotfinal[order(plotfinal$mgmtvar,plotfinal$metric),]
plotfinal$rowid <- 1:nrow(plotfinal)

xloc.mgmtvars <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), mean)$x
xloc.divide <- aggregate(plotfinal$rowid, list(plotfinal$mgmtvar), max)$x + 0.5
xloc.divide <- xloc.divide[-length(xloc.divide)]

x <- c(min(plotfinal$rowid):max(plotfinal$rowid))

plot(plotfinal$pred~x, ylim=c(0,1), xlim=c(min(plotfinal$rowid),max(plotfinal$rowid)), pch=plotfinal$pch, bg=as.character(plotfinal$col), cex=1.8, xaxt="n", xlab="", ylab="", las=1, bty="n", col.axis="white", fg="white", col="white")
abline(v=xloc.divide, lty=3, lwd=1.5, col="white")
text(xloc.mgmtvars, par("usr")[3]-0.05, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \napplied","mowing \nreduced", "grazing \napplied", "grazing \nreduced","fertiliser/\npesticides \nreduced","nest \nprotection \napplied","predator \ncontrol \napplied","water \napplied", "water \nreduced"), col="white")
arrows(x, plotfinal$pred, x, plotfinal$lwr, angle=90, length=0.05, col="white")
arrows(x, plotfinal$pred, x, plotfinal$upr, angle=90, length=0.05, col="white")
abline(h=0.05, lty=3, lwd=2, col="white")


legend("topleft", legend=pch$metric, pch=pch$pch, pt.bg=as.character(pch$col), pt.cex=1.2, bty="n", xpd=TRUE, inset=c(0.01,-0.2), text.col="white", col="white")


mtext("Intervention", side=1, outer=TRUE, at=0.35, line=-1, cex=1.2, col="white")
mtext("Predicted probability of success \n (significant positive impact)", side=2, outer=TRUE, cex=1.2, line=1.5, col="white")


dev.off()

