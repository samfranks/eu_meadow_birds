#################################################################
#
#     Step 4: EU meadow birds meta-analysis - models evaluating effect size
#
#################################################################

# Samantha Franks
# 27 April 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme","multcomp")

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
dat1 <- subset(dat0, select=c("reference","lit.type","score","country","study.length","habitat","species","overall.metric","metric","specific.metric","stan.metric","sample.size","analysis2","success","stan.calc","stan.metric.before","stan.metric.after","stan.effect.size","sig","effect.dir",mgmtvars))

# subset dataset to include records with effect sizes measured (whatever the significance of the effect)
dat2 <- subset(dat1, stan.effect.size!="" & stan.effect.size!="#DIV/0!")
dat2$stan.effect.size <- as.numeric(as.character(dat2$stan.effect.size))
dat2 <- droplevels(dat2)
hist(dat2$stan.effect.size)

###*** Long right-tail on the distribution, some very large positive effects

dat <- dat2

#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}

#=================================  ANALYSIS  ===============================

# can combine standardised abundance metrics - effect sizes comparable
# combine chick survival standardised metrics - effect sizes comparable

dat$new.stan.metric <- ifelse(grepl("chick survival", dat$stan.metric), "chick survival", ifelse(grepl("number of", dat$stan.metric), "abundance", as.character(dat$stan.metric)))

# metrics to test - other metrics have sample sizes which are too small
metrics <- c("abundance","multiplicative yearly slope","nest survival (Mayfield)")

###----  Nuisance variables, all metrics pooled ----###

nui.dat <- unique(dat[,c("reference","study.length","sample.size","analysis2","lit.type","stan.effect.size")])

# use unique dataset only to test effect sizes (replication because of multiple records of same study)
m.nui1 <- lme(stan.effect.size ~ study.length + sample.size + analysis2 + lit.type, random = ~1|reference, data=nui.dat)
summary(m.nui1)

# no significant effects of any variable on whole dataset (metrics pooled)

###----  Nuisance variables, individual metric subsets ----###

nui.dat <- unique(dat[,c("reference","study.length","sample.size","analysis2","lit.type","stan.effect.size","new.stan.metric")])

m.nui2 <- list()
m.nui2.dat <- list()

for (i in 1:length(metrics)) {
  
  print(metrics[i])
  mdat <- subset(nui.dat, new.stan.metric==metrics[i])
  mdat <- droplevels(mdat)
  m.nui2.dat[[i]] <- mdat
  
  print(nrow(mdat))
  
  m.nui2[[i]] <- lme(stan.effect.size ~ study.length + sample.size + analysis2 + lit.type, random = ~ 1|reference, data=mdat)
}

names(m.nui2) <- metrics
lapply(m.nui2, summary)

# abundance: no significant effects of any nuisance variables
# abundance change: no significant effects of any nuisance variables
# nest survival: significant effect of lit.type (negative effect of primary literature)

# tried to include 'score' as an additional variable but creates singularities in the model



#---------------------   Effect size on abundance, abundance change, nest survival metrics in relation to management intervention --------------------

mod.high <- list()
mod.spec <- list()
moddat <- list()

for (i in 1:length(metrics)) {
  
  print(metrics[i])
  mdat <- subset(dat, new.stan.metric==metrics[i])
  
  #   # check for sample sizes of different management treatment groups, remove small ones
  #   checksample <- table(mdat$mgmt)
  #   removelevels <- names(checksample)[checksample < 5]
  out <- list()
  for(j in 1:length(mgmtvars)) {
    out[[j]] <- table(dat[,mgmtvars[j]], dat$reference)
  }
  names(out) <- mgmtvars
  out
  
  # DFs for predator control (applied and reduced), fertpest, and water reduced are very small because they are each only tested in a single study so the model can't estimate the effect of those interventions, separately from the effect of the random effect study (checked this with Ali)
  # not gaining anything useful from this so remove these interventions and/or levels from the model
  if (metrics[i]=="abundance") {
    mdat <- subset(mdat, water!="reduced")
    mdat <- subset(mdat, stan.effect.size < 16) # remove large outlier
  }
  if (metrics[i]=="multiplicative yearly slope") {
    mdat <- subset(mdat, stan.effect.size > -0.4) # remove large outlier
    mdat <- subset(mdat, mowing!="reduced")
  }
  if (metrics[i]=="nest survival (Mayfield)") mdat <- subset(mdat, species!="curlew" & species!="dunlin" & AE.level!="basic") # control for too few studies assessing curlew or dunlin
  
  # removelevels <- c("predator.control applied", "predator.control reduced", "water reduced")
  
  # if (length(removelevels>0)) mdat <- subset(mdat, !grepl(paste(removelevels,collapse="|"), mdat$mgmt)) else mdat <- mdat
  
  mdat <- droplevels(mdat)
  moddat[[i]] <- mdat
  
  print(nrow(mdat))
  
  if (metrics[i]=="abundance") {
    
    mod.high[[i]] <- lme(stan.effect.size ~ AE.level + reserve.desig + species, random = ~1|reference, data=mdat) # no predator control, fertpest records that aren't none in abundance model
    mod.spec[[i]] <- lme(stan.effect.size ~ mowing + grazing + nest.protect + water + species, random = ~1|reference, data=mdat) # no predator control, fertpest records that aren't none in abundance dataset
    
    print(summary(mod.high[[i]]))
    print(summary(mod.spec[[i]]))
    
  }
  
  if (metrics[i]=="multiplicative yearly slope") {
    
    mod.high[[i]] <- lme(stan.effect.size ~ AE.level + reserve.desig + species, random = ~1|reference, data=mdat) # no predator control, fertpest records that aren't none in abundance model
    print(summary(mod.high[[i]]))
    
    mod.spec[[i]] <- c("for specific intervention model, adding more than one intervention creates singularities or leads to having studies which only test a single intervention and no other, creating weird DF estimation")
    
  }
  
  if (metrics[i]=="nest survival (Mayfield)") {
    
    mod.high[[i]] <- lme(stan.effect.size ~ AE.level + species, random = ~1|reference, data=mdat)
    print(summary(mod.high[[i]]))
    
    # mod.high[[i]] <- c("for high-level model, only enough data to robustly estimate parameter for AE higher; small DFs for AE basic and reserves")
    
    mod.spec[[i]] <- c("for specific intervention model, adding more than one intervention creates singularities or leads to having studies which only test a single intervention and no other, creating weird DF estimation")
    
  }
  
}

names(mod.high) <- metrics
names(mod.spec) <- metrics
names(moddat) <- metrics

mod.high

x <- glht(mod.high[[2]], linfct = mcp(AE.level = "Tukey"))
summary(x)

#---------------------   Effect size on nest/chick survival measures pooled --------------------

# run effect size for all productivity measures pooled (apart from proportion pairs with chicks since it's a very different metric), but with factor variable controlling for type of productivity measure if interventions have different effects on nests vs chicks etc

mdat.prod.high <- subset(dat, metric=="productivity" & new.stan.metric!="proportion pairs with chicks" & stan.effect.size < 30 & species!="curlew" & species!="dunlin" & AE.level!="basic")
mdat.prod.high <- droplevels(mdat.prod.high)
mdat.prod.high$new.stan.metric <- as.factor(mdat.prod.high$new.stan.metric)

# DF for reserve.desig too small even when in the model without AE.level
mod.prod.high <- lme(stan.effect.size ~ AE.level + species + new.stan.metric, random= ~1|reference, data=mdat.prod.high)
summary(mod.prod.high)

mdat.prod.spec <- subset(dat, metric=="productivity" & new.stan.metric!="proportion pairs with chicks" & stan.effect.size < 30 & species!="curlew" & species!="dunlin" & mowing!="applied" & grazing!="reduced")
mdat.prod.spec <- droplevels(mdat.prod.spec)
mdat.prod.spec$new.stan.metric <- as.factor(mdat.prod.spec$new.stan.metric)

mod.prod.spec <- lme(stan.effect.size ~ mowing + grazing + nest.protect + species + new.stan.metric, random= ~1|reference, data=mdat.prod.spec)
summary(mod.prod.spec)


#---------------------   Output model results  --------------------

setwd(outputwd)
sink(paste("1_effect size analysis.txt", sep=" "))

cat("\n########==========  1a) effect size of high-level management on abundance, abundance change, and nest survival \n", sep="\n")
print(lapply(mod.high, summary))

cat("\n########==========  1b) effect size of specific management on abundance, abundance change, and nest survival ==========########\n", sep="\n")
print(lapply(mod.spec, summary))

cat("\n########==========  1c) effect size of high-level management on pooled productivity measures (nest survival, chick survival, fledglings/pair) ==========########\n", sep="\n")
print(summary(mod.prod.high))

cat("\n########==========  1d) effect size of specific management on pooled productivity measures (nest survival, chick survival, fledglings/pair) ==========########\n", sep="\n")
print(summary(mod.prod.spec))
sink()

### Save individual interventions models
saveRDS(mod.high, file=paste(workspacewd, "models_1a.rds", sep="/"))
saveRDS(mod.spec, file=paste(workspacewd, "models_1b.rds", sep="/"))
saveRDS(mod.prod.high, file=paste(workspacewd, "models_1c.rds", sep="/"))
saveRDS(mod.prod.spec, file=paste(workspacewd, "models_1d.rds", sep="/"))

### Save dataset for models
saveRDS(moddat, file=paste(workspacewd, "model dataset_1a-b.rds", sep="/"))
saveRDS(mdat.prod.high, file=paste(workspacewd, "model dataset_1c.rds", sep="/"))
saveRDS(mdat.prod.spec, file=paste(workspacewd, "model dataset_1d.rds", sep="/"))

### Notes on model asumptions and outputs
# 3 very obvious outliers for abundance change analysis, all from the same study (#62 on snipe) - exclude these as very obvious outlier residuals
# 1 outlier for abundance analysis
# residual normality for abundance and nest survival is not amazing but not awful
# heterogeneity of residuals for nest survival is generally ok apart from extra large variance for fitted values ~ 2.05, but no clear pattern
# heterogeneity of residuals for abundance generally ok apart from extra large variance for fitted values ~ 3.8, but no obvious pattern


### Notes on model asumptions and outputs
# very obvious outliers for productivity analysis with effect size > 30 - exclude from analysis
# normality not bad, bit of a long right tail
# heterogeneity also not bad, bit of a large variance spread for some fitted values ~ 2.2


# see http://stats.stackexchange.com/questions/101020/how-to-interpret-multiple-factors-in-model-output-in-r for interpreting parameter estimates when there are multiple factors in a model
# The (intercept) indicates the value of the reference category. You have two categorical variables, so you have two reference categories. Your reference categories are LandUseHigh and an unspecified level of Type_LU (which I assume you know). So the value of the Estimate in the (intercept) row is the predicted mean for those study units in both of those categories when all the continuous covariates are equal to 00. Again, because you don't have an interaction term, the value of LandUseHigh when Type_LU is conserva, private, or state is the sum of the estimate for the intercept plus the estimate for the appropriate level of Type_LU.

# the reference level for all mgmt variables is 'none', so the model provides estimates of either applied or applied/reduced for each mgmt type relative to the baseline, which is all mgmt vars='none' and black-tailed godwit


#=================================  PLOT MODEL OUTPUTS  ===============================

setwd(outputwd)


### Read individual interventions models
mod.high <- readRDS(file=paste(workspacewd, "models_1a.rds", sep="/"))
mod.spec <- readRDS(file=paste(workspacewd, "models_1b.rds", sep="/"))
mod.prod.high <- readRDS(file=paste(workspacewd, "models_1c.rds", sep="/"))
mod.prod.spec <- readRDS(file=paste(workspacewd, "models_1d.rds", sep="/"))

### Read dataset
moddat <- readRDS(file=paste(workspacewd, "model dataset_1a-b.rds", sep="/"))
mdat.prod.high <- readRDS(file=paste(workspacewd, "model dataset_1c.rds", sep="/"))
mdat.prod.spec <- readRDS(file=paste(workspacewd, "model dataset_1d.rds", sep="/"))


# Ben Bolker: R-sig-mixed-models post 20 Feb 2010 https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q1/003336.html
# and GLMM wiki FAQs
# Ben Bolker: Do note (as commented in the FAQ) that this only accounts for the uncertainty of the fixed effects conditional on the estimates of the random-effect variances and BLUPs/conditional modes

# fm1 <- lme(distance ~ age*Sex, random = ~ 1 + age | Subject,
#            data = Orthodont)
# plot(Orthodont)
# 
# newdat <- expand.grid(age=c(8,10,12,14), Sex=c("Male","Female"))
# newdat$pred <- predict(fm1, newdat, level = 0)
# 
# Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
# predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
# newdat$SE <- sqrt(predvar) # for confidence intervals
# newdat$SE2 <- sqrt(predvar+fm1$sigma^2) # for prediction intervals, add the residual variance
# 
# library(ggplot2)
# pd <- position_dodge(width=0.4)
# ggplot(newdat,aes(x=age,y=pred,colour=Sex))+
#   geom_point(position=pd)+
#   geom_linerange(aes(ymin=pred-2*SE,ymax=pred+2*SE),
#                  position=pd)
# 
# ## prediction intervals
# ggplot(newdat,aes(x=age,y=pred,colour=Sex))+
#   geom_point(position=pd)+
#   geom_linerange(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),
#                  position=pd)


#----------   PLOTS OF 1a) EFFECTS OF HIGH-LEVEL INTERVENTIONS   ---------------

####     Abundance, abundance change, nest survival    ####
plotdat <- list()

for (i in 1:length(mod.high)) {
  
  # random effects not needed for new prediction
  
  plotmod <- mod.high[[i]] # model to plot results
  origdat <- moddat[[i]] # original dataset
  
  unique.mgmtvars <- unique(origdat[,mgmtvars]) # unique combination of mgmtvars appearing in the original dataset
  
  newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))
  
  newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))
  
  Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
  predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
  newdat$SE <- sqrt(predvar)
  newdat$lwr <- newdat$pred - (1.96*newdat$SE)
  newdat$upr <- newdat$pred + (1.96*newdat$SE)
  newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)
  
  fits <- newdat[,c("AE.level","reserve.desig","species","pred","SE","lwr","upr")]
  fits <- unique(fits)
  
  sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.level=fits$AE.level, reserve.desig=fits$reserve.desig), mean)
  
  # nest survival model doesn't include reserve.desig variable, so aggregate across reserve.desig as well as species
  if (metrics[i]=="nest survival (Mayfield)") {
    fits <- newdat[,c("AE.level","species","pred","SE","lwr","upr")]
    fits <- unique(fits)
    sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.level=fits$AE.level), mean)
  }
  
  plotdat[[i]] <- data.frame(sum.fits, metric=names(mod.high)[i])  
  
  
}

####  Pooled productivity measures effect size predictions  ####

plotmod <- mod.prod.high # model to plot results
origdat <- mdat.prod.high # original dataset

unique.mgmtvars <- unique(origdat[,mgmtvars]) # unique combination of mgmtvars appearing in the original dataset

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

# need to replicate across the different productivity measures as well
newdat <- data.frame(newdat[rep(seq_len(nrow(newdat)), times=3),], new.stan.metric=rep(levels(as.factor(origdat$new.stan.metric)), each=nrow(newdat)))

newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))

Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar)
newdat$lwr <- newdat$pred - (1.96*newdat$SE)
newdat$upr <- newdat$pred + (1.96*newdat$SE)
newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)

fits <- newdat[,c("AE.level","species","new.stan.metric","pred","SE","lwr","upr")]
fits <- unique(fits)

sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.level=fits$AE.level), mean)

plotdat[[4]] <- data.frame(sum.fits, metric="pooled productivity")



###-------- Output plot --------###
png("1a_high level intervention effect size.png", res=300, height=12, width=15, units="in", pointsize=20)

par(mfrow=c(2,2))

for (i in 1:length(metrics)) {
  
  par(mar=c(6,5,2,2))
  
  plotsub <- plotdat[[i]]
  
  x <- c(1:nrow(plotsub))
  
  if (metrics[i]=="abundance") {
    
    plotsub$pch <- c(0,1,2,15,16,17)
    
    plot(plotsub$pred~x, pch=plotsub$pch, cex=1.5, ylim=c(min(plotsub$lwr)*1.1,max(plotsub$upr)*1.1), xaxt="n", xlab="", ylab="Predicted effect size", las=1, bty="n")
    arrows(x, plotsub$pred, x, plotsub$lwr, angle=90, length=0.05)
    arrows(x, plotsub$pred, x, plotsub$upr, angle=90, length=0.05)
    abline(h=0, lty=3, lwd=2)
    # axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
    axis(1, x, labels=rep("",nrow(plotsub)), tick=TRUE)
    text(x, par("usr")[3]*1.1, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=0.8)
    text(c(2,5), par("usr")[3]*1.6, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/ \n designation", "nature reserve/ \n designation"), font=2)
    title(metrics[i])
  }
  
  if (metrics[i]=="multiplicative yearly slope") {
    
    plotsub$pch <- c(0,1,2,15,16,17)
    
    plot(plotsub$pred~x, pch=plotsub$pch, cex=1.5, ylim=c(min(plotsub$lwr)*1.5,max(plotsub$upr)*1.2), xaxt="n", xlab="", ylab="Predicted effect size", las=1, bty="n")
    arrows(x, plotsub$pred, x, plotsub$lwr, angle=90, length=0.05)
    arrows(x, plotsub$pred, x, plotsub$upr, angle=90, length=0.05)
    abline(h=0, lty=3, lwd=2)
    #     axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE)
    axis(1, x, labels=rep("",nrow(plotsub)), tick=TRUE)
    text(x, par("usr")[3]*1.15, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=0.8)
    text(c(2,5), par("usr")[3]*1.8, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/ \n designation", "nature reserve/ \n designation"), font=2)
    title(metrics[i])
  }
  
  if (metrics[i]=="nest survival (Mayfield)") {
    
    x <- c(2:3)
    
    plotsub$pch <- c(0,2)
    
    plot(plotsub$pred~x, pch=plotsub$pch, cex=1.5, ylim=c(min(plotsub$lwr)*1.3,max(plotsub$upr)*1.2), xlim=c(1.5,3.5), xaxt="n", xlab="", ylab="Predicted effect size", las=1, bty="n")
    arrows(x, plotsub$pred, x, plotsub$lwr, angle=90, length=0.05)
    arrows(x, plotsub$pred, x, plotsub$upr, angle=90, length=0.05)
    abline(h=0, lty=3, lwd=2)
    #     axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE)
    axis(1, x, labels=rep("",nrow(plotsub)), tick=TRUE)
    text(x, par("usr")[3]*1.15, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","higher-level \n AES"), cex=0.8)
    title(metrics[i])
  }
}

### Plot pooled productivity ###

par(mar=c(6,5,2,2))

plotsub <- plotdat[[4]]

x <- c(1:nrow(plotsub))
x <- c(2:3)

plotsub$pch <- c(0,2)

plot(plotsub$pred~x, pch=plotsub$pch, cex=1.5, ylim=c(min(plotsub$lwr)*1.3,max(plotsub$upr)*1.2), xlim=c(1.5,3.5), xaxt="n", xlab="", ylab="Predicted effect size", las=1, bty="n")
arrows(x, plotsub$pred, x, plotsub$lwr, angle=90, length=0.05)
arrows(x, plotsub$pred, x, plotsub$upr, angle=90, length=0.05)
abline(h=0, lty=3, lwd=2)
#     axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE)
axis(1, x, labels=rep("",nrow(plotsub)), tick=TRUE)
text(x, par("usr")[3]*1.15, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","higher-level \n AES"), cex=0.8)
title("pooled productivity metrics")


dev.off()



#----------   PLOTS OF 1b) EFFECTS OF SPECIFIC INTERVENTIONS   ---------------

####     Abundance only   ####

plotdat.spec <- list()

for (i in 1:1) {
  
  # random effects not needed for new prediction
  
  plotmod <- mod.spec[[i]] # model to plot results
  origdat <- moddat[[i]] # original dataset
  
  unique.mgmtvars <- unique(origdat[,mgmtvars]) # unique combination of mgmtvars appearing in the original dataset
  
  newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))
  
  newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))
  
  Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
  predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
  newdat$SE <- sqrt(predvar)
  newdat$lwr <- newdat$pred - (1.96*newdat$SE)
  newdat$upr <- newdat$pred + (1.96*newdat$SE)
  newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)
  
  fits <- newdat[,c("mowing","grazing","nest.protect","water","species","pred","SE","lwr","upr")]
  fits <- unique(fits)
  
  sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(mowing=fits$mowing, grazing=fits$grazing, nest.protect=fits$nest.protect, water=fits$water), mean)
  
  # ordered from worst to best combinations of interventions
  sum.fits <- sum.fits[order(sum.fits$pred),]
  
  # sum.fits <- sum.fits[order(sum.fits$mowing, sum.fits$grazing, sum.fits$nest.protect, sum.fits$water),]
  
  plotdat.spec[[i]] <- data.frame(sum.fits, metric=names(mod.spec)[i])  
  
  
}



###-------- Output plot --------###
png("1b_specific intervention effect size.png", res=300, height=12, width=18, units="in", pointsize=20)

for (i in 1:length(plotdat.spec)) {
  
  par(mar=c(6,5,2,2))
  
  plotsub <- plotdat.spec[[i]]
  
  x <- c(1:nrow(plotsub))
  
  plot(plotsub$pred~x, pch=16, cex=1.5, ylim=c(min(plotsub$lwr)*1.1,max(plotsub$upr)*1.1), xaxt="n", xlab="", ylab="Predicted effect size", las=1, bty="n")
  arrows(x, plotsub$pred, x, plotsub$lwr, angle=90, length=0.05)
  arrows(x, plotsub$pred, x, plotsub$upr, angle=90, length=0.05)
  abline(h=0, lty=3, lwd=2)
  # axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
  axis(1, x, labels=rep("",nrow(plotsub)), tick=TRUE)
  text(x, par("usr")[3]*1.1, srt = 0, pos=1, xpd = TRUE, cex=0.8,
       labels=c("mowing reduced +\ngrazing reduced +\nwater applied",
                "grazing applied",
                "mowing reduced",
                "mowing reduced +\nnest protection",
                "no interventions",
                "nest protection",
                "mowing reduced +\nwater applied",
                "water applied"))
  title("abundance")
}

dev.off()

