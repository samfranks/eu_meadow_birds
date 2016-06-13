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

#=================================  LOAD FUNCTIONS =================================

### Ben Bolker's function for calculating CIs on predictions from a merMod object and plotting the results from his RPubs GLMM worked examples
# http://rpubs.com/bbolker/glmmchapter
# by specifying re.form=NA we're saying that we want the population-level prediction, i.e. setting the random effects to zero and getting a prediction for an average (or unknown) group
# Computing confidence intervals on the predicted values is relatively easy if we're willing to completely ignore the random effects, and the uncertainty of the random effects
# this easy method produces similar width CIs to using the bootMer function in lme4, perhaps slightly wider CIs in some cases

# can change to alpha=0.16, approximately equal to 84% CIs
easyPredCI.lme <- function(model,newdata,alpha=alphalevel) {
  
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

mod <- list()
moddat <- list()

for (i in 1:length(metrics)) {
  
  print(metrics[i])
  mdat <- subset(dat, new.stan.metric==metrics[i])
  
  #   # check for sample sizes of different management treatment groups, remove small ones
  #   checksample <- table(mdat$mgmt)
  #   removelevels <- names(checksample)[checksample < 5]
  
  # DFs for predator control (applied and reduced), water reduced are very small because they are each only tested in a single study so the model can't estimate the effect of those interventions, separately from the effect of the random effect study (checked this with Ali)
  # not gaining anything useful from this so remove these interventions
  if (metrics[i]=="abundance") {
    mdat <- subset(mdat, predator.control!="applied" & predator.control!="reduced" & fertpest!="reduced" & water!="reduced")
    mdat <- subset(mdat, stan.effect.size < 16) # remove large outlier
  }
  if (metrics[i]=="multiplicative yearly slope") mdat <- subset(mdat, stan.effect.size > -0.4) # remove large outlier
  if (metrics[i]=="nest survival (Mayfield)") mdat <- subset(mdat, species!="curlew" & species!="dunlin") # control for too few studies assessing curlew or dunlin
  
  # removelevels <- c("predator.control applied", "predator.control reduced", "water reduced")
  
  # if (length(removelevels>0)) mdat <- subset(mdat, !grepl(paste(removelevels,collapse="|"), mdat$mgmt)) else mdat <- mdat
  mdat <- droplevels(mdat)
  moddat[[i]] <- mdat
  
  print(nrow(mdat))

  if (metrics[i]=="abundance") mod[[i]] <- lme(stan.effect.size ~ AE.level + reserve.desig + mowing + grazing +  nest.protect + water + species, random = ~1|reference, data=mdat) # no predator control, fertpest records that aren't none in abundance model
  print(summary(mod[[i]]))
  
  # see http://stats.stackexchange.com/questions/101020/how-to-interpret-multiple-factors-in-model-output-in-r for interpreting parameter estimates when there are multiple factors in a model
  # The (intercept) indicates the value of the reference category. You have two categorical variables, so you have two reference categories. Your reference categories are LandUseHigh and an unspecified level of Type_LU (which I assume you know). So the value of the Estimate in the (intercept) row is the predicted mean for those study units in both of those categories when all the continuous covariates are equal to 00. Again, because you don't have an interaction term, the value of LandUseHigh when Type_LU is conserva, private, or state is the sum of the estimate for the intercept plus the estimate for the appropriate level of Type_LU.
  
  # the reference level for all mgmt variables is 'none', so the model provides estimates of either applied or applied/reduced for each mgmt type relative to the baseline, which is all mgmt vars='none' and black-tailed godwit
  
}


############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################


#------------  Reshape data from wide to long, bringing interventions under one column -----------------

dat.wide <- dat

# remove levels of 'none' - we only want to assess effect size for use
dat <- gather(dat.wide, key=intervention, value=mgmtlevel, AE:water)
dat <- subset(dat, mgmtlevel!="none")

# combine interventions and their levels into a single variable
dat$mgmt <- paste(dat$intervention, dat$mgmtlevel)

#=================================  ANALYSIS  ===============================

# can combine standardised abundance metrics - effect sizes comparable
# combine chick survival standardised metrics - effect sizes comparable

dat$new.stan.metric <- ifelse(grepl("chick survival", dat$stan.metric), "chick survival", ifelse(grepl("number of", dat$stan.metric), "abundance", as.character(dat$stan.metric)))

#------------------------------ Test effect of nuisance variables on success for the full dataset -------------------------

# metrics to test - other metrics have sample sizes which are too small
metrics <- c("abundance","multiplicative yearly slope","nest survival (Mayfield)")

###----  Nuisance variables, all metrics pooled ----###

nui.dat <- unique(dat[,c("reference","study.length","sample.size","analysis2","lit.type","stan.effect.size")])

# use unique dataset only to test effect sizes (replication because of multiples of same study, same effect sizes, but with different interventions used)
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


mod <- list()
moddat <- list()

for (i in 1:length(metrics)) {
  
  print(metrics[i])
  mdat <- subset(dat, new.stan.metric==metrics[i])
  
#   # check for sample sizes of different management treatment groups, remove small ones
#   checksample <- table(mdat$mgmt)
#   removelevels <- names(checksample)[checksample < 5]
  
  # DFs for predator control (applied and reduced), water reduced are very small because they are each only tested in a single study so the model can't estimate the effect of those interventions, separately from the effect of the random effect study (checked this with Ali)
  # not gaining anything useful from this so remove these interventions
  if (metrics[i]=="abundance") {
    mdat <- subset(mdat, mgmt!="predator.control applied" & mgmt!="predator.control reduced" & mgmt!="water reduced")
    mdat <- subset(mdat, stan.effect.size < 16) # remove large outlier
  }
  if (metrics[i]=="multiplicative yearly slope") mdat <- subset(mdat, stan.effect.size > -0.4) # remove large outlier
  if (metrics[i]=="nest survival (Mayfield)") mdat <- subset(mdat, species!="curlew" & species!="dunlin") # control for too few studies assessing curlew or dunlin
  
    # removelevels <- c("predator.control applied", "predator.control reduced", "water reduced")
  
  # if (length(removelevels>0)) mdat <- subset(mdat, !grepl(paste(removelevels,collapse="|"), mdat$mgmt)) else mdat <- mdat
  mdat <- droplevels(mdat)
  mdat$mgmt <- as.factor(mdat$mgmt)
  moddat[[i]] <- mdat
  
  print(nrow(mdat))
  print(table(mdat$mgmt))
  print(table(mdat$mgmt,mdat$reference))
  
  mod[[i]] <- lme(stan.effect.size ~ mgmt + species, random = ~1|reference, data=mdat)
  print(summary(mod[[i]]))
  
}

### Notes on model asumptions and outputs
# 3 very obvious outliers for abundance change analysis, all from the same study (#62 on snipe) - exclude these as very obvious outlier residuals
# 1 outlier for abundance analysis
# residual normality for abundance and nest survival is not amazing but not awful
# heterogeneity of residuals for nest survival is generally ok apart from extra large variance for fitted values ~ 2.05, but no clear pattern
# heterogeneity of residuals for abundance generally ok apart from extra large variance for fitted values ~ 3.8, but no obvious pattern

####    Effect size on nest/chick survival measures pooled    ####
# run effect size for all productivity measures pooled (apart from proportion pairs with chicks since it's a very different metric), but with factor variable controlling for type of productivity measure if interventions have different effects on nests vs chicks etc
mdat.prod <- subset(dat, metric=="productivity" & new.stan.metric!="proportion pairs with chicks" & stan.effect.size < 30)
mod.prod <- lme(stan.effect.size ~ mgmt + new.stan.metric, random= ~1|reference, data=mdat.prod)
summary(mod.prod)

### Notes on model asumptions and outputs
# very obvious outliers for productivity analysis with effect size > 30 - exclude from analysis
# normality not bad, bit of a long right tail
# heterogeneity also not bad, bit of a large variance spread for some fitted values ~ 2.2

#---------------------   Test for all contrasts in effect sizes between intervention types  --------------------

# for making pairwise comparisons between the different treatment groups, use glht() in multcomp package - nice description of how given here: http://mindingthebrain.blogspot.co.uk/2013/04/multiple-pairwise-comparisons-for.html

# currently, parameter tables give the difference for each management's estimate with the reference level (AE applied)
# we want tests comparing all of the effects of the different management types against each other - which ones are significantly better than others (not just which are sig different than the reference level)?

### Abundance ###

contrast.matrix <- rbind(
  `mgmtAE.level basic vs. mgmtAE.level higher` = c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  `Time:Diet1 vs. Time:Diet3` = c(0, 0, 0, 0, 0, 0, 1, 0),
  `Time:Diet1 vs. Time:Diet4` = c(0, 0, 0, 0, 0, 0, 0, 1),
  `Time:Diet2 vs. Time:Diet3` = c(0, 0, 0, 0, 0, -1, 1, 0),
  `Time:Diet2 vs. Time:Diet4` = c(0, 0, 0, 0, 0, -1, 0, 1),
  `Time:Diet3 vs. Time:Diet4` = c(0, 0, 0, 0, 0, 0, -1, 1))


amod <- lm(breaks ~ tension, data = warpbreaks)
### set up all-pair comparisons for factor `tension'
### using a symbolic description (`type' argument
### to `contrMat()')
com1 <- glht(amod, linfct = mcp(tension = "Tukey"))
### alternatively, describe differences symbolically
com2 <- glht(amod, linfct = mcp(tension = c("M - L = 0",
                                    "H - L = 0",
                                    "H - M = 0")))


x <- glht(mod[[1]], linfct = mcp(mgmt = "Tukey"))

#=================================  PLOT MODEL OUTPUTS  ===============================

# Ben Bolker: R-sig-mixed-models post 20 Feb 2010 https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q1/003336.html
# and GLMM wiki FAQs
# Ben Bolker: Do note (as commented in the FAQ) that this only accounts for the uncertainty of the fixed effects conditional on the estimates of the random-effect variances and BLUPs/conditional modes

fm1 <- lme(distance ~ age*Sex, random = ~ 1 + age | Subject,
           data = Orthodont)
plot(Orthodont)

newdat <- expand.grid(age=c(8,10,12,14), Sex=c("Male","Female"))
newdat$pred <- predict(fm1, newdat, level = 0)

Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) # for confidence intervals
newdat$SE2 <- sqrt(predvar+fm1$sigma^2) # for prediction intervals, add the residual variance

library(ggplot2)
pd <- position_dodge(width=0.4)
ggplot(newdat,aes(x=age,y=pred,colour=Sex))+
  geom_point(position=pd)+
  geom_linerange(aes(ymin=pred-2*SE,ymax=pred+2*SE),
                 position=pd)

## prediction intervals
ggplot(newdat,aes(x=age,y=pred,colour=Sex))+
  geom_point(position=pd)+
  geom_linerange(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),
                 position=pd)


plotdat <- list()

for (i in 1:length(mod)) {
  
  # random effects not needed for new prediction
  
  plotmod <- mod[[i]] # model to plot results
  origdat <- moddat[[i]] # original dataset
  
  newdat <- with(origdat, data.frame(stan.effect.size, mgmt, species)) # create new dataset with species 
  
  newdat 
  
  newdat$pred <- predict(plotmod, level=0, newdat)
  
  Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
  predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
  newdat$SE <- sqrt(predvar)
  newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)
  
  fits <- data.frame(pred=pred$fit, se.fit=pred$se.fit, lci=(pred$fit - (1.96*pred$se.fit)), uci=(pred$fit + (1.96*pred$se.fit)))
  
  
  pred.CI <- easyPredCI(mod[[i]], moddat[[i]])
  
  fits <- data.frame(pred,pred.CI,lit.type=moddat[[i]][,"lit.type"],mgmtvar=paste(mgmtvars[i], moddat[[i]][,mgmtvars[i]]))
  unique.fits <- unique(fits)
  
  plotdat[[i]] <- aggregate(unique.fits[,c("pred","lwr","upr")], by=list(mgmtvar=unique.fits$mgmtvar), mean)
  
}

plotfinal <- do.call(rbind, plotdat)

###-------- Output plot --------###
png("0a_overall model results.png", res=300, height=12, width=20, units="in", pointsize=20)
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
write.csv(plotfinal[,c("pred","lwr","upr","mgmtvar")], "0a_overall probabilities and CIs.csv", row.names=FALSE)




easyPredCI <- function(model,newdata,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  
  model <- mod[[i]]
  newdata <- moddat[[i]]
  
  pred0 <- predict(model,newdata=newdata)
  
  pred <- predict(model, level=0, newdata=newdata)
  
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
cpred1.CI <- easyPredCI(cmod_lme4_L,pframe)
