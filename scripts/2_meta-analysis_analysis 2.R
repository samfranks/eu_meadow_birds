#################################################################
#
#     Step 2: EU meadow birds meta-analysis - models evaluating success
#
#################################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================


#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme")

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

source(paste(scriptswd, "2_meta-analysis_data preparation.R", sep="/"))



#=================================  ANALYSIS  2 ===============================


#-----------------   HIGH-LEVEL INTERVENTION SUCCESS  ------------------



####---- Method 1: Run models using only those records which evaluate either AES and/or nature reserves ----####



# this tells us the probability of success IF THESE INTERVENTIONS WERE ATTEMPTED (subtly different than below Method 2, which tells us the probability of success if they were used vs not used)
# structured the data to combine AE and reserves into a single variable because otherwise lack of the none/none category produces rank deficiency in the model and it drops a coefficient

# > table(mdat$AE,mdat$reserve.desig)
# 
# none applied
# none       0      87
# applied  190      72

mdat <- subset(dat, high.int.used==1)
mdat <- subset(mdat, species!="ruff" & species!="dunlin")
mdat <- droplevels(mdat)

newlevels <- data.frame(AE=c("applied","none","applied"), reserve.desig=c("none","applied","applied"), AE.reserve=c("AE-no reserve","no AE-reserve","AE-reserve"))

mdat <- merge(mdat, newlevels, by=c("AE","reserve.desig"))

mdat$AE.reserve <- relevel(mdat$AE.reserve, ref="no AE-reserve")
table(mdat$AE.reserve)

newlevels2 <- data.frame(AE.level=c("basic","higher","none","basic","higher"), reserve.desig=c("none","none","applied","applied","applied"), AE.reserve2=c("AE.basic-no reserve","AE.higher-no reserve","no AE-reserve","AE.basic-reserve","AE.higher-reserve"))

mdat <- merge(mdat, newlevels2, by=c("AE.level","reserve.desig"))

mdat$AE.reserve2 <- relevel(mdat$AE.reserve2, ref="no AE-reserve")
table(mdat$AE.reserve2)

###--- Model ---###

m.high <- glmer(success ~ AE.reserve + species + (1|reference), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.high)

library(multcomp)

setwd(outputwd)
sink(paste("model output_2a_method 1.txt", sep=" "))
cat("\n########==========  Success of higher-level interventions combined ==========########\n", sep="\n")
print(summary(m.high))

cat("\n###----  Tukey contrasts  ---###\n", sep="\n")
print(summary(glht(m.high, linfct = mcp(AE.reserve="Tukey"))))

sink()

### Save individual interventions models
saveRDS(m.high, file=paste(workspacewd, "models_2a_method 1.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2a_method 1.rds", sep="/"))




####---- Method 2: Run AE*nature reserve models where none/none categories are included ----####



# this tells us about probability of success for cases where these interventions are used vs not used
# we want to know whether there is a greater prob of success when AE is applied together with reserves
# i.e. is success of AE applied-reserve applied greater than AE applied-reserve none or AE none-reserve applied? (we don't really care about AE none-reserve none)

mdat <- subset(dat, species!="ruff")
mdat <- droplevels(mdat)

m.high <- glmer(success ~ AE*reserve.desig + species + (1|reference), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.high)


setwd(outputwd)
sink(paste("model output_2a_method 2.txt", sep=" "))
cat("\n########==========  Success of higher-level interventions combined ==========########\n", sep="\n")
print(summary(m.high))

cat("\n###----  Significance of interaction term  ---###\n", sep="\n")
print(drop1(m.high, scope = ~AE:reserve.desig, test="Chisq"))

sink()

### Save individual interventions models
saveRDS(m.high, file=paste(workspacewd, "models_2a_method 2.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2a_method 2.rds", sep="/"))



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
m.spec <- glmer(success ~ mowing + grazing + fertpest + nest.protect + predator.control + water + species + (1|reference), data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
summary(m.spec)

setwd(outputwd)
sink(paste("model output_2b.txt", sep=" "))

cat("\n########==========  Success of specific interventions combined ==========########\n", sep="\n")
print(summary(m.spec))

cat("\n###----  Significance of each intervention (in the presence of others)  ---###\n", sep="\n")
print(drop1(m.spec, scope = ~mowing, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~grazing, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~fertpest, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~nest.protect, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~predator.control, test="Chisq"))
cat("\n")
print(drop1(m.spec, scope = ~water, test="Chisq"))

sink()

### Save individual interventions models
saveRDS(m.spec, file=paste(workspacewd, "models_2b.rds", sep="/"))

### Save dataset for 0b models
saveRDS(mdat, file=paste(workspacewd, "model dataset_2b.rds", sep="/"))






#=================================  PLOT OUTPUTS  ===============================

setwd(outputwd)


#-------------  HIGHER LEVEL INTERVENTIONS - METHOD 1 -------------

###----  Produce plotting dataset predictions ----###

# Model
plotmod <- readRDS(file=paste(workspacewd, "models_2a_method 1.rds", sep="/")) # model to plot results
# Model dataset
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

###---- Output plot ----###

png("2a_high level combination intervention success_method 1.png", res=300, height=12, width=12, units="in", pointsize=20)

par(mar=c(6,6.5,2,3))

x <- c(1:nrow(plotdat))

plotdat$pch <- c(16,16,16)
# plotdat$pch <- c(1,2,15,16,17)

plot(plotdat$pred~x, pch=plotdat$pch, cex=1.5, ylim=c(0,1), xlim=c(0.8,3.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3]*2, srt = 0, pos=1, xpd = TRUE, labels=c("AES only","nature reserve only", "AES + nature reserve"), cex=1)
# text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=4.5)

dev.off()



#-------------  HIGHER LEVEL INTERVENTIONS - METHOD 2 -------------

###----  Produce plotting dataset predictions ----###

# Model
plotmod <- readRDS(file=paste(workspacewd, "models_2a_method 2.rds", sep="/")) # model to plot results
# Model dataset
origdat <- readRDS(file=paste(workspacewd, "model dataset_2a_method 2.rds", sep="/")) # original dataset

unique.mgmtvars <- unique(origdat[,c("AE","reserve.desig")])

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

pred <- predict(plotmod, newdat, type="response", re.form=NA)
pred.CI <- easyPredCI(plotmod, newdat)
fits <- data.frame(newdat, pred, pred.CI)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(AE=fits$AE, reserve.desig=fits$reserve.desig), mean)

plotdat <- sum.fits
plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE),]
plotdat <- plotdat[-1,]

###---- Output plot ----###

setwd(outputwd)
png("2a_high level combination intervention success_method 2.png", res=300, height=12, width=12, units="in", pointsize=20)

par(mar=c(6,6.5,2,3))

x <- c(1:nrow(plotdat))

plotdat$pch <- c(16,16,16)
# plotdat$pch <- c(1,2,15,16,17)

plot(plotdat$pred~x, pch=plotdat$pch, cex=1.5, ylim=c(0,1), xlim=c(0.8,3.2), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3]*2, srt = 0, pos=1, xpd = TRUE, labels=c("AES only","nature reserve only", "AES + nature reserve"), cex=1)
# text(x, par("usr")[3]*2, srt = 0, pos=1, xpd = TRUE, labels=c("no AES or reserves","AES only","nature reserve only", "AES + nature reserve"), cex=1)
# text(x, par("usr")[3]*1.2, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(x, par("usr")[3]*1.5, srt = 0, pos=1, xpd = TRUE, labels=c("no AES","basic-level \n AES","higher-level \n AES"), cex=1)
# text(c(2,5), par("usr")[3]*4, srt = 0, pos=1, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=4.5)

dev.off()


#-------------  SPECIFIC INTERVENTIONS - main intervention effects -------------

###----  Produce plotting dataset predictions ----###

### Read model
plotmod <- readRDS(file=paste(workspacewd, "models_2b.rds", sep="/"))
# Read model dataset
origdat <- readRDS(file=paste(workspacewd, "model dataset_2b.rds", sep="/")) # original dataset

unique.mgmtvars <- unique(origdat[,mgmtvars[4:9]]) # unique combination of mgmtvars appearing in the original dataset

# unique.mgmtvars <- data.frame(mowing=c("applied","reduced","none","none","none","none","none","none","none"), grazing=c("none","none","applied","reduced","none","none","none","none","none"), fertpest=c("none","none","none","none","applied","reduced","none","none","none"), nest.protect=c("none","none","none","none","none","none","applied","none","none"), predator.control=c("none","none","none","none","none","none","none","applied","none"), water=c("none","none","none","none","none","none","none","none","water"))

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

pred <- predict(plotmod, newdat, type="response", re.form=NA)
pred.CI <- easyPredCI(plotmod, newdat)
fits <- data.frame(newdat, pred, pred.CI)

# produce mean population level prediction for all intervention combinations across species
sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(mowing=fits$mowing, grazing=fits$grazing, fertpest=fits$fertpest, nest.protect=fits$nest.protect, predator.control=fits$predator.control, water=fits$water), mean)
plotdat <- sum.fits
plotdat <- plotdat[order(plotdat$mowing, plotdat$grazing, plotdat$fertpest, plotdat$nest.protect, plotdat$predator.control, plotdat$water),]

# count up number of interventions applied at once and add sum to plotdat
interventions.only <- plotdat[,1:6]
num.interventions <- ifelse(interventions.only=="none",0,1)
num.interventions.sum <- apply(num.interventions, 1, sum)

plotdat <- data.frame(plotdat, num.interventions.sum)

# add variable showing whether success is 'significant' (i.e. lwr does not overlap 0.05)
plotdat$sig <- ifelse(plotdat$lwr > 0.05, "Y", "N")

# order by 'significance' and magnitude of success
plotdat <- plotdat[order(plotdat$sig, plotdat$pred),]

saveRDS(plotdat, file=paste(workspacewd, "2b_all combinations plotting dataset.rds", sep="/"))


# produce overall effect of intervention level, averaged across predicted values of other intervention levels, and also averaged across all species
avg.fits <- list()
for (i in 1:6) {
  temp <- aggregate(sum.fits[,c("pred","lwr","upr")], by=list(level=sum.fits[,i]), mean)
  temp <- data.frame(temp, intervention=names(sum.fits[i]))
  avg.fits[[i]] <- temp[-which(temp[,1]=="none"),]
}
plotdat <- do.call(rbind, avg.fits)

saveRDS(plotdat, file=paste(workspacewd, "2b_average intervention plotting dataset.rds", sep="/"))
write.csv(plotdat, "2b_average intervention probabilities and CIs.csv", row.names=FALSE)



###---- Output plot - combination effects ----###

setwd(outputwd)
png("2b_specific combination intervention success.png", res=300, height=15, width=30, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(1,8,1,2))

plotdat <- readRDS(file=paste(workspacewd, "2b_all combinations plotting dataset.rds", sep="/"))

x <- c(1:nrow(plotdat))-0.5 # gives an x axis to plot against

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(0,1), xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(1,26))
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05) # add error bars, lwr and upr to each prediction
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2) # add a 'significance' line (what is the threshold for 'success'?)
abline(v=max(which(plotdat$sig=="N")), lty=3, lwd=2) # add a line dividing 'successful' vs 'unsuccessful' intervention combos
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=0, xpd=TRUE)
title(ylab="Predicted probability of success \n (significant positive impact) ", cex.lab=1.5, font=2, line=3)


### Create a 'table' of intervention combinations to display below the plot showing the predicted success ###

y <- length(mgmtvars[4:9]):1 # how many interventions there are (will be labels down the y-axis of the table starting at the top and working down)

# create the table of coordinates for the table (centres of the grid cells for the table), which is 28 across (the number of different intervention combinations) x 6 down (the number of types of interventions)
x <- c(1:nrow(plotdat))
new.x <- rep(x, each=max(y))
new.y <- rep(y, times=max(x))
tab <- data.frame(x=new.x,y=new.y)
tab <- tab[order(tab$y, decreasing=TRUE),]
labs <- plotdat[,1:6]

labs.long <- gather(labs, intervention, level, mowing:water) # convert the interventions from wide to long format using tidyr

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
axis(2, y+0.5, labels=c("mowing","grazing","fertiliser/\npesticides", "nest protection","predator control", "water"), las=1, cex.axis=1, font=2,tick=FALSE) # draw the labels for the rows, from top to bottom
text(tab$x-0.5, tab$y+0.5, labels=ifelse(tab.filled$level=="applied", "\U2191", ifelse(tab.filled$level=="reduced", "\U2193", tab.filled$level)), cex=2) # fill in the values of the grid cells, but if an intervention was applied then use a unicode 'up' arrow, and if it was reduced than use a down arrow
text(tab$x-0.5, 0.5, labels=total.interventions$sum, cex=1) # add the total number of interventions to the bottom 'row' of the table


dev.off()



###---- Output plot - averaged effects ----###

setwd(outputwd)
png("2b_specific combination intervention success_average.png", res=300, height=12, width=26, units="in", pointsize=20)

plotdat <- readRDS(file=paste(workspacewd, "2b_average intervention plotting dataset.rds", sep="/"))

par(mar=c(6,6,1,2))

x <- c(1:nrow(plotdat)) # gives an x axis to plot against

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(0,1), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05) # add error bars, lwr and upr to each prediction
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
axis(1, x, labels=rep("",nrow(plotdat)), tick=TRUE)
text(x, par("usr")[3]*2, srt = 0, pos=1, xpd = TRUE, labels=c("mowing \napplied","mowing \nreduced","grazing \napplied","grazing \nreduced","fertiliser/pesticides \napplied","fertiliser/pesticides \nreduced","nest protection \napplied","predator control \napplied","water \napplied", "water \nreduced"), cex=1)
title(xlab="Intervention (controlling for other interventions used)", cex.lab=1.5, font=2, line=4, xpd=TRUE)
title(ylab="Predicted probability of success \n (significant positive impact) ", cex.lab=1.5, font=2, line=3)

dev.off()




#-------------  SPECIFIC INTERVENTIONS - combination intervention effects -------------#

###----  Produce plotting dataset predictions ----###

### Read model
plotmod <- readRDS(file=paste(workspacewd, "models_2b.rds", sep="/"))
# Read model dataset
origdat <- readRDS(file=paste(workspacewd, "model dataset_2b.rds", sep="/")) # original dataset

unique.mgmtvars <- unique(origdat[,mgmtvars[4:9]]) # unique combination of mgmtvars appearing in the original dataset

newdat <- data.frame(unique.mgmtvars[rep(seq_len(nrow(unique.mgmtvars)), times=length(levels(origdat$species))),], species=rep(levels(origdat$species), each=nrow(unique.mgmtvars)))

pred <- predict(plotmod, newdat, type="response", re.form=NA)
pred.CI <- easyPredCI(plotmod, newdat)
fits <- data.frame(newdat, pred, pred.CI)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","lwr","upr")], by=list(mowing=fits$mowing, grazing=fits$grazing, fertpest=fits$fertpest, nest.protect=fits$nest.protect, predator.control=fits$predator.control, water=fits$water), mean)

plotdat <- sum.fits
# plotdat <- plotdat[order(plotdat$mowing, plotdat$grazing, plotdat$fertpest, plotdat$nest.protect, plotdat$predator.control, plotdat$water),]

# count up number of interventions applied at once and add sum to plotdat
interventions.only <- plotdat[,1:6]
num.interventions <- ifelse(interventions.only=="none",0,1)
num.interventions.sum <- apply(num.interventions, 1, sum)

plotdat <- data.frame(plotdat, num.interventions.sum)

# add variable showing whether success is 'significant' (i.e. lwr does not overlap 0.05)
plotdat$sig <- ifelse(plotdat$lwr > 0.05, "Y", "N")

# order by 'significance' and magnitude of success
plotdat <- plotdat[order(plotdat$sig, plotdat$pred),]



###---- Output plot ----###

setwd(outputwd)
png("2b_specific combination intervention success.png", res=300, height=15, width=30, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(1,8,1,2))

x <- c(1:nrow(plotdat))-0.5 # gives an x axis to plot against

plot(plotdat$pred~x, pch=16, cex=1.5, ylim=c(-0.3,1), xaxt="n", xlab="", ylab="", las=1, bty="n", xlim=c(1,26))
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05) # add error bars, lwr and upr to each prediction
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2) # add a 'significance' line (what is the threshold for 'success'?)
abline(v=max(which(plotdat$sig=="N")), lty=3, lwd=2) # add a line dividing 'successful' vs 'unsuccessful' intervention combos
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=0, xpd=TRUE)
title(ylab="Predicted probability of success \n (significant positive impact) ", cex.lab=1.5, font=2, line=3)


### Create a 'table' of intervention combinations to display below the plot showing the predicted success ###

y <- length(mgmtvars[4:9]):1 # how many interventions there are (will be labels down the y-axis of the table starting at the top and working down)

# create the table of coordinates for the table (centres of the grid cells for the table), which is 28 across (the number of different intervention combinations) x 6 down (the number of types of interventions)
x <- c(1:nrow(plotdat))
new.x <- rep(x, each=max(y))
new.y <- rep(y, times=max(x))
tab <- data.frame(x=new.x,y=new.y)
tab <- tab[order(tab$y, decreasing=TRUE),]
labs <- plotdat[,1:6]

labs.long <- gather(labs, intervention, level, mowing:water) # convert the interventions from wide to long format using tidyr

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
axis(2, y+0.5, labels=c("mowing","grazing","fertiliser/\npesticides", "nest protection","predator control", "water"), las=1, cex.axis=1, font=2,tick=FALSE) # draw the labels for the rows, from top to bottom
text(tab$x-0.5, tab$y+0.5, labels=ifelse(tab.filled$level=="applied", "\U2191", ifelse(tab.filled$level=="reduced", "\U2193", tab.filled$level)), cex=2) # fill in the values of the grid cells, but if an intervention was applied then use a unicode 'up' arrow, and if it was reduced than use a down arrow
text(tab$x-0.5, 0.5, labels=total.interventions$sum, cex=1) # add the total number of interventions to the bottom 'row' of the table


dev.off()



#-------------  HIGHER LEVEL INTERVENTIONS - METRIC  -------------#

###----  Produce plotting dataset predictions ----###

plotmod <- m.high.metric # model to plot results
origdat <- dat.high.metric # original dataset

unique.mgmtvars <- unique(origdat$AE.reserve)

newdat <- data.frame(AE.reserve=rep(unique.mgmtvars, times=length(levels(origdat$species))), species=rep(levels(origdat$species), each=length(unique.mgmtvars)))
newdat <- data.frame(newdat[rep(seq_len(nrow(newdat)), times=length(levels(as.factor(origdat$new.metric)))),], new.metric=rep(levels(as.factor(origdat$new.metric)), each=nrow(newdat)))

newdat$pred <- as.numeric(predict(plotmod, level=0, newdat))

Designmat <- model.matrix(eval(eval(plotmod$call$fixed)[-2]), newdat[-which(names(newdat) %in% "pred")])
predvar <- diag(Designmat %*% plotmod$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar)
newdat$lwr <- newdat$pred - (1.96*newdat$SE)
newdat$upr <- newdat$pred + (1.96*newdat$SE)
newdat$SE2 <- sqrt(predvar + plotmod$sigma^2)

fits <- newdat[,c("AE.reserve","species","new.metric","pred","SE","lwr","upr")]
fits <- unique(fits)

# produce mean population level prediction for interventions across species
sum.fits <- aggregate(fits[,c("pred","SE","lwr","upr")], by=list(AE.reserve=fits$AE.reserve, new.metric=fits$new.metric), mean)

plotdat <- sum.fits

plotdat <- merge(plotdat, unique(origdat[,c("AE.level","reserve.desig","AE.reserve")]), by="AE.reserve")
plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE.level, plotdat$new.metric),]


# set up colours and plotting points for plot
cols <- rev(grey(seq(from=0,to=1,length.out = 5)))
AEreserve.combo <- unique(dat[,c("AE.level","reserve.desig")])
AEreserve.combo <- AEreserve.combo[order(AEreserve.combo$reserve.desig,AEreserve.combo$AE.level),][-1,]
colpch <- data.frame(cols, AEreserve.combo)
# colpch <- data.frame(cols, AE.level=rep(levels(plotdat$AE.level), times=2), reserve.desig=rep(levels(plotdat$reserve.desig), each=3))
pchs <- data.frame(pch=c(21,22,23), new.metric=levels(plotdat$new.metric))

plotdat <- merge(plotdat, colpch, by=c("AE.level","reserve.desig"))
plotdat <- merge(plotdat, pchs, by=c("new.metric"))

plotdat <- plotdat[order(plotdat$reserve.desig, plotdat$AE.level, plotdat$new.metric),]


###---- Output plot ----###

setwd(outputwd)
png("2c_high level combination intervention success_metric.png", res=300, height=10, width=22, units="in", pointsize=20)

par(mar=c(7,6,2,2))

x <- c(1:nrow(plotdat))

plot(plotdat$pred~x, pch=plotdat$pch, bg=as.character(plotdat$cols), cex=1.5, ylim=c(-0.5,1.5), xaxt="n", xlab="", ylab="", las=1, bty="n")
arrows(x, plotdat$pred, x, plotdat$lwr, angle=90, length=0.05)
arrows(x, plotdat$pred, x, plotdat$upr, angle=90, length=0.05)
abline(h=0.05, lty=3, lwd=2)
# abline(v=12.5, lty=3, lwd=2)
# axis(1, x, labels=rep(c("no AES","basic-level \n AES","higher-level \n AES"), times=2), tick=TRUE, cex.axis=0.8)
axis(1, seq(from=2, by=length(levels(plotdat$new.metric)), length.out=5), labels=rep("", 5), tick=TRUE)
text(seq(from=2, by=length(levels(plotdat$new.metric)), length.out=5), par("usr")[3]*1.15, srt = 0, pos=1, xpd = TRUE, labels=c("basic-level AES\n no nature reserve","higher-level AES\n no nature reserve", "no AES \n nature reserve", "basic-level AES\n nature reserve", "higher-level AES\n nature reserve"), cex=1)
# text(c(4,16), par("usr")[3]*2, srt = 0, pos=4, xpd = TRUE, labels=c("no nature reserve/designation", "nature reserve/designation"), font=2, cex=1.2)
title(ylab="Predicted probability of success \n (significant positive impact)", cex.lab=1.5, font=2, line=3)
title(xlab="Intervention combination", cex.lab=1.5, font=2, line=5.5)

legend(min(x)-0.5,1.2, legend=unique(plotdat$new.metric), pch=plotdat$pch, pt.cex=1.2, bty="n", xpd=TRUE)

dev.off()





######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################



# #-----------------   HIGH-LEVEL INTERVENTION SUCCESS by metric  ------------------
# 
# mdat <- subset(dat, high.int.used==1)
# 
# mdat <- subset(mdat, new.metric!="survival" & new.metric!="recruitment" & new.metric!="productivity")
# mdat <- droplevels(mdat)
# 
# newlevels <- data.frame(unique(mdat[,c("AE.level","reserve.desig")]), AE.reserve=c("AE.basic-no reserve","AE.higher-no reserve","no AE-reserve","AE.basic-reserve","AE.higher-reserve"))
# mdat <- merge(mdat, newlevels, by=c("AE.level","reserve.desig"))
# mdat$AE.reserve <- relevel(mdat$AE.reserve, ref="no AE-reserve")
# 
# m.high.metric <- lme(success ~ AE.reserve*new.metric + species, random = ~1|reference, data=mdat)
# 
# summary(m.high.metric)
# 
# setwd(outputwd)
# sink(paste("model output_2c.txt", sep=" "))
# cat("\n########==========  Success of higher-level interventions combined, by metric ==========########\n", sep="\n")
# print(summary(m.high.metric))
# 
# sink()
# 
# ### Save individual interventions models
# saveRDS(m.high.metric, file=paste(workspacewd, "models_2c.rds", sep="/"))
# 
# ### Save dataset for 0b models
# saveRDS(mdat, file=paste(workspacewd, "model dataset_2c.rds", sep="/"))
# 
# # had to remove productivity from the metrics evaluated because of singularities in model

# #-----------------   HIGH-LEVEL INTERVENTION SUCCESS by habitat  ------------------
# 
# out <- list()
# for(i in 1:length(mgmtvars)) {
#   out[[i]] <- table(dat$newhabitat, dat[,mgmtvars[i]])
# }
# names(out) <- mgmtvars
# out
#
# $AE
# 
# none applied
# arable       55      46
# pastoral    205     208
# unenclosed   71       8
# 
# $AE.level
# 
# none basic higher
# arable       55    11     35
# pastoral    205   132     76
# unenclosed   71     8      0
# 
# $reserve.desig
# 
# none applied
# arable      101       0
# pastoral    274     139
# unenclosed   59      20

# m.high.hab <- lme(success ~ AE*newhabitat + reserve.desig*newhabitat + species, random = ~1|reference, data=dat)

##### model doesn't run because of singularities (probably due to low-no unenclosed habitats using AES and no arable habitats in reserves) ####

#-----------------   SPECIFIC INTERVENTION SUCCESS by metric  ------------------

##### model doesn't run because of singularities  ####

# # remove categories with too few records, which make the model collapse
# mdat <- subset(dat, new.metric!="survival" & new.metric!="recruitment" & new.metric!="occupancy" & new.metric!="abundance change")
# mdat <- subset(mdat, grazing!="reduced" & mowing!="reduced" & fertpest!="applied" & water!="reduced")
# 
# # subset data to only use records where specific interventions were tested
# mdat <- subset(mdat, spec.int.used==1)
# mdat <- droplevels(mdat)
# 
# 
# # identify which categories have low numbers
# out <- list()
# for(i in 1:length(mgmtvars)) {
#   out[[i]] <- table(mdat$new.metric, mdat[,mgmtvars[i]])
# }
# names(out) <- mgmtvars
# out
# 
# # m.spec.metric <- lme(success ~ new.metric*mowing + new.metric*grazing + new.metric*fertpest + new.metric*nest.protect + new.metric*predator.control + new.metric*water + species, random = ~1|reference, data=mdat)




#-----------------   SPECIFIC INTERVENTION SUCCESS by habitat  ------------------

##### model doesn't run because of singularities ####