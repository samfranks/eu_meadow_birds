##################################################################################
#
#     Step 2: EU meadow birds meta-analysis - SUMMARY OF DATA IN REVIEW DATABASE
#
##################################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016

# =================================  SET LOGIC STATEMENTS  ====================



# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

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
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


# =================================  LOAD DATA  ===============================

source(paste(scriptswd, "source_model data preparation.R", sep="/"))



# =================== Study frequency of unique interventions/intervention combinations ====================

library(plyr)

# select only the unique combinations of interventions tested by a study (i.e. some studies may assess interventions across multiple species/metrics)
# use count() in plyr to count the number of occurrences of unique combinations of variable values in a dataset
subdat <- unique(subset(dat, select=c("reference",mgmtvars)))
temp <- count(subdat, vars = c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water"))

write.csv(temp, paste(outputwd, "intervention combination frequency by study.csv", sep="/"), row.names=FALSE)

# =================================  SUMMARY STATISTICS  ===============================

# total number of studies in dataset
num.studies <- length(unique(dat$reference))

# ---------- Proportion of studies by literature type (grey or primary) ----------

# summarize proportions of studies of different lit types
# create based on a unique dataset of reference id and lit.type
litsum <- table(unique(dat[,c("reference","lit.type")])$lit.type)
litsum.prop <- litsum/num.studies


png(paste(outputwd, "summary_proportion of studies by lit type.png", sep="/"), res=300, height=12, width=12, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(litsum.prop, space=0.1, las=1, col="grey90", ylim=c(0,1), xaxt="n")
text(x, par("usr")[3]-0.04, srt = 0, adj=1, xpd = TRUE, labels = c(names(litsum.prop)))
title(xlab="Literature type", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, litsum.prop+0.02, litsum) # sample sizes for each lit type

dev.off()

# ---------- Proportion of studies by country ----------

# summarize proportions of studies from different countries
# create based on a unique dataset of reference id and country
countrysum <- table(unique(dat[,c("reference","country")])$country)
countrysum.prop <- countrysum/num.studies


png(paste(outputwd, "summary_proportion of studies by country + species.png", sep="/"), res=300, height=20, width=16, units="in", pointsize=20)

# png(paste(outputwd, "summary_proportion of studies by country.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20)

par(mfrow=c(2,1))

par(mar=c(6,5,2,1))
x <- barplot(countrysum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.6), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c(names(countrysum.prop)[-which(names(countrysum.prop) %in% "United Kingdom")], "United \n Kingdom"))
title(xlab="Country", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, countrysum.prop+0.02, countrysum) # sample sizes for each country

# dev.off()


# ---------- Proportion of studies by species ----------

speciessum <- table(unique(dat[,c("reference","species")])$species)
speciessum.prop <- speciessum/num.studies

# png(paste(outputwd, "summary_proportion of studies by species.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(speciessum.prop, space=0.1, las=1, col="grey90", ylim=c(0,1), xaxt="n")
text(x, par("usr")[3]-0.03, srt = 0, pos=1, xpd = TRUE, labels = c( "black-tailed \n godwit",names(speciessum.prop)[-which(names(speciessum.prop) %in% "black-tailed godwit")]))
title(xlab="Species", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, speciessum.prop+0.02, speciessum) # sample sizes for each species

dev.off()


# ---------- Proportion of studies by habitat ----------

d1.hab <- dat
d1.hab$newhabitat <- factor(d1.hab$newhabitat, c("arable","pastoral","unenclosed"), ordered=TRUE)
habsum <- table(unique(d1.hab[,c("reference","newhabitat")])$newhabitat)
habsum.prop <- habsum/num.studies

png(paste(outputwd, "summary_proportion of studies by habitat.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(habsum.prop, space=0.1, las=1, col="grey90", ylim=c(0,1), xaxt="n")
text(x, par("usr")[3]-0.04, srt = 0, xpd = TRUE, labels = c("arable","pastoral","unenclosed"))
title(xlab="Habitat", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, habsum.prop+0.02, habsum) # sample sizes for each habitat

dev.off()


# ---------- Proportion of studies by overall metric ----------

# create an ordered factor of metrics for this summary (particularly for productivity so it's nest, chick, nest+chick)
# don't want to change this factor to ordered for the analysis though, because it will be treated differently than a regular factor in the model

d1.metric <- dat
d1.metric$new.metric <- factor(d1.metric$new.metric, c("abundance/occupancy","abundance/occupancy change","productivity","recruitment","survival"), ordered=TRUE)

metricsum <- table(unique(d1.metric[,c("reference","new.metric")])$new.metric)
metricsum.prop <- metricsum/num.studies

png(paste(outputwd, "summary_proportion of studies by metric.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(metricsum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.8), xaxt="n")
text(x, par("usr")[3]-0.04, srt = 30, adj=1, xpd = TRUE, labels = c("abundance/\noccupancy","abundance/\noccupancy change","productivity","recruitment","survival"))
title(xlab="Study metric", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, metricsum.prop+0.02, metricsum) # sample sizes for each metric

dev.off()

# ---------- Proportion of studies by specific metric (for productivity)----------

# create an ordered factor of metrics for this summary (particularly for productivity so it's nest, chick, nest+chick)

d1.metric <- dat

# pool abundance/occupancy and change metrics
# show different productivity metrics separately
d1.metric$prod.metric <- with(d1.metric, ifelse(d1.metric$overall.metric=="abundance" | d1.metric$overall.metric=="occupancy", "abundance/occupancy", ifelse(d1.metric$overall.metric=="abundance change" | d1.metric$overall.metric=="occupancy change", "abundance/occupancy change", as.character(overall.metric))))

# create ordered factor
d1.metric$prod.metric <- factor(d1.metric$prod.metric, c("abundance/occupancy","abundance/occupancy change","productivity (nest level)", "productivity (chick level)", "productivity (nest + chick)","recruitment","survival"), ordered=TRUE)

metricsumprod <- table(unique(d1.metric[,c("reference","prod.metric")])$prod.metric)
metricsumprod.prop <- metricsumprod/num.studies

png(paste(outputwd, "summary_proportion of studies by productivity metric.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(c(metricsumprod.prop[1:2], metricsum.prop["productivity"], metricsumprod.prop[3:7]), space=0.1, las=1, col="grey90", ylim=c(0,0.8), xaxt="n")
text(x, par("usr")[3]-0.01, srt = 0, pos=1, xpd = TRUE, labels = c("abundance/\noccupancy","abundance/\noccupancy \nchange", "productivity \n(all)", "productivity \n(nest)", "productivity \n(chick)", "productivity \n(nest+chick)","recruitment","survival"))
title(xlab="Study metric", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, c(metricsumprod.prop[1:2]+0.02, metricsum.prop["productivity"]+0.02, metricsumprod.prop[3:7]+0.02), c(metricsumprod[1:2], metricsum["productivity"], metricsumprod[3:7])) # sample sizes for each metric

dev.off()

# ---------- Proportion of studies examining the effect of an intervention type ----------

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

png(paste(outputwd, "summary_proportion of studies by intervention.png", sep="/"), res=300, height=12, width=16, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(intervensum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.6), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c("AES","site\nprotection", "mowing","grazing","fertiliser/\npesticides","nest\nprotection","predator\ncontrol","water\nmanagement"))
title(xlab="Intervention", font=2, cex.lab=1.2, line=4.5)
title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
text(x, intervensum.prop+0.02, intervensum) # sample sizes for each intervention type

dev.off()


# ---------- Types of management interventions comprised within AES types ----------


AES.level <- c("basic","higher")

for (j in 1:length(AES.level)) {
  
  ### filter dataset to only those studies which test AES effects
  subdat <- subset(dat, AE.level==AES.level[j])
  
  num.studies <- length(unique(subdat$reference))
  
  
  ### determine how many AES = BASIC studies evaluated the effect of the intervention (repeat for each intervention)
  # create blank objects
  intervensum <- numeric()
  intervensum.prop <- numeric()
  intervensum.level <- list()
  intervensum.level.prop <- list()
  
  # mgmtvars to evluate, minus AE.level
  eval.mgmtvars <- mgmtvars[-which(mgmtvars %in% c("AE","AE.level","reserve.desig"))]
  
  for (i in 1:length(eval.mgmtvars)) {
    
    # number of unique cases (i.e. unique studies) where mgmtvar level != 'none'
    x <- unique(subdat[,c("reference",eval.mgmtvars[i])]) # unique references and levels of the intervention
    y <- x[x[eval.mgmtvars[i]] != "none",] # remove the cases where intervention not evaluated
    intervensum[i] <- length(unique(y$reference)) # the number of studies that evaluated the intervention
    intervensum.prop[i] <- intervensum[i]/num.studies
    intervensum.level[[i]] <- table(y[,eval.mgmtvars[i]])
    intervensum.level.prop[[i]] <- intervensum.level[[i]]/num.studies
  }
  
  names(intervensum) <- eval.mgmtvars
  names(intervensum.prop) <- eval.mgmtvars
  names(intervensum.level) <- eval.mgmtvars
  names(intervensum.level.prop) <- eval.mgmtvars
  
  # count number of unique studies that evaluate no interventions as part of AES, add to plotting dataset
  no.interventions.tested <- subset(subdat, mowing=="none" & grazing=="none" & fertpest=="none" & water=="none" & nest.protect=="none" & predator.control=="none")
  no.interventions.tested <- unique(no.interventions.tested[,c("reference",eval.mgmtvars)])
  
  intervensum <- c(nrow(no.interventions.tested), intervensum)
  intervensum.prop <- c(nrow(no.interventions.tested)/num.studies, intervensum.prop)
  
  # remove AE level from being plotted
  # intervensum.prop <- intervensum.prop[-which(names(intervensum.prop) %in% "AE.level")]
  # intervensum <- intervensum[-which(names(intervensum) %in% "AE.level")]
  
  png(paste(outputwd, "/summary_proportion of management interventions used by AES level ", AES.level[j], ".png", sep=""), res=300, height=12, width=16, units="in", pointsize=20)
  
  par(mar=c(6,5,2,1))
  x <- barplot(intervensum.prop, space=0.1, las=1, col="grey90", ylim=c(0,1), xaxt="n")
  text(x, par("usr")[3]-0.02, srt = 0, pos=1, xpd = TRUE, labels = c("no\ninterventions\nevaluated","mowing","grazing","agrochemicals","nest\nprotection","predator\ncontrol","water\nmanagement"))
  title(xlab="Intervention", font=2, cex.lab=1.2, line=4.5)
  title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
  text(x, intervensum.prop+0.02, intervensum) # sample sizes for each intervention type
  
  if (AES.level[j]=="basic") mtext("a)", side=3, adj=0, line=0.5)
  if (AES.level[j]=="higher") mtext("b)", side=3, adj=0, line=0.5)
  
  
  dev.off()
  
}


# ---------- Summary statistics text output  ----------

setwd(outputwd)
sink(paste("summary statistics output.txt", sep=" "))

cat("\n\n####### Literature summary #######\n")
print(litsum.prop)
cat("\n\n####### Country summary #######\n")
print(countrysum.prop)
cat("\n\n####### Intervention summary #######\n")
print(intervensum.prop)
cat("\n\n####### Species summary #######\n")
print(speciessum.prop)
cat("\n\n####### Metric summary #######\n")
print(metricsum.prop)
cat("\n\n####### Metric summary (specific productivity metrics) #######\n")
print(metricsumprod.prop)
cat("\n\n####### Habitat summary #######\n")
print(habsum.prop)

sink()



# #---------- Proportion of studies by intervention level  ----------
# 
# ### plot levels of different management types attempted (e.g. applied/reduced)
# level.prop.all <- as.data.frame(do.call(rbind, intervensum.level.prop))
# 
# level.prop.all$level1 <- level.prop.all$basic
# level.prop.all$level2 <- level.prop.all$higher
# # level.prop.all$level2 <- ifelse(level.prop.all$basic==level.prop.all$higher, 0, level.prop.all$higher)
# 
# level.prop.all <- level.prop.all[-which(rownames(level.prop.all) %in% c("AE","reserve.desig","nest.protect")),]
# 
# level.all <- as.data.frame(do.call(rbind, intervensum.level))
# level.all <- level.all[-which(rownames(level.all) %in% c("AE","reserve.desig","nest.protect")),]
# 
# level.sum <- apply(level.all,1,sum) # sum the number of studies which examine impact of different levels of each intervention (the same study could look at 2 levels of the same intervention, so totals will not equal above, which is the proportion of studies examining each intervention type; in fact, this is only the case for one study that looks at both basic and higher level AES effects)
# 
# # level.prop.all <- level.prop.all[-which(rownames(level.prop.all) %in% "AE.level"),]
# 
# (level.prop.all <- t(level.prop.all[,c("level1","level2")]))
# 
# 
# png(paste(outputwd, "summary_proportion of studies by intervention_level.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)
# 
# par(mar=c(6,5,2,1))
# x <- barplot(level.prop.all, las=1, col=c("grey90","grey30"), beside=FALSE, ylim=c(0,0.6), xaxt="n")
# text(x, par("usr")[3]-0.02, srt = 30, adj=1, xpd = TRUE, labels = c("AES level","mowing","grazing","fertiliser/ \n pesticides","predator \n control","water \n management"))
# title(xlab="Management intervention", font=2, cex.lab=1.2, line=4.5)
# title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
# text(x, apply(level.prop.all,2,sum)+0.02, level.sum) # sample sizes for each intervention type
# 
# dev.off()
# 
# #---------- Proportion of studies examining the effect of an intervention type, by species ----------
# 
# ### determine how many studies evaluated the effect of the intervention (repeat for each intervention)
# # create blank objects
# intervensum <- list()
# intervensum.prop <- list()
# 
# # mgmtvars to evluate, minus AE.level
# eval.mgmtvars <- mgmtvars[-which(mgmtvars %in% "AE.level")]
# 
# for (i in 1:length(eval.mgmtvars)) {
#   
#   # number of unique cases (i.e. unique studies) where mgmtvar level != 'none'
#   x <- unique(dat[,c("reference","species",eval.mgmtvars[i])]) # unique references and levels of the intervention
#   y <- x[x[eval.mgmtvars[i]] != "none",] # remove the cases where intervention not evaluated
#   intervensum[[i]] <- table(y$species)
# }
# 
# intervensum.all <- do.call(cbind, intervensum)
# colnames(intervensum.all) <- eval.mgmtvars
# 
# intervensum.prop.all <- intervensum.all/num.studies
# 
# # creates grey-scale colour vector for plotting, but randomly shuffled so darks don't end up next to each other
# set.seed(3)
# # colourvec <- sample(grey(seq(from=0,to=1,length.out = 8)), 8)
# colourvec <- grey(seq(from=1,to=0,length.out = 8))
# 
# png(paste(outputwd, "summary_proportion of studies by intervention by species.png", sep="/"), res=300, height=12, width=20, units="in", pointsize=20)
# 
# par(mar=c(6,5,2,1))
# x <- barplot(intervensum.prop.all, beside=TRUE, las=1, col=colourvec, ylim=c(0,0.4), xaxt="n")
# text(apply(x,2,mean), par("usr")[3]-0.02, srt = 0, xpd = TRUE, labels = c("AES","nature reserve/ \n designation", "mowing","grazing","fertiliser/ \n pesticides","nest \n protection","predator \n control","water \n management"))
# title(xlab="Management intervention", font=2, cex.lab=1.2, line=3)
# title(ylab=paste("Proportion of total studies (n=", num.studies, ")", sep=""), font=2, cex.lab=1.2, line=3)
# legend("topright",legend=rownames(intervensum.prop.all), pch=rep(22,8), col="black",pt.bg=colourvec, cex=1, pt.cex=1.5, bty="n")
# # text(x, intervensum.prop.all+0.02, intervensum.all) # sample sizes for each intervention type
# 
# dev.off()

