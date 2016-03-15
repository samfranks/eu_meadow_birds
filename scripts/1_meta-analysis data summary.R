#######################################################
#
#     EU meadow birds meta-analysis
#
#######################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal")

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


#=================================  LOAD & CLEAN DATA  ===============================

d0 <- read.csv(paste(datawd, "meadow birds data extraction template_final_primary.csv", sep="/"), header=TRUE, skip=1)

#------- Meta-data reference for studies -------------
# create a meta-data reference file for studies with reference numbers, reference name, summary, country, region
metadat0 <- unique(d0[,c("reference.number","reference","one.sentence.summary","score","country","region1","region2")])

#------- Clean dataset -----------

# columns required
d0.1 <- d0[,c("reference.number","record.number","country","region1","habitat","habitat1","habitat2","start.year","end.year","type.of.study","species","assemblage","agri.environment","basic.agri.environment", "targeted.agri.environment..wader.specific.or.higher.level.", "site.protection...nature.reserve","site.protection...designation", "mowing","grazing","fertilizer","herbicides...pesticides","nest.protection...agricultural.activities","nest.protection...predation..enclosures.or.exclosures.",     "ground.water.management..drainage.inhibited.","wet.features...surface.water.management","predator.control","other.mgmt", "management.notes","overall.metric","specific.metric","unit.of.analysis","sample.size","analysis.type.1","analysis.type.2","analysis.type.details","significant.effect..Y.N..U.","direction.of.effect..positive...negative...none.","reference.metric.before.management","metric.after.management","log.effect.size.", "effect.size..where.reported.","values.obtained.from.plot.")]

# rename to easier variables 
d0.2 <- d0.1
names(d0.2) <- c("reference","record","country","region1","habitat","habitat1","habitat2","start.year","end.year","study.type","species","assemblage","AE","basic.AE","higher.AE","reserve","designation","mowing","grazing","fertilizer","pesticide","nest.protect.ag","nest.protect.predation","groundwater.drainage","surface.water","predator.control","other.mgmt","mgmt.notes","overall.metric","specific.metric","analysis.unit","sample.size","analysis1","analysis2","analysis3","sig","effect.dir","metric.before","metric.after","log.effect.size","percent.change.effect.size","values.from.plot")

# management intervention variables
mgmtvars <- c("AE","basic.AE","higher.AE","reserve","designation","mowing","grazing","fertilizer","pesticide","nest.protect.ag","nest.protect.predation","groundwater.drainage","surface.water","predator.control","other.mgmt")

d0.3 <- d0.2

# recode manamgement vars as characters to be able to use string substitution find and replace to create generic applied, restricted, removed levels for all management types
d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, as.character)

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
  gsub("applied site scale", "applied", x, fixed=FALSE)
})

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("applied landscape scale", "applied", x, fixed=FALSE)
  })

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("restricted site scale", "restricted", x)
  })

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("restricted landscape scale", "restricted", x)
  })

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("removed site scale", "removed", x)
  })

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("removed landscape scale", "removed", x)
  })

# plug 'none' into all the blanks where management intervention not used
for (i in 1:length(mgmtvars)) {
  d0.3[d0.3[,mgmtvars[i]]=="",mgmtvars[i]] <- "none"
}

# change management vars back to factors for analysis
d0.4 <- d0.3
# d0.4[,mgmtvars] <- apply(d0.4[,mgmtvars], 2, function(x) as.factor(x)) # this line won't convert back to factors for some reason!
for (i in 1:length(mgmtvars)) {
  d0.4[,mgmtvars[i]] <- as.factor(d0.4[,mgmtvars[i]])
}
summary(d0.4)

### Add some additional grouping variables and success variable

# success variable defined as 1 = significant positive effect, 0 = neutral or negative effect
d0.4$success <- ifelse(d0.4$sig=="Y" & d0.4$effect.dir=="positive", 1, 0) # success variable (defin)

# group fertilizer and pesticides into single variable
d0.4$fertpest <- ifelse(d0.4$fertilizer=="applied" | d0.4$pesticide=="applied", 
                        "applied",
                        ifelse(d0.4$fertilizer=="restricted" | d0.4$pesticide=="restricted", 
                               "restricted", 
                               ifelse(d0.4$fertilizer=="removed" | d0.4$pesticide=="removed", 
                                      "removed", "none"))
                        )
# group groundwater.drainage and surface.water into single variable meaning 'more water'
# restricted/removed groundwater drainage equates to more water (same as applying surface water)
# combinations of drainage/surface water in dataset
unique(d0.4[,c("groundwater.drainage","surface.water")])
d0.4$water <- ifelse(d0.4$groundwater.drainage=="restricted" | d0.4$groundwater.drainage=="removed" & d0.4$surface.water=="applied", "applied", ifelse(d0.4$groundwater.drainage=="restricted" | d0.4$groundwater.drainage=="removed", "applied", ifelse(d0.4$surface.water=="applied", "applied", ifelse(d0.4$groundwater.drainage=="applied","restricted","none"))))

# group nest protection (predation) with predator control (more sensible than grouping it with nest protection for agriculture given predation measures are more likely to go together)
unique(d0.4[,c("nest.protect.ag","nest.protect.predation","predator.control")])
d0.4$predation.reduction <- ifelse(d0.4$nest.protect.predation=="applied" | d0.4$predator.control=="applied", "applied", ifelse(d0.4$predator.control=="restricted", "restricted", ifelse(d0.4$predator.control=="removed", "removed","none")))

# group reserves and site designations
d0.4$reserve.desig <- ifelse(d0.4$reserve=="applied" | d0.4$designation=="applied", "applied", "none")

# create a AE-level variable (with basic and higher as levels) for analysis 1a
# if no info was provided on type of AES, then assume it was basic rather than higher-level or targetted
d0.4$AE.level <- ifelse(d0.4$higher.AE=="applied", "higher", ifelse(d0.4$AE=="none", "none", "basic"))

# calculate study duration variable
d0.4$study.length <- d0.4$end.year - d0.4$start.year + 1

#------------- Definitive dataset --------------

# final dataset for analysis 
d1 <- d0.4

# new set of management variables
mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect.ag","predation.reduction","water")



#=================================  SUMMARY STATISTICS  ===============================

#---------- Proportion of studies by country ----------

# summarize proportions of studies from different countries
# create based on a unique dataset of reference id and country
countrysum <- table(unique(d1[,c("reference","country")])$country)
countrysum.prop <- countrysum/62


png(paste(outputwd, "summary_proportion of studies by country.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(countrysum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.6), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 45, adj=1, xpd = TRUE, labels = c(names(countrysum.prop)[-which(names(countrysum.prop) %in% "United Kingdom")], "United \n Kingdom"))
title(xlab="Country", font=2, cex.lab=1.2, line=4.5)
title(ylab="Proportion of total studies (n=62)", font=2, cex.lab=1.2, line=3)
text(x, countrysum.prop+0.02, countrysum) # sample sizes for each country

dev.off()


#---------- Proportion of studies by species ----------

speciessum <- table(unique(d1[,c("reference","species")])$species)
speciessum.prop <- speciessum/62

png(paste(outputwd, "summary_proportion of studies by species.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(speciessum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.8), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 45, adj=1, xpd = TRUE, labels = c( "black-tailed \n godwit",names(speciessum.prop)[-which(names(speciessum.prop) %in% "black-tailed godwit")]))
title(xlab="Species", font=2, cex.lab=1.2, line=4.5)
title(ylab="Proportion of total studies (n=62)", font=2, cex.lab=1.2, line=3)
text(x, speciessum.prop+0.02, speciessum) # sample sizes for each species

dev.off()

#---------- Proportion of studies by overall metric ----------

# create an ordered factor of metrics for this summary (particularly for productivity so it's nest, chick, nest+chick)
# don't want to change this factor to ordered for the analysis though, because it will be treated differently than a regular factor in the model

d1.metric <- d1
d1.metric$overall.metric <- factor(d1.metric$overall.metric, c("abundance","abundance change","occupancy","occupancy change","productivity (nest level)","productivity (chick level)","productivity (nest + chick)","recruitment","survival"), ordered=TRUE)

metricsum <- table(unique(d1.metric[,c("reference","overall.metric")])$overall.metric)
metricsum.prop <- metricsum/62

png(paste(outputwd, "summary_proportion of studies by metric.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(metricsum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.8), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 45, adj=1, xpd = TRUE, labels = c("abundance","abundance \n change","occupancy","occupancy \n change","productivity \n (nest-level)","productivity \n (chick-level)","productivity \n (nest + chick-level)", "recruitment","adult survival"), cex=0.8)
title(xlab="Study metric", font=2, cex.lab=1.2, line=4.5)
title(ylab="Proportion of total studies (n=62)", font=2, cex.lab=1.2, line=3)
text(x, metricsum.prop+0.02, metricsum) # sample sizes for each metric

dev.off()

#---------- Proportion of studies examining the effect of an intervention type ----------

### determine how many studies evaluated the effect of the intervention (repeat for each intervention)
# create blank objects
intervensum <- numeric()
intervensum.prop <- numeric()

# mgmtvars to evluate, minus AE.level
eval.mgmtvars <- mgmtvars[-which(mgmtvars %in% "AE.level")]

for (i in 1:length(eval.mgmtvars)) {
  
  # number of unique cases (i.e. unique studies) where mgmtvar level != 'none'
  x <- unique(d1[,c("reference",eval.mgmtvars[i])]) # unique references and levels of the intervention
  y <- x[x[eval.mgmtvars[i]] != "none",] # remove the cases where intervention not evaluated
  intervensum[i] <- nrow(y) # the number of studies that evaluated the intervention
  intervensum.prop[i] <- intervensum[i]/62
  
}

names(intervensum) <- eval.mgmtvars
names(intervensum.prop) <- eval.mgmtvars

png(paste(outputwd, "summary_proportion of studies by intervention.png", sep="/"), res=300, height=12, width=15, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(intervensum.prop, space=0.1, las=1, col="grey90", ylim=c(0,0.6), xaxt="n")
text(x, par("usr")[3]-0.02, srt = 45, adj=1, xpd = TRUE, labels = c("AES","nature reserve/ \n designation", "mowing","grazing","fertiliser/ \n pesticides","nest \n protection","predation \n reduction","water \n management"))
title(xlab="Management intervention", font=2, cex.lab=1.2, line=4.5)
title(ylab="Proportion of total studies (n=62)", font=2, cex.lab=1.2, line=3)
text(x, intervensum.prop+0.02, intervensum) # sample sizes for each intervention type

dev.off()

#---------- Proportion of studies examining the effect of an intervention type, by species ----------

### determine how many studies evaluated the effect of the intervention (repeat for each intervention)
# create blank objects
intervensum <- list()
intervensum.prop <- list()

# mgmtvars to evluate, minus AE.level
eval.mgmtvars <- mgmtvars[-which(mgmtvars %in% "AE.level")]

for (i in 1:length(eval.mgmtvars)) {
  
  # number of unique cases (i.e. unique studies) where mgmtvar level != 'none'
  x <- unique(d1[,c("reference","species",eval.mgmtvars[i])]) # unique references and levels of the intervention
  y <- x[x[eval.mgmtvars[i]] != "none",] # remove the cases where intervention not evaluated
  intervensum[[i]] <- table(y$species)
}

intervensum.all <- do.call(cbind, intervensum)
colnames(intervensum.all) <- eval.mgmtvars

intervensum.prop.all <- intervensum.all/62

# creates grey-scale colour vector for plotting, but randomly shuffled so darks don't end up next to each other
set.seed(3)
# colourvec <- sample(grey(seq(from=0,to=1,length.out = 8)), 8)
colourvec <- grey(seq(from=1,to=0,length.out = 8))

png(paste(outputwd, "summary_proportion of studies by intervention by species.png", sep="/"), res=300, height=12, width=20, units="in", pointsize=20)

par(mar=c(6,5,2,1))
x <- barplot(intervensum.prop.all, beside=TRUE, las=1, col=colourvec, ylim=c(0,0.4), xaxt="n")
text(apply(x,2,mean), par("usr")[3]-0.02, srt = 0, xpd = TRUE, labels = c("AES","nature reserve/ \n designation", "mowing","grazing","fertiliser/ \n pesticides","nest \n protection","predation \n reduction","water \n management"))
title(xlab="Management intervention", font=2, cex.lab=1.2, line=3)
title(ylab="Proportion of total studies (n=62)", font=2, cex.lab=1.2, line=3)
legend("topright",legend=rownames(intervensum.prop.all), pch=rep(22,8), col="black",pt.bg=colourvec, cex=1, bty="n")
# text(x, intervensum.prop.all+0.02, intervensum.all) # sample sizes for each intervention type

dev.off()


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################




