############################################################################################
#
#     Step 1: EU meadow birds meta-analysis - DATA PREPARATION FROM EXTRACTED DATABASE
#
############################################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016

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
outputwd <- paste(parentwd, "output/revision Dec 2016", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


#=================================  LOAD & CLEAN DATA  ===============================

# d0 <- read.csv(paste(datawd, "meadow birds data extraction template_final_primary.csv", sep="/"), header=TRUE, skip=1)

d0 <- read.csv(paste(datawd, "Meadow birds data extraction template_primary and grey_standardized_FINAL.csv", sep="/"), header=TRUE)

#------- Meta-data reference for studies -------------
# create a meta-data reference file for studies with reference numbers, reference name, summary, country, region
metadat0 <- unique(d0[,c("reference.number","reference","literature.type","one.sentence.summary","score","country","region1","region2")])

#------- Clean dataset -----------

# columns required
d0.1 <- d0[,c("reference.number","record.number","literature.type","score","country","region1","habitat","habitat1","habitat2","start.year","end.year","type.of.study","species","assemblage","agri.environment","basic.agri.environment", "targeted.agri.environment..wader.specific.or.higher.level.", "site.protection...nature.reserve","site.protection...designation", "mowing","grazing","fertilizer","herbicides...pesticides","nest.protection...agricultural.activities","nest.protection...predation..enclosures.or.exclosures.",     "ground.water.management..drainage.inhibited.","wet.features...surface.water.management","predator.control","other.mgmt", "management.notes","overall.metric","specific.metric","reference.metric.before.management","metric.after.management","standardized.metric","standardisation.calculation","stand..reference.metric.before.management","stand..metric.after.management", "stand..effect.size","significant.effect..Y.N..U.","direction.of.effect..positive...negative...none.","unit.of.analysis","sample.size","analysis.type.1","analysis.type.2","analysis.type.details","log.effect.size.", "effect.size..where.reported.","values.obtained.from.plot.")]

# rename to easier variables 
d0.2 <- d0.1
names(d0.2) <- c("reference","record","lit.type","score","country","region1","habitat","habitat1","habitat2","start.year","end.year","study.type","species","assemblage","AE","basic.AE","higher.AE","reserve","designation","mowing","grazing","fertilizer","pesticide","nest.protect.ag","nest.protect.predation","groundwater.drainage","surface.water","predator.control","other.mgmt","mgmt.notes","overall.metric","specific.metric","metric.before","metric.after","stan.metric","stan.calc","stan.metric.before","stan.metric.after","stan.effect.size","sig","effect.dir","analysis.unit","sample.size","analysis1","analysis2","analysis3","log.effect.size","percent.change.effect.size","values.from.plot")

# management intervention variables
mgmtvars <- c("AE","basic.AE","higher.AE","reserve","designation","mowing","grazing","fertilizer","pesticide","nest.protect.ag","nest.protect.predation","groundwater.drainage","surface.water","predator.control","other.mgmt")

### exlude studies 2 and 36
# 2: remove this reference (Kruk et al. 1997) as it doesn't really measure a population or demographic metric
# 36: remove this reference (Kleijn et al. 2004) as it pools an assessment of conservation across multiple species

d0.2 <- subset(d0.2, reference!=36) # remove this reference (Kruk et al. 1997) as it doesn't really measure a population or demographic metric
d0.2 <- subset(d0.2, reference!=2) # remove this reference (Kleijn et al. 2004) as it pools an assessment of conservation across multiple species
d0.2 <- droplevels(d0.2)

d0.3 <- d0.2

# recode manamgement vars as characters to be able to use string substitution find and replace to create generic applied, restricted, removed levels for all management types
d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, as.character)

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
  gsub("applied site scale", "applied", x)
})

d0.3[,mgmtvars] <- apply(d0.3[,mgmtvars], 2, function(x) {
    gsub("applied landscape scale", "applied", x)
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

# recode sample size as small, medium, large
d0.3$sample.size <- ifelse(d0.3$sample.size=="small (< 30)", "small", ifelse(d0.3$sample.size=="medium (30-100)", "medium", "large"))



# redefine dataset
d0.4 <- d0.3

# # change management vars back to factors for analysis
# # d0.4[,mgmtvars] <- apply(d0.4[,mgmtvars], 2, function(x) as.factor(x)) # this line won't convert back to factors for some reason!
# for (i in 1:length(mgmtvars)) {
#   d0.4[,mgmtvars[i]] <- as.factor(d0.4[,mgmtvars[i]])
# }
# summary(d0.4)

#---------- Add some additional grouping variables -----------

# group fertilizer and pesticides into single variable
d0.4$fertpest <- ifelse(d0.4$fertilizer=="applied" | d0.4$pesticide=="applied", "applied", ifelse(d0.4$fertilizer=="restricted" | d0.4$pesticide=="restricted", "restricted", ifelse(d0.4$fertilizer=="removed" | d0.4$pesticide=="removed", "removed", "none")))

# group groundwater.drainage and surface.water into single variable meaning 'more water'
# restricted/removed groundwater drainage equates to more water (same as applying surface water)
# combinations of drainage/surface water in dataset
unique(d0.4[,c("groundwater.drainage","surface.water")])
d0.4$water <- ifelse(d0.4$groundwater.drainage=="restricted" | d0.4$groundwater.drainage=="removed" & d0.4$surface.water=="applied", "applied", ifelse(d0.4$groundwater.drainage=="restricted" | d0.4$groundwater.drainage=="removed", "applied", ifelse(d0.4$surface.water=="applied", "applied", ifelse(d0.4$groundwater.drainage=="applied","restricted","none"))))

# group nest protection (predation and agricultural) variables together
unique(d0.4[,c("nest.protect.ag","nest.protect.predation")])
d0.4$nest.protect <- ifelse(d0.4$nest.protect.predation=="applied" | d0.4$nest.protect.ag=="applied", "applied","none")

# # group nest protection (predation) with predator control (more sensible than grouping it with nest protection for agriculture given predation measures are more likely to go together)
# unique(d0.4[,c("nest.protect.ag","nest.protect.predation","predator.control")])
# d0.4$predation.reduction <- ifelse(d0.4$nest.protect.predation=="applied" | d0.4$predator.control=="applied", "applied", ifelse(d0.4$predator.control=="restricted", "restricted", ifelse(d0.4$predator.control=="removed", "removed","none")))

# group reserves and site designations
d0.4$reserve.desig <- ifelse(d0.4$reserve=="applied" | d0.4$designation=="applied", "applied", "none")

# create a AE-level variable (with basic and higher as levels) for analysis 1a
# if no info was provided on type of AES, then assume it was basic rather than higher-level or targetted
d0.4$AE.level <- ifelse(d0.4$higher.AE=="applied", "higher", ifelse(d0.4$AE=="none", "none", "basic"))

# calculate study duration variable
d0.4$study.length <- d0.4$end.year - d0.4$start.year + 1

# add some overall metrics which lump all productivity metrics, all abundance metrics, all occupancy metrics
d0.4$metric <- ifelse(grepl("productivity", d0.4$overall.metric), "productivity", ifelse(grepl("abundance", d0.4$overall.metric), "abundance", ifelse(grepl("recruitment", d0.4$overall.metric), "recruitment", ifelse(grepl("survival", d0.4$overall.metric), "survival", "occupancy"))))


#------------- Change the predator.control level for studies 5 & 10 ---------------

# these 2 studies both deal with the effects of a halt in predator control/game-keepering on grouse moors and the impacts on wader populations
# kind of a reverse of what the conservation measure would normally be (control applied), so reverse the level of predator control to 'applied' and change the direction of the effect (but obviously leave the significance)
# create 5 new records for these studies (2 and 3 each), then add them to the dataset WITH THEIR EFFECT SIZES REMOVED so there is no confusion

temp <- d0.4[d0.4$reference=="5" | d0.4$reference=="10",]

newtemp <- temp

# change predator control to applied
newtemp$predator.control <- "applied"

# change positives to negatives and vice versa
newtemp$effect.dir <- ifelse(newtemp$effect.dir=="positive","negative","positive")

newtemp$metric.before <- temp$metric.after
newtemp$metric.after <- temp$metric.before
newtemp$stan.metric.before <- temp$stan.metric.after
newtemp$stan.metric.after <- temp$stan.metric.before
newtemp$stan.effect.size <- (newtemp$stan.metric.after - newtemp$stan.metric.before)/abs(newtemp$stan.metric.before)

newtemp[,c("log.effect.size","percent.change.effect.size")] <- ""


# remove the original records from the dataset and add these new ones in
d0.4 <- d0.4[-which(d0.4$reference %in% c("5","10")),]
d0.4 <- rbind(d0.4, newtemp)

#------------ Add the success/failure/outcome variables --------------

# success variable defined as 1 = significant positive effect, 0 = neutral or negative effect
d0.4$success <- ifelse(d0.4$sig=="Y" & d0.4$effect.dir=="positive", 1, 0) # success variable

# failure variable defined as 1 = significant negative effect, 0 = neutral or positive effect
d0.4$failure <- ifelse(d0.4$sig=="Y" & d0.4$effect.dir=="negative", 1, 0) # failure variable

# outcome variable: -1 = significant negative, 0 = no effect, 1 = significant positive
d0.4$outcome <- ifelse(d0.4$sig=="Y" & d0.4$effect.dir=="positive", 1, ifelse(d0.4$sig=="Y" & d0.4$effect.dir=="negative", -1, 0)) # success variable) # success variable

#------------- Definitive dataset --------------

# final dataset for analysis 
d1 <- d0.4

# new set of management variables
mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# convert removed or restricted levels of the management vars (all but AE.level) to a single level of removed/restricted
# use find and replace with gsub
d1[,mgmtvars] <- apply(d1[,mgmtvars], 2, function(x) {
  gsub("removed", "reduced", x)
})

d1[,mgmtvars] <- apply(d1[,mgmtvars], 2, function(x) {
  gsub("restricted", "reduced", x)
})

### Save definitive dataset
saveRDS(d1, file=paste(workspacewd, "meadow birds analysis dataset_full.rds", sep="/"))
write.table(d1, file=paste(datawd, "meadow birds analysis dataset_full.txt", sep="/"), row.names=FALSE, quote=FALSE, sep="\t")
write.csv(d1, file=paste(datawd, "meadow birds analysis dataset_full.csv", sep="/"), row.names=FALSE)


