#######################################################
#
#     EU meadow birds meta-analysis
#
#######################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4")

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

dat0 <- readRDS(paste(workspacewd, "meadow birds analysis dataset.rds", sep="/"))

mgmtvars2 <- mgmtvars[-which(mgmtvars %in% "AE")]

# subset dataset for analysis to desired columns only
dat <- subset(dat0, select=c("reference","country","study.length","habitat","species","overall.metric","sample.size","analysis2",mgmtvars2))

# recode management variables as factors instead of characters and specify the reference level


#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars2)) {
  dat[,mgmtvars2[i]] <- as.factor(dat[,mgmtvars2[i]])
  dat[,mgmtvars2[i]] <- relevel(dat[,mgmtvars2[i]], ref="none")
}

#=================================  ANALYSIS  ===============================

# 1a) success of AES and nature reserves
# 1b) success of specific management interventions
# 2a) success of AES and nature reserves by species
# 2b) success of specific management interventions by species
# 3a) success of AES and nature reserves by different metrics (abundance vs reproductive success)
# 3b) success of specific management interventions by different metrics (abundance vs reproductive success)
# 4a) success of AES and nature reserves by reproductive success metrics (nest vs chick)
# 4b) success of specific management interventions by reproductive success metrics (nest vs chick)
# 5a) success of AES and nature reserves by abundance metrics (abundance vs abundance change)
# 5b) success of specific management interventions by abundance metrics (abundance vs abundance change)
# 6a) success of AES and nature reserves by habitat and species
# 6b) success of specific management interventions by habitat and species

#------------------------------ 1a) overall success of AES and nature reserves -------------------------

# code AE.level as factor and relevel to make 'none' the reference level
dat$AE.level <- as.factor(dat$AE.level)
dat$AE.level <- relevel(dat$AE.level, ref="none")

m <- glmer(success ~ AE.level + reserve.desig + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))



#######################################################
#######################################################
#######################################################
#######################################################
#######################################################




