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

mgmtvars <- c("AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# subset dataset for analysis to desired columns only
dat <- subset(dat0, select=c("reference","country","study.length","habitat","species","overall.metric","sample.size","analysis2","success",mgmtvars))


#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}

#=================================  ANALYSIS  ===============================

# 0a) success of individual management types overall
# 0b) success of individual management types by species
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

### Nuisance explanatory variables and random effects
# ‘nuisance explanatory variables: study duration (continuous), sample size (categorical: small, medium, large), multivariate/univariate (categorical)
# Random effects: country, study id - habitat will be included later on

#------------------------------ 0a) success of individual management types -------------------------

m <- list()

print(Sys.time())

for (i in 1:length(mgmtvars)) {
  
#   mdat <- dat[dat[,mgmtvars[i]]!="none",]
#   mdat <- droplevels(mdat)
  
  mdat <- dat
  
  modform <- as.formula(paste("success ~ ", mgmtvars[i], " + study.length + analysis2 + sample.size + (1|reference) + (1|country)", sep=""))
  m[[i]] <- glmer(modform, data=mdat, family=binomial, control=glmerControl(optimizer="bobyqa"))
  
}

print(Sys.time())


#------------------------------ 1a) overall success of AES and nature reserves -------------------------
m <- list()

# subset to include data where AES or reserve management was attempted
subdat <- subset(dat, AE.level!="none" | reserve.desig!="none")
subdat <- droplevels(subdat)

print(Sys.time())

m[[1]] <- glmer(success ~ AE.level + reserve.desig + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())

#------------------------------ 1b) success of specific management interventions -------------------------

print(Sys.time())

m[[2]] <- glmer(success ~ mowing + grazing + fertpest + nest.protect + predator.control + water + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())

#------------------------------ 2a) success of AES and nature reserves by species -------------------------

# subset dataset to remove dunline and ruff, not enough data for these species
subdat <- subset(dat, species!="dunlin" & species!="ruff")
subdat <- droplevels(subdat)

print(Sys.time())

m[[3]] <- glmer(success ~ AE.level*species + reserve.desig*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=subdat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())


#------------------------------ 2b) success of specific management interventions by species -------------------------

# subset dataset to remove dunline and ruff, not enough data for these species
subdat <- subset(dat, species!="dunlin" & species!="ruff")
subdat <- droplevels(subdat)

print(Sys.time())

m[[3]] <- glmer(success ~ AE.level*species + reserve.desig*species + study.length + analysis2 + sample.size + (1|reference) + (1|country), data=subdat, family=binomial, control=glmerControl(optimizer="bobyqa"))

print(Sys.time())


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################




