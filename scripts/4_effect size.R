#################################################################
#
#     Step 4: EU meadow birds meta-analysis - models evaluating effect size
#
#################################################################

# Samantha Franks
# 27 April 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","raster","sp","rgeos","rgdal","lme4","car","blme","tidyr","nlme")

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
dat1 <- subset(dat0, select=c("reference","lit.type","country","study.length","habitat","species","overall.metric","metric","specific.metric","stan.metric","sample.size","analysis2","success","stan.calc","stan.metric.before","stan.metric.after","stan.effect.size","sig","effect.dir",mgmtvars))

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

#---------------------   Effect size in relation to management intervention --------------------


mod <- list()
moddat <- list()

for (i in 1:length(metrics)) {
  
  print(metrics[i])
  mdat <- subset(dat, new.stan.metric==metrics[i])
  
  # check for sample sizes of different management treatment groups, remove small ones
  checksample <- table(mdat$mgmt)
  removelevels <- names(checksample)[checksample < 5]
  mdat <- subset(mdat, !grepl(paste(removelevels,collapse="|"), mdat$mgmt))
  mdat <- droplevels(mdat)
  moddat[[i]] <- mdat
  
  print(nrow(mdat))
  print(table(mdat$mgmt))
  
  mod[[i]] <- lme(stan.effect.size ~ mgmt, random = ~1|reference, data=mdat)
  print(summary(mod[[i]]))
  
}
m <- lme(stan.effect)

# for making pairwise comparisons between the different treatment groups, use glht() in multcomp package - nice description of how given here: http://mindingthebrain.blogspot.co.uk/2013/04/multiple-pairwise-comparisons-for.html

