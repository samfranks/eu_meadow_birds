#################################################################
#
#     Step 2: EU meadow birds meta-analysis - models evaluating success
#
#################################################################

# Samantha Franks
# 11 March 2016


#=================================  SET LOGIC STATEMENTS  ====================



#=================================  LOAD PACKAGES =================================



#=================================  LOAD FUNCTIONS =================================



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
# outputwd <- paste(parentwd, "output", sep="/")
workspacewd <- paste(parentwd, "workspaces", sep="/")

options(digits=6)


#=================================  LOAD DATA  ===============================

dat0 <- readRDS(paste(workspacewd, "meadow birds analysis dataset_full.rds", sep="/"))

mgmtvars <- c("AE","AE.level","reserve.desig","mowing","grazing","fertpest","nest.protect","predator.control","water")

# subset dataset for analysis to desired columns only
dat <- subset(dat0, select=c("reference","lit.type","score","country","study.length","habitat","species","overall.metric","metric","sample.size","analysis2","success","failure","outcome",mgmtvars))

### Combine habitat categories ###

# combine habitat categories even further into 3 different categories
# arable: includes arable and mixed arable/pastoral
# pastoral: includes pastoral and mixed pastoral/unenclosed
# unenclosed
dat$newhabitat <- ifelse((dat$habitat %in% c("arable","mixed arable/pastoral")), "arable", ifelse((dat$habitat %in% c("pastoral","mixed pastoral/unenclosed")), "pastoral", "unenclosed"))

### Identify studies where specific intervention is evaluated ###

# identify studies where a specific intervention is evaluated (could be in combination with others, and also with higher-level interventions)
# if any specific measure is used, then dat$spec.int.used=1, otherwise if all specific measures are 'none', then spec.int.used=0
dat$spec.int.used <- with(dat, ifelse(mowing=="none" & grazing=="none" & fertpest=="none" & nest.protect=="none" & predator.control=="none" & water=="none", 0, 1))


### Identify studies where high level interventions are evaluated ###

# identify studies where a high level intervention is evaluated (could be in combination with others, and also with specific interventions)
# if any higher level measure is used, then dat$high.int.used=1, otherwise if all high level measures are 'none', then high.int.used=0
dat$high.int.used <- with(dat, ifelse(AE.level=="none" & reserve.desig=="none", 0, 1))


### Identify abundance change studies ###
# dat$new.metric <- with(dat, ifelse(dat$overall.metric=="abundance change", "abundance change", metric))
dat$new.metric <- with(dat, ifelse(dat$overall.metric=="abundance" | dat$overall.metric=="occupancy", "abundance/occupancy", ifelse(dat$overall.metric=="abundance change" | dat$overall.metric=="occupancy change", "abundance/occupancy change", metric)))

### Identify abundance/occupancy studies ###
dat$biased.metric <- with(dat, ifelse(dat$overall.metric=="abundance" | dat$overall.metric=="occupancy", "Y", "N"))



#------------  Recode management variables as factors for analysis and make 'none' the reference level -----------------

for (i in 1:length(mgmtvars)) {
  dat[,mgmtvars[i]] <- as.factor(dat[,mgmtvars[i]])
  dat[,mgmtvars[i]] <- relevel(dat[,mgmtvars[i]], ref="none")
  print(levels(dat[,mgmtvars[i]]))
}



