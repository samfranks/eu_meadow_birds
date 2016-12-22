#################################################################
#
#     Step 3: EU meadow birds meta-analysis - NUISANCE VARIABLE ANALYSIS
#
#################################################################

# Samantha Franks
# 11 March 2016
# 22 Dec 2016


# ----------- NOTES ABOUT ANALYSIS -----------------

# significance is given by Wald t tests (default for summary.glmer()) - these are not the greatest measure of significance according to the GLMM wiki, so use the results from the LRT tests instead which compare models with and without each dropped covariate
# the only significant effect of a nuisance variable is literature type (when 'score' is not included as a variable)
# primary literature study is more likely to be unsuccessful than successful
# controlling for interaction between lit.type and score removes significant effect of lit.type, suggesting that variance in the data explained by literature type could be accounted for by the quality of the analysis
# including score only, it's not significant
# including a variable "biased metric" which includes the abundance/occupancy problem suggests no different in success between unbiased vs biased metrics

# proportion of scores for each literature type
#             good   medium     poor
# grey    0.224299 0.485981 0.289720
# primary 0.405350 0.465021 0.129630



# =================================  SET LOGIC STATEMENTS  ====================



# =================================  LOAD PACKAGES =================================

list.of.packages <- c("MASS","reshape","dplyr","lme4","car","blme","tidyr","nlme")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only=TRUE)


# =================================  LOAD FUNCTIONS =================================



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



# ================================== Test effect of nuisance variables on success for the full dataset ===========================

vars <- c("study.length","sample.size","analysis2","lit.type","score","biased.metric")

m.global <- glmer(success ~ study.length + sample.size + analysis2 + lit.type + score + biased.metric + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))

m <- list()
for (i in 1:length(vars)) {
  m[[i]] <- drop1(m.global, scope = as.formula(paste("~", vars[i], sep="")), test="Chisq")
}


# m.nui0 <- glmer(success ~ study.length + sample.size + analysis2 + lit.type + score + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# summary(m.nui0)
# 
# m.nui1 <- glmer(success ~ study.length + sample.size + analysis2 + lit.type*score + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# summary(m.nui1)
# 
# m.nui2 <- glmer(success ~ study.length + sample.size + analysis2 + score + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# summary(m.nui2)
# 
# m.nui3 <- glmer(success ~ study.length + sample.size + analysis2 + lit.type*score + biased.metric + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# summary(m.nui3)
# 
# m.nui4 <- glmer(success ~ study.length + sample.size + analysis2 + biased.metric + (1|reference), data=dat, family=binomial, control=glmerControl(optimizer="bobyqa"))
# summary(m.nui4)


# --------- OUTPUT MODEL RESULTS --------

setwd(outputwd)
sink(paste("model output_nuisance variables.txt", sep=" "))

cat("\n########==========  Nuisance variables - Likelihood ratio tests for each variable in global model ==========########\n", sep="\n")
print(summary(m.global))

cat("\n###---  Likelihood Ratio Tests ---###\n", sep="\n")
print(m)

sink()



# ---------- OUTPUT TABLE of nuisance variable LRT results (dplyr/piping way) -------

lapply(m, function(x) {data.frame(df=x[,"Df"][2], LRT=x[,"LRT"][2], pval=x[,"Pr(Chi)"][2]) %>% return
}) %>%
  do.call(rbind, .) %>%
  data.frame(vars) %>%
  write.csv(paste(outputwd, "model output_nuisance variables table.csv", sep="/"), row.names=FALSE)


# # ---------- OUTPUT TABLE of nuisance variable LRT results (old school way) -------
# 
# # extract df, chi-square, and p-values from each dropped term into a data.frame
# # output data.frame to a csv file
# # old way
# mout.list <- lapply(m, function(x) {
#   modout <-data.frame(df=x[,"Df"][2], LRT=x[,"LRT"][2], pval=x[,"Pr(Chi)"][2])
#   return(modout)
# })
# 
# mout.all <- do.call(rbind, mout.list)
# mout.all <- data.frame(mout.all, vars) # append name of variable dropped to its relevant output
# write.csv(mout.all, paste(outputwd, "model output_nuisance variables table.csv", sep="/"), row.names=FALSE)






