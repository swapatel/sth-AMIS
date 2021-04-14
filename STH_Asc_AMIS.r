###############################################
# Projections using transmission model and geostatistical map of STH
#
#  Using algorithm written by Retkute et al. 2020
#  "Integrating geostatistical maps and transmission models using 
# multiple impotance sampling
#  Modified by SPatel
#
### Requires: 1 map data files, 2 prior files, python code with parameter file, AMIS source, MDA reference file

### for CLUSTER
###############################################
igrp = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
rm(list=setdiff(ls(), "igrp"))
disease <- "Asc"
print(igrp)

library(tmvtnorm)
library(mnormt)
library("mclust")
library(ggplot2)
library(ggpubr)
library(reticulate)
require(gridExtra)
source('~/STH/Ascaris/AMIS_functions.R')
source('~/STH/Ascaris/AMIS_source.R')
source('~/STH/Ascaris/Prior/Ascaris_prior.R')
source('~/STH/Ascaris/create_pyfiles.R')
source('~/STH/Ascaris/qa_ECDF_funcs.R')
plot_folder <- "output/qual_assess_ecdf"

# load map and scen ref data (cluster)
map_data <- read.csv('~/STH/Ascaris/mapdata_ascaris_prev2018.csv')
scen_grp_ref <- read.csv('~/STH/Ascaris/scen_grp_ref_ascaris.csv')

prior_file <- '~/STH/Ascaris/Prior/RawData.csv' # raw data for prior

# pick out IUs
idx_for_grp <- which(scen_grp_ref$group == igrp)
IUs_in_grp <- scen_grp_ref[idx_for_grp, 'IU_ID2']
iscen <- unique(scen_grp_ref[idx_for_grp, "scen"])
if(length(iscen)>1) print("Error in scenario")
idx_for_map <- which(map_data$IU_ID %in% IUs_in_grp)
map_data_for_grp <- map_data[idx_for_map, ]
  
# identify folders for output, need to premake folders and files
folder <- sprintf("output/scen%g/", iscen)
python_file <- sprintf("main_scen%g_grp%g.py", iscen, igrp)
run_py_file <- "run_model.py" 
run_py_file2 <- 'rerun_model.py' # for resampling
# make sure this is consistent with main.py
prevalence_output <- sprintf("output/OutputPrevKKSAC_scen%g_grp%g.csv", iscen, igrp)
prevalence_output2 <- sprintf("output/OutputPrevMHISAC_scen%g_grp%g.csv", iscen, igrp) 
inputRk <- sprintf("output/InputRk_scen%g_grp%g.csv", iscen, igrp) 
create_main_pyfile(python_file, iscen, igrp)

source_python(python_file)

############## AMIS and MAP parameters ############
n.IUs <- length(IUs_in_grp)
n.map.sampl <- 1000
ESS.R<-200 # Desired effective sample 
delta<-.05 # delta value (width for the Radon-Nikodym derivative) 
n.param<-2
maxT<-100 # max number of iterations
NN<-100  # Number of parameter sets in each iteration
N<-rep(NN,maxT)

mean.prev <- sapply(1:n.IUs, function(x) mean(as.numeric(map_data_for_grp[x, 2:(n.map.sampl+1)])))

median.prev <- sapply(1:n.IUs, function(x) median(as.numeric(map_data_for_grp[x, 2:(n.map.sampl+1)])))

# Set distribution for proposal: Student's t distribution
proposal=mvtComp(df=3); mixture=mclustMix(); 
dprop <- proposal$d
rprop <- proposal$r

# Set prior distribution: 
d <- read.csv(prior_file)   # raw data for making prior

param<-matrix(NA, ncol=n.param+2, nrow=sum(N))  # Matrix for seed, parameter values, prevalence 
Sigma <- list(NA, 10*maxT)
Mean<-list(NA, 10*maxT)
PP<-list(NA,maxT)
GG<-list(NA,maxT)


###################################################################	
#          Iteration 1. 
####################################################################	
t<-1  # Iteration 
tmp<-rprop0(N[t])    #N[t] random draws of parameters from prior
x <- tmp[[1]]  # R0
y <- tmp[[2]]  # k
seed <- c(1:N[1])
input_params <- cbind(seed, x, y)  
write.csv(input_params, file=inputRk, row.names=FALSE)

print(Sys.time())
source_python(run_py_file)
print(Sys.time())

res <- read.csv(prevalence_output)
ans <- res[,dim(res)[2]]
n.samples <- sum(N[1:(t)])

param[1:n.samples,1]<-seed
param[1:n.samples,2]<-x
param[1:n.samples,3]<-y
param[1:n.samples,4]<-ans

weight_matrix <- matrix(NA, n.samples, n.IUs)
first_weight <- rep(1, N[1])
sample_prev <- param[1:n.samples,4]
for(i in 1:n.IUs){
  data_prev <- map_data_for_grp[i, 2:(n.map.sampl+1)]
  second_weight_for_IU <- get_second_weight(sample_prev, data_prev, delta, first_weight)
  weight_for_IU <- second_weight_for_IU*first_weight
  if(sum(weight_for_IU)>0) weight_for_IU<-weight_for_IU/sum(weight_for_IU)
  weight_matrix[1:n.samples, i] <-  weight_for_IU
}

ess <- calculate_ESS(weight_matrix)
ESS<-matrix(ess, nrow=1, ncol=n.IUs)

###################################################################	
#          Iteration 2+
####################################################################	
stop<-0
while(stop==0){
  
  t<-t+1
  cat(c("Iteration: ", t,", min(ESS): ", min(ess),"\n"))
  
  wh<-which(ess>=ESS.R)
  W1<-t(weight_matrix); W1[wh,]<-0  # check
  
  w1<- c(colSums(W1))  # check
  
  J<-sample(1:sum(N[1:(t-1)]), NN, prob= w1, replace=T)
  xx<-param[J,2:3] 
  clustMix <- mixture(xx)
  
  G <- clustMix$G
  cluster <- clustMix$cluster
  
  ### Components of the mixture
  ppt <- clustMix$alpha
  muHatt <- clustMix$muHat
  varHatt <- clustMix$SigmaHat
  GG[[t-1]]<-G
  G1<-0; G2<-G
  if(t>2) {
    G1<-sum(sapply(1:(t-2), function(a) GG[[a]]))
    G2<-sum(sapply(1:(t-1), function(a) GG[[a]]))
  }
  for(i in 1:G){
    Sigma[[i+G1]] <- varHatt[,,i]
    Mean[[i+G1]] <- muHatt[i,]	
    PP[[i+G1]]<-ppt[i]   ### scale by number of points
  }
  
  ### Sample new from the mixture...
  ans<-c(); x<-c(); y<-c()
  print("start sampling")
  print(Sys.time())
  while(length(x)<N[t]){
    compo <- sample(1:G,1,prob=ppt) 
    x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
    new.param<-as.numeric(x1)
    if(dprop0(new.param[1],new.param[2])>0 & new.param[1]>=1){
      x<-c(x, new.param[1])
      y<-c(y, new.param[2])
    }
  }
  print("done sampling")
  print(Sys.time())
  
  seed <- c((max(seed)+1): (max(seed)+N[t]))
  input_params <- cbind(seed, x, y)  
  write.csv(input_params, file=inputRk, row.names=FALSE)
  
  source_python(run_py_file) # model outputs to file
  print(Sys.time())
  res <- read.csv(prevalence_output)
  ans <- res[,dim(res)[2]]
  n.samples <- sum(N[1:t])
  
  param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),1]<-seed
  param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),2]<-x
  param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),3]<-y
  param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),4]<-ans
  
  prop.val <- sapply(1:n.samples,function(b)  sum(sapply(1:G2, function(g) PP[[g]] * dprop(param[b,2:3],mu= Mean[[g]], Sig=Sigma[[g]]))) + dprop0(param[b,2], param[b,3])) # assumes all iterations have same number of samples
  
  first_weight <- sapply(1:n.samples, function(b) dprop0(param[b,2], param[b,3])/prop.val[b])   # prior/proposal
  
  # get second weight
  weight_matrix <- matrix(NA, n.samples, n.IUs)
  for(i in 1:n.IUs){
    sample_prev <- param[1:n.samples,4]
    data_prev <- map_data_for_grp[i, 2:(n.map.sampl+1)]
    second_weight_for_IU <- get_second_weight(sample_prev, data_prev, delta, first_weight)
    weight_for_IU <- second_weight_for_IU*first_weight
    if(sum(weight_for_IU)>0) weight_for_IU<-weight_for_IU/sum(weight_for_IU)
    weight_matrix[1:n.samples, i] <-  weight_for_IU
  }

  ess <- calculate_ESS(weight_matrix)
  cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))
  ESS<-rbind(ESS, as.numeric(ess))
  weight_all <- c(rowSums(weight_matrix))
  if(min(ess)>=ESS.R) stop<-1
  if(t>= maxT) stop<-1
}

print(t)
param_weight <- cbind(param[1:n.samples,], weight_all, weight_matrix)
write.csv(param_weight,file=paste0("paramWW_scen", iscen, "_grp", igrp, ".csv"), row.names = FALSE)

## quality assurance plots
ecdf_df <- qa_plots(param_weight, map_data_for_grp, plot_folder, iscen, igrp)
write.csv(ecdf_df, file=paste0(plot_folder, "/ecdf_scen", iscen, "_grp", igrp, ".csv"), row.names = FALSE)

write.csv(IUs_in_grp, file=paste0(folder, "IUs_scen", iscen, "_grp", igrp, ".csv"))

###################END OF AMIS 

print("checkpoint 1")

#############################################
### RESAMPLING ###########################
###########################################
PickleFile = sprintf('output/OutputPickle_scen%g_grp%g.p', iscen, igrp)
MDAFilePath = sprintf('MDAfiles/STH_MDA_scen%g.csv', iscen)

# if resampling done separately, need to read in files
# paramFile = paste0("paramWW_scen", iscen, "_grp", igrp, ".csv")
# IU_file <- paste0(folder, "IUs_scen", iscen, "_grp", igrp, ".csv")
# param_weight <- read.csv(paramFile)
# IUs_in_grp <- read.csv(IU_file)
# IUs_in_grp <- IUs_in_grp[,2]
# n.samples <- dim(param_weight)[1]
# n.IUs <- length(IUs_in_grp)

n.fill <- 5 # number of rows before weights start
desired_samples <- 200
first_ten <- floor(seq(1, desired_samples, length.out=10))

for(i in 1:n.IUs){
  print(i)
  w <- param_weight[,(n.fill+i)]
  resamples <- sample(1:n.samples, desired_samples, replace=TRUE, prob=w)
  print(length(unique(resamples)))
  
  # sort based on prevalence
  reparamWW <- data.frame(param_weight[resamples, 1:(n.fill-1)]) #4th is prev
  names(reparamWW)[n.fill-1] <- "prev"
  temp <- order(reparamWW$prev)
  order_params <- c(temp[first_ten], c(1:desired_samples)[-first_ten])  
  # end sort 
  
  print("checkpoint 2")
  reparamWW <- reparamWW[order_params, 1:3]
  write.csv(reparamWW, file=paste0(folder, "Input_Rk_", disease, "_", IUs_in_grp[i], ".csv"), row.names=FALSE) #outputs to scen folder
  
  ##### Simulate in python ######
  write.csv(reparamWW, file=inputRk, row.names=FALSE) 
  source_python(run_py_file2)
  file.rename(from=prevalence_output, to=paste0(folder, "PrevKKSAC", disease, "_", IUs_in_grp[i], ".csv"))
  file.rename(from=prevalence_output2, to=paste0(folder, "PrevMHISAC", disease, "_", IUs_in_grp[i], ".csv"))
  
  file.rename(from=PickleFile, to=paste0(folder, disease, "_", IUs_in_grp[i], ".p"))
  file.copy(from=MDAFilePath, to=paste0(folder, "STH_MDA_", IUs_in_grp[i], ".csv"))
  ##### END #####################
}