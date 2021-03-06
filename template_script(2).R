###### Set working directory to the desired location
setwd("C:/Users/Marc/Desktop/Community Ecology Research")
setwd("X:/My desktop/Community Ecology Research")
##################
rm(list=ls())
### Clear environment
# rm(list = ls()) # /!\ Make sure you saved everything important before
#### install packages
install.packages('picante')
install.packages('DDD')
install.packages('pspline')
install.packages('qpcR')
#### load packages
library('picante')
library('DDD')
#library('RPANDA') # RPANDA function fit_env() does not seem to work
library('pspline')
library('qpcR')
##################
# source RPANDA functions instead
source("fit_bd_positive.R")
source("fit_env_bd.R")
source("likelihood_bd.R")
source("Phi.R")
source("Psi.R")
source("integrate.R")

##### load datasets
# Temperature data
env_data<-read.table('PhanerozoicTemperature.txt', header = T)
df<-100 # degrees of freedom for the smooth spline function in fit_env_bd()
# Phylogenetic data
load('TempFamilies_dataset.Rdata')
names<-names(phylogenies)
###################

#### set parameters for ML functions
# fit_bd() & fit_env_bd()
cond_RPANDA <-'crown'
meth<-"Nelder-Mead" # optimisation method
# DD_ML()
ddmodel <- 1 # linear dependence in speciation rate with K
tol <- rep(1E-6,3) # tolerance
methode = 'analytical' # solving method
cond_DDD <- 1 # conditioning

#####################

### set output datasets
res<-c() # empty list to append through the loop
# default results if optimisation fails
fit_error<-list(
  'model' = 'birth death',
  'LH' = -1,
  'aicc'= -1,
  'lamb_par'=-1,
  'K_par'=-1
)
#######################

### Iterate through phlyogenies to fit all 4 models
for(i in 1:length(phylogenies[1:3])){
  
  # Pick current phylogeny in the list
  phylo<-phylogenies[[i]]$tree
  tot_time<-max(node.age(phylo)$ages) # extract age
  # Extract number of taxa and sampling fraction
  nobs<-Ntip(phylo)
  ntotal<- as.numeric(phylogenies[[i]]$totalsp)
  f <- nobs / ntotal
  
  # BCST (Pure Birth)
  print("BCST")
  f.lamb<-function(t,y){y[1]}
  f.mu<-function(t,y){0}
  lamb_par<-c(0.1)
  mu_par<-c()
  cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  fit_BCST<-fit_error
  try(fit_BCST<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,
                         cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,
                       expo.mu=expo.mu,fix.mu=fix.mu,cond=cond_RPANDA))
  print(fit_BCST)
  
  # BTimeVar LIN
  print("BTimeVar LIN")
  f.lamb<-function(t,y){y[1]+y[2]*t}
  f.mu<-function(t,y){0}
  lamb_par<-c(0.1,0.01)
  mu_par<-c()
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  fit_BTimeVar_LIN<-fit_error
  try(fit_BTimeVar_LIN<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,
                               cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,
                               expo.mu=expo.mu,fix.mu=fix.mu,cond=cond_RPANDA))
  print(fit_BTimeVar_LIN)
  
  # BTempVar LIN
  print("BTempVar LIN")
  f.lamb<-function(t,x,y){y[1]+y[2]*x}
  f.mu<-function(t,x,y){0}
  lamb_par<-c(0.1,0.01)
  mu_par<-c()
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  fit_BTempVar_LIN<-fit_error
  try(fit_BTempVar_LIN<-fit_env_bd(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,
                                   cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,
                                   expo.mu=expo.mu,fix.mu=fix.mu,cond=cond_RPANDA))
  print(fit_BTempVar_LIN)
  
  # BDivDep LIN
  brts <- as.numeric(branching.times(phylo))
  K_ini <- ntotal # assumes the clade is currently at equilibrium diversity
  # DD_ML() parameters
  initparsopt <- c(0.1,K_ini)+1E-6
  idparsopt <-c(1,3)
  idparsfix <-2
  parsfix<-0+1E-6
  p <- length(initparsopt)
  missnumspec <- ntotal - nobs
  
  fit_BDivDep_LIN<-fit_error
  DD_res<- try(dd_ML(brts,initparsopt = initparsopt, idparsopt = idparsopt,
                     idparsfix = idparsfix, parsfix = parsfix,
                     ddmodel = ddmodel,missnumspec = missnumspec,
                     cond = cond_DDD,tol = tol, methode = methode))
  # convert DD_ML() output to RPANDA output format
  fit_BDivDep_LIN$LH <- DD_res$loglik
  fit_BDivDep_LIN$aicc <- 2*(-DD_res$loglik) + 2*p + (2*p*(p+1)) / (nobs-p-1)
  fit_BDivDep_LIN$lamb_par <- DD_res$lambda
  fit_BDivDep_LIN$K_par <- DD_res$K
  
  print(fit_BDivDep_LIN)
  
  ###### Collect results ########
  ### Set results matrix
  results<-matrix(NA,4,8)
  colnames(results)<-c("Models","Parameters","logL","AICc","AICw","Lambda","Alpha","K")
  results[,1]<-c("BCST","BTimeVar_LIN","BTempVar_LIN","BDivDep_LIN")
  
  ### Fill results matrix
  
  #Nb of parameters
  results[1,2]<-1
  results[2:4,2]<-2
  
  #logL
  results[1,3]<-fit_BCST$LH
  results[2,3]<-fit_BTimeVar_LIN$LH
  results[3,3]<-fit_BTempVar_LIN$LH
  results[4,3]<-fit_BDivDep_LIN$LH
  
  #AICc
  results[1,4]<-fit_BCST$aicc
  results[2,4]<-fit_BTimeVar_LIN$aicc
  results[3,4]<-fit_BTempVar_LIN$aicc
  results[4,4]<-fit_BDivDep_LIN$aicc
  
  #AICw
  allAIC<-c(results[1:4,4])
  AICw<-akaike.weights(as.numeric(allAIC))
  results[,5]<-round(AICw$weights, digits = 2)
  
  #Lambda0
  results[1,6]<-fit_BCST$lamb_par[1]
  results[2,6]<-fit_BTimeVar_LIN$lamb_par[1]
  results[3,6]<-fit_BTempVar_LIN$lamb_par[1]
  results[4,6]<-fit_BDivDep_LIN$lamb_par[1]
  
  #Alpha
  results[2,7]<-fit_BTimeVar_LIN$lamb_par[2]
  results[3,7]<-fit_BTempVar_LIN$lamb_par[2]

  #K
  results[4,8]<-floor(fit_BDivDep_LIN$K_par)
  
  ###########
  
  # rearrange matrix to put best model on top
  results <- results[order(as.numeric(results[,"AICc"])),]
  print(results)
  
  # Set up final output for the clade
  resi<-list("Clade"=names[i],
             "Clade_age"=tot_time,
             "Taxon_sampling"=nobs,
             "Clade_size"=ntotal,
             "Sampling_fraction"=f,
             "Results"=results
             )
  # Append results to results for all clades
  res<-c(resi,resi)
  
  # Save at R dataset format
  saveRDS(res,file = 'fits_Temp_DD.rds')
  
  #############################
}

