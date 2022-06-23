####################################################################################################################################################################################################
#Code for running complex survey simulations from: "An Augmented Likelihood Approach for the Discrete Proportional Hazards Model Using Auxiliary and Validated Outcome Data - with Application to HCHS/SOL Study"
####################################################################################################################################################################################################

#We need to load in our own functions for the log-likelihood/gradient and some functions from Gu et al. 2015 that help our calculations 
#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('U:/Paper 2/Code/RcppFunc.cpp')
source('U:/Paper 2/Code/Paper 2 Likelihood FINAL.R')

options(survey.adjust.domain.lonely=TRUE)
options(survey.lonely.psu="adjust")

library(pracma)
library(tidyr)
library(MASS)
library(reshape2)
library(numDeriv)

#Function for getting rid of any 1's after first true positive
#Get rid of any 1's after first true positive
after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}


#Set sensitivity and specificity
sensitivity<-0.80
specificity<-0.90

#Setting to vary: proportion of missingness in gold standard
prop_m<-0

#What (average) sample size, censoring rate, and covariate distirbutions do we want to consider
#meanN<-1000
meanN<-10000

#myCR<-0.9
#myCR<-0.7
myCR<-0.5

#covardist<-"Normal"
covardist<-"Gamma"

#Starting value for beta
beta1_0<-0.5

### Assign population name
popname = paste0('TargetPopulationData',covardist,'_CR',myCR,'_N',meanN)
sampname = 'SampleData'

#Number of samples drawn 
idx=1:430
S=1000

### Load in data
files = paste0('U:/',popname,'/',sampname,"_",str_pad(1:S,nchar(S),pad=0),".RData")

#Create an empty matrix for saving values 
ansmat<-NULL

#Set seed for reproducibility
set.seed(1020)

for(i in idx){  
  load(file=files[i])
  samp <- as.data.table(samp)
  
  N<-dim(samp)[1]
  
  #calculate mean CR for this sample
  truecens<-mean(samp$true_result_4==0)
  
  #put samp in long form
  samp_long1<-reshape(data = as.data.frame(samp),
                      idvar = "ID",
                      varying = list(true_result=c("true_result_1","true_result_2","true_result_3","true_result_4"),result=c("result_1","result_2","result_3","result_4")),
                      direction="long",
                      v.names = c("true_result","result"),
                      sep="_")
  
  #order long dataset by each subject's ID
  samp_long <- samp_long1[order(samp_long1$ID),]
  
  #Save gold standard (GS) data - just saving off each individual's ID and GS binary indicator (delta) which is GS_vis4 here
  GS_data<-samp_long[!duplicated(samp_long$ID),c("ID","GS_vis4")]
  
  #Save vector of weights for this dataset
  weights<-as.numeric(unlist(samp_long[!duplicated(samp_long$ID),c("bghhsub_s2")]))
  
  #Delete all positives after first one for error-prone variable, Y* (called "result" here)
  keep<-unlist(tapply(samp_long$result,samp_long$ID,after_first_pos))
  datafinal_1<-samp_long[keep,]
  
  #Consider that a certain proportion of gold standard outcomes are randomly missing (MCAR)
  mcar<-runif(N,0,1)
  #create an indicator variable called GS_vis4_mis indicting whether each subject missing GS or not
  GS_data$GS_vis4_mis<-ifelse(mcar<prop_m,NA,GS_data$GS_vis4)
  NMiss_GS<-table(GS_data$GS_vis4_mis,exclude=FALSE)[3]
  
  #Create final data sets for the proposed method analysis and gold standard analysis
  datafinal<-merge(datafinal_1,GS_data,by="ID")
  datafinal_GS<-merge(samp_long,GS_data[,c("ID","GS_vis4_mis")],by="ID")
  
  #Create formula: just a single covariate x in the outcome model
  formula=result~myX
  
  #Create ID, testtime, result variables for analysis
  id <- eval(substitute(datafinal$ID), datafinal, parent.frame())
  time <- eval(substitute(datafinal$time), datafinal, parent.frame())
  result <- eval(substitute(result), datafinal, parent.frame())
  ord <- order(id, time)
  
  #Make sure data is ordered properly before we calculate the sum of the likelihood ocmponents 
  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    datafinal <- datafinal[ord, ]}
  utime <- sort(unique(time))
  timen0 <- (time != 0)
  
  #Create D and C matrices, as described in the methods section of the main manuscript (this is where we use functions from Gu et al. 2015)
  Dm <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
             specificity, 1)
  Cm <- cmat(id[timen0], time[timen0], result[timen0], sensitivity,
             specificity, 1)
  
  #Save off J, number of visits
  J <- ncol(Dm) - 1
  nsub <- nrow(Dm)
  
  #Create our xmatrix
  Xmat <- model.matrix(formula, data = datafinal)[, -1, drop = F]
  beta.nm <- colnames(Xmat)
  nbeta <- ncol(Xmat)
  uid <- getrids(id, nsub)
  
  #Reduce xmatrix to one row per person
  Xmat <- Xmat[uid, , drop = F]
  
  #Assign gold standard indicator, GSdelta, the visit time for the gold standard, and the weights 
  #For most of our simulations, we assume validation time is time 4, the last auxiliary data visit time
  GSdelta <- datafinal[uid,"GS_vis4_mis"]
  GSVis<-rep(4,nsub)
  
  #Create interval censored data
  IC_data<-datafinal_GS[!duplicated(datafinal_GS$ID),c("result","true_result","myX","GS_vis4","GS_vis4_mis","BGid", "strat", "bghhsub_s2")]
  #Now create final dataset for standard IC analysis, omitting anyone with missing GS
  IC_GS_datafinal<-IC_data[complete.cases(IC_data$GS_vis4_mis),]
  
  #Now create data with gold standard observed every 4 years (no error-prone auxiliary)
  IC_data4year<-datafinal_GS[,c("ID","time","true_result","myX","GS_vis4_mis","BGid", "strat", "bghhsub_s2")]
  #Get rid of anyone in gold standard, observe every year dataset who is missing visit 4
  IC_data4year_complete<-IC_data4year[complete.cases(IC_data4year$GS_vis4_mis),]
  
  #Only keep positive gold standard after first positive for case where observed all four years
  keepafterpos<-unlist(tapply(IC_data4year_complete$true_result,IC_data4year_complete$ID,after_first_pos))
  datafinal_4vis<-IC_data4year_complete[keepafterpos,]
  
  #Now create survey designs using survey package.
  #Create 3 designs based on (1) long data, error-prone result, (2) long data, true result, (3) one row per subject, true result (GS)
  #samp_design_long_E = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=datafinal)
  samp_design_long_T = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=datafinal_4vis)
  samp_design_reg = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=IC_GS_datafinal)
  samp_design_reg_complete = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=IC_data)
  
  #Fit model: grouped time survival approach using  gold standard outcome observed at all 4 years
  fit_truth4vis<-svyglm(true_result~time+myX,family=quasibinomial(link="cloglog"),data=datafinal_4vis,design=samp_design_long_T)
  fitsum_truth4Vis<-summary(fit_truth4vis)
  beta1est_truth4Vis<-fitsum_truth4Vis$coefficients["myX",1]
  beta1se_truth4Vis<-fitsum_truth4Vis$coefficients["myX",2]
  
  #Fit model: interval standard, gold standard where gold standard only observed at year 4
  fit_truth4Year<-svyglm(GS_vis4~myX,family=quasibinomial(link="cloglog"),data=IC_GS_datafinal,design=samp_design_reg)
  fitsum_truth4Year<-summary(fit_truth4Year)
  beta1est_truth4Year<-fitsum_truth4Year$coefficients["myX",1]
  beta1se_truth4Year<-fitsum_truth4Year$coefficients["myX",2]
  
  #Create initial parameters for optim for proposed method
  lami <- log(-log(seq(1, 0.1, length.out = J + 1)[-1]))
  lami <- c(lami[1], diff(lami))
  tosurv <- function(x) exp(-exp(cumsum(x)))
  lowlam <- c(-Inf, rep(0, J - 1))
  lowerLBFGS <- c(rep(-Inf, nbeta),lowlam)
  upperLBFGS <- c(Inf,Inf,Inf,Inf,Inf)
  parmi <- c(beta1_0,lami)

  #Begin optimization for WEIGHTED log likelihood. Getting convergence issues in optim with L-BFGS-B so try just BFGS.
  proposed_fit<-nlminb(start=parmi, objective=log_like_proposed, gradient= gradient_proposed,
                           lower = lowerLBFGS, upper = upperLBFGS,nsub=nsub,J=J,nbeta=nbeta,Dm=Dm,
                            Cm=Cm,Xmat=Xmat,GSdelta=GSdelta,GSVis=GSVis,weights=weights,
                           purpose="SUM")
  result_var_w<-solve(pracma::hessian(log_like_proposed,proposed_fit$par,nsub=nsub,J=J,nbeta=nbeta,Dm=Dm,
                                      Cm=Cm,Xmat=Xmat,GSdelta=GSdelta,GSVis=GSVis,weights=weights,
                                      purpose="SUM"))
  
  beta1est_aux<-  proposed_fit$par[1]
  beta1se_aux<-sqrt(result_var_w[1,1])
  
  #Calculate Sandwich Variance: Proposed Estimator
  #Be sure to provided matrix of unweighted estimating equation contributions (here let weights =1 )
  U_prop1<-gradient_proposed(proposed_fit$par,nsub=nsub,J=J,nbeta=nbeta,Dm=Dm,
                             Cm=Cm,Xmat=Xmat,GSdelta=GSdelta,GSVis=GSVis,weights=rep(1,nsub),
                             purpose="INDIVIDUAL")
  infl1<- U_prop1%*%result_var_w
  mySandVar<- vcov(svytotal(infl1,samp_design_reg_complete))
  sandVar<-  sqrt(diag(mySandVar))[1]
  
  #Combine all results
  ansmat<-rbind(ansmat,c(beta1est_aux,beta1se_aux,sandVar,
                         beta1est_truth4Year,beta1se_truth4Year,beta1est_truth4Vis,beta1se_truth4Vis,truecens))
  
  print(i)
}

ansmat<-data.frame(ansmat)
names(ansmat)<-c("beta1est_aux","beta1se_aux_naive","beta1se_sand",
                 "beta1est_truth4Year","beta1se_truth4Year",
                 "beta1est_truthAll4Vis","beta1se_truthAll4Vis","CR")

#Function for assessing performance
Data_output_singlecov_median <- function(betaest,betase,betatrue){
  
  ASE<-median(betase)
  ESE<-sd(betaest)
  MAD<-mad(betaest)
  
  #Median percent bias corrected = (estimated-target)/target
  bias<-(betaest-betatrue)/betatrue
  median_percent_bias<-median(bias)*100
  
  #Coverage Probability Calculation
  beta_1_l95<-betaest-qnorm(0.975)*(betase)
  beta_1_r95<-betaest+qnorm(0.975)*(betase)
  
  beta_1_ci<-as.data.frame(cbind(beta_1_l95,beta_1_r95))
  
  # check the proportion of intervals containing the parameter
  CP1<-mean(apply(beta_1_ci, 1, findInterval, x = betatrue) == 1)
  
  results<-cbind(median_percent_bias,ASE,MAD,CP1)
  
  return(results)
}

#Proposed
Data_output_singlecov_median(ansmat$beta1est_aux,ansmat$beta1se_sand,log(1.5))

#Gold Standard, Interval Censored
Data_output_singlecov_median(ansmat$beta1est_truth4Year,ansmat$beta1se_truth4Year,log(1.5))

#Gold Standard all four years, Grouped Survival
Data_output_singlecov_median(ansmat$beta1est_truthAll4Vis,ansmat$beta1se_truthAll4Vis,log(1.5))
