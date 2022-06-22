####################################################################################################################################################################################################
#Code for running simulations from: "An Augmented Likelihood Approach for the Discrete Proportional Hazards Model Using Auxiliary and Validated Outcome Data - with Application to HCHS/SOL Study"
####################################################################################################################################################################################################

#We need to load in our own functions for the log-likelihood/gradient and some functions from Gu et al. 2015 that help our calculations 
#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('U:/Paper 2/Code/RcppFunc.cpp')
source('U:/Paper 2/Code/Paper 2 Likelihood FINAL.R')

#If running code using Shell Script in Linux, can loop through different arguments to vary 
# args=(commandArgs(TRUE))
# 
# ## args is now a list of character vectors
# ## Check to see if any arguments were passed.
# if(length(args)==0){
#   print("No arguments supplied")
#   ## supply default values
#   prop_m<-0
#   N<-0
#   myCR<-0
# } else{
#   prop_m<-as.numeric(args[1])
#   N=as.numeric(args[2])
#   myCR=as.numeric(args[3])
# 
# }


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

######################################################################
################Begin data generating mechanism#######################
######################################################################

#Settings to vary: proportion of missingness in gold standard, sample size, censoring rate for latent true event time
prop_m<-0
N<-1000
myCR<-0.9

#Set sensitivity and specificity
sensitivity<-0.80
specificity<-0.90

#Vector of visit times for auxiliary data
testtimes<-c(1,2,3,4)
ntest <- length(testtimes)

#Data is going to be in long form: creates 4 IDs per subject 1-1000
ID <- rep(1:N, each = ntest)
time <- rep(testtimes, times = N)

#Pick baseline distribution for event times
#dist<-"Weib"
dist<-"Exp"

#Pick covariate distribution
#covardist<-"Normal"
covardist<-"Gamma"

#Choose baseline parameter to generate event times based on average censoring rate we want and event time distribution
blambda<-ifelse(dist=="Exp",ifelse(myCR==0.5,0.17,ifelse(myCR==0.7,0.08,ifelse(myCR==0.9,0.023,NA))),
       ifelse(dist=="Weib",ifelse(myCR==0.5,0.39,ifelse(myCR==0.7,0.19,ifelse(myCR==0.9,0.052,NA)))))

#Set "shape" for event time generation
shape<-0.4 #shape=1 implies constant (exponential) baseline hazard

#Set a reasonable value for our true beta_1 and true beta_2
beta_1<-log(1.5)

NSIM<-1000

#Create an empty matrix for saving values 
ansmat<-TrueCensRate<-NULL

#Set seed for reproducibility
set.seed(1020)

for(iter in 1:NSIM){
  time <- rep(testtimes, times = N)

  x_1<-if (covardist=="Gamma") {
    rgamma(N, shape=1/5, scale = 1)
  }else if (covardist=="Normal") {
      rnorm(N, 1/5, 1)
    }else (NA)


  #Uses exponential or Weibull distribution to generate  lambda, then time to event of interest (event time = ET)
  lambda<-if (dist=="Exp") {blambda*exp(x_1*beta_1)} else {(blambda)^(-1/shape)/exp((x_1*beta_1)/shape)}
  ET<-if (dist=="Exp") {rexp(N, lambda)} else {rweibull(N, shape, lambda)}

  #put covariate in long form
  x_1<- x_1[ID]

  #Repeats event time for each person for total number of visits (data in long form)
  ET <- ET[ID]

  #Indicator for whether visit time >= actual event time
  occur <- time >= ET
  true_result<-as.numeric(occur)

  #If visit time >= actual event time, then we flip sensitivity coin.
  #Otherwise, flip specificity coin.
  probs <- ifelse(occur, sensitivity, 1 - specificity)

  #generates (N)*(#visit times) points of 0 or 1, with probability vector "probs" above
  #Having a probability vector means each number in that vector is cycled through as the probability at that particular coin toss
  result <- rbinom(length(occur), 1, probs)

  #Create our data set
  data1<-data.frame(ID,time,result,true_result,x_1)

  #Now delete any positives after first positive self-report/auxiliary data outcome
  keep<-unlist(tapply(data1$result,data1$ID,after_first_pos))
  datafinal_1<-data1[keep,]

  #Put data in wide form so we can see gold standard results
  data_wide_trueresult <- reshape2::dcast(data1, ID ~ time, value.var="true_result")
  truecens<-mean(data_wide_trueresult$`4`==0)
  GS_vis4_dat<-data_wide_trueresult[,c("ID","4")]
  colnames(GS_vis4_dat)<-c("ID","GS_vis4")

  #Create a variable indicating which interval the event first occurs in
  data_wide_trueresult$interval_occur<-ifelse(data_wide_trueresult$`1`==1,1,
                                              ifelse(data_wide_trueresult$`1`==0 & data_wide_trueresult$`2`==1,2,
                                                     ifelse(data_wide_trueresult$`1`==0 & data_wide_trueresult$`2`==0 & data_wide_trueresult$`3`==1,3,
                                                            ifelse(data_wide_trueresult$`1`==0 & data_wide_trueresult$`2`==0 & data_wide_trueresult$`3`==0 & data_wide_trueresult$`4`==1,4,5))))
  EventsDist<-table(data_wide_trueresult$interval_occur)

  #Consider that a certain proportion of gold standard outcomes are randomly missing (MCAR)
  mcar<-runif(N,0,1)
  GS_vis4_dat$GS_vis4_mis<-ifelse(mcar<prop_m,NA,GS_vis4_dat$GS_vis4)

  #Create final data sets for the proposed method analysis and gold standard analysis
  datafinal<-merge(datafinal_1,GS_vis4_dat,by="ID")
  datafinal_GS<-merge(data1,GS_vis4_dat,by="ID")

  #Create formula: just a single covariate x in the outcome model
  formula=result~x_1

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
  #For a random sample, weights are 1 for everyone
  GSdelta <- datafinal[uid,"GS_vis4_mis"]
  GSVis<-rep(4,nsub)
  weights<-rep(1,nsub)
  
  #Create interval censored data
  IC_data<-datafinal_GS[uid,c("true_result","x_1","GS_vis4_mis")]
  
  #Create data assuming gold standard observed every 4 years
  IC_data4year<-datafinal_GS[,c("ID","time","true_result","x_1","GS_vis4_mis")]

  #Keep complete data only in case we are missing gold standard
  IC_data_complete<-IC_data[complete.cases(IC_data),]
  IC_data4year_complete<-IC_data4year[complete.cases(IC_data),]

  myid.uni <- unique(datafinal_GS$ID)
  a<-length(myid.uni)
  last <- c()

  #Take last row of data in long form to be GS indicator data
  for (i in 1:a) {
    temp<-subset(datafinal_GS, ID==myid.uni[i])
    if (dim(temp)[1] > 1) {
      last.temp<-temp[dim(temp)[1],]
    }
    else {
      last.temp<-temp
    }
    last<-rbind(last, last.temp)
  }

  last2<-last[complete.cases(GS_vis4_dat$GS_vis4_mis),]

  #Only keep positive gold standard after first positive for case where observed all four years
  keepafterpos<-unlist(tapply(IC_data4year_complete$true_result,IC_data4year_complete$ID,after_first_pos))
  datafinal_4vis<-IC_data4year_complete[keepafterpos,]


  #Fit model: grouped time survival approach using  gold standard outcome observed at all 4 years
  fit_truth4vis<-glm(true_result~time+x_1,family=binomial(link="cloglog"),data=datafinal_4vis)
  fitsum_truth4Vis<-summary(fit_truth4vis)
  beta1est_truth4Vis<-fitsum_truth4Vis$coefficients["x_1",1]
  beta1se_truth4Vis<-fitsum_truth4Vis$coefficients["x_1",2]

  #Fit model: interval standard, gold standard where gold standard only observed at year 4
  fit_truth4Year<-glm(true_result~x_1,family=binomial(link="cloglog"),data=last2)
  fitsum_truth4Year<-summary(fit_truth4Year)
  beta1est_truth4Year<-fitsum_truth4Year$coefficients["x_1",1]
  beta1se_truth4Year<-fitsum_truth4Year$coefficients["x_1",2]

  #Create initial parameters for optim for proposed method
  lami <- log(-log(seq(1, 0.1, length.out = J + 1)[-1]))
  lami <- c(lami[1], diff(lami))
  tosurv <- function(x) exp(-exp(cumsum(x)))
  lowlam <- c(-Inf, rep(0, J - 1))
  lowerLBFGS <- c(rep(-Inf, nbeta),lowlam)
  parmi <- c(beta1_0,lami)

  #fit proposed model
  #functions "log_like_proposed" and "gradient_proposed" come from file "PROPOSED_AUGMENTEDLIKELIHOOD_FUNCTIONS.R" that 
  proposed_fit<-optim(par=parmi, log_like_proposed,gradient_proposed,lower = lowerLBFGS,upper=c(Inf,Inf,Inf,Inf,Inf),method = "L-BFGS-B",nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights,
            purpose="SUM")
  result_var<-solve(pracma::hessian(log_like_proposed,x=proposed_fit$par))
  beta1est_aux<-  proposed_fit$par[1]
  beta1se_aux<-sqrt(result_var[1,1])

  ansmat<-rbind(ansmat,c(beta1est_aux,beta1se_aux,beta1est_truth4Year,beta1se_truth4Year,beta1est_truth4Vis,beta1se_truth4Vis,truecens,EventsDist))

  print(iter)
}

ansmat<-data.frame(ansmat)
names(ansmat)<-c("beta1est_aux","beta1se_aux","beta1est_truth4Year","beta1se_truth4Year","beta1est_truthAll4Vis","beta1se_truthAll4Vis","CR","int1","int2",
                 "int3","int4","int5")

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
Data_output_singlecov_median(ansmat$beta1est_aux,ansmat$beta1se_aux,log(1.5))

#Gold Standard, Interval Censored
Data_output_singlecov_median(ansmat$beta1est_truth4Year,ansmat$beta1se_truth4Year,log(1.5))

#Gold Standard all four years, Grouped Survival
Data_output_singlecov_median(ansmat$beta1est_truthAll4Vis,ansmat$beta1se_truthAll4Vis,log(1.5))





