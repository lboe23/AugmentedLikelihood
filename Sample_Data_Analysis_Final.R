#############################################################################################################################################################################################################
#Code for running a sample data analysis from: "An Augmented Likelihood Approach for the Discrete Proportional Hazards Model Using Auxiliary and Validated Outcome Data - with Application to HCHS/SOL Study"
#############################################################################################################################################################################################################

#We need to load in our own functions for the log-likelihood/gradient and some functions from Gu et al. 2015 that help our calculations 
#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('U:/Paper 2/Code/RcppFunc.cpp')
source('U:/Paper 2/Code/Paper 2 Likelihood FINAL.R')


library(survey)
library(dplyr)
library(MASS)
library(stringr)
library(data.table)

after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}

#Set sensitivity and specificity - correct values needed for analysis
sensitivity<-0.61
specificity<-0.98

#Set proportion of missingness in gold standard
prop_m<-0.20

#Load in sample simulated data
  load(file=paste0('U:/Paper 2/SampleData.RData'))
 
  summary(samp$xstarstar)
  N<-dim(samp)[1]

#Fitting calibration model
    samp.solnas <- samp[(solnas==T),]
    lm.lsodi <- glm(xstarstar ~  xstar+z1+z2,data=samp.solnas) #univariable RC
    x.lsodi <- model.matrix(lm.lsodi) #X
    xtx.lsodi <- t(x.lsodi)%*%x.lsodi #X`X`
    ixtx.lsodi <- solve(xtx.lsodi) #(X`X)^-1
    samp[,xhat := predict(lm.lsodi,newdata=samp,'response')]

#put samp in long form
  samp_long1<-reshape(data = samp,
                      idvar = "ID",
                      varying = list(true_result=c("true_result_1","true_result_2","true_result_3","true_result_4","true_result_5","true_result_6",
                                                   "true_result_7","true_result_8"),result=c("result_1","result_2","result_3","result_4","result_5",
                                                                                             "result_6","result_7","result_8")),
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
  set.seed(2548)
  mcar<-runif(N,0,1)

#create an indicator variable called GS_vis4_mis indicting whether each subject missing GS or not
  GS_data$GS_vis4_mis<-ifelse(mcar<prop_m,NA,GS_data$GS_vis4)

#Create two FINAL datasets: datafinal, in long form with one row per visit and datafinal_GS, with 1 row per person for interval censored analysis only
  datafinal<-merge(datafinal_1,GS_data,by="ID")
  datafinal_GS<-merge(samp_long,GS_data[,c("ID","GS_vis4_mis")],by="ID")

#Write down formula for survival model: outcome is error-prone result Y* (result), covariate is x (X1_cov)
  formula=result~xhat+z1+z2

#Now, let's make sure our data is ordered properly before we begin calculating sum of the likelihood components.
  id <- eval(substitute(datafinal$ID), datafinal, parent.frame())
  time <- eval(substitute(datafinal$time), datafinal, parent.frame())
  result <- eval(substitute(result), datafinal, parent.frame())
  ord <- order(id, time)
  if (is.unsorted(ord)) {
    id <- id[ord]
    time <- time[ord]
    result <- result[ord]
    datafinal <- datafinal[ord, ]}
  utime <- sort(unique(time))
  timen0 <- (time != 0)

#Calculate D matrix and C matrix for log likelihood
  Dm <- dmat(id[timen0], time[timen0], result[timen0], sensitivity,
             specificity, 1)
  Cm <- cmat(id[timen0], time[timen0], result[timen0], sensitivity,
             specificity, 1)

#Save off J, number of visits
  J <- ncol(Dm) - 1
  nsub <- nrow(Dm)

#Create Xmatrix, which is our model matrix with one column per covariate and currently is in long form (one row per visit)
  Xmat <- model.matrix(formula, data = datafinal)[, -1, drop = F]
  beta.nm <- colnames(Xmat)
  nbeta <- ncol(Xmat)
  uid <- getrids(id, nsub)

#Create unique Xmat with only one row per person
  Xmat <- Xmat[uid, , drop = F]

#Create unique vector with only one row per person indicating whether each person had GS or not
#This will be used to calculate likelihood for proposed estimator...if else based on if GSdelta=NA, 0 or 1.
  GSdelta <- datafinal[uid,"GS_vis4_mis"]
  GSVis<-rep(4,nsub)
  noweights<-rep(1,nsub)
  
#Create interval censored data - unique ID's only so data has N rows
  IC_data<-datafinal_GS[!duplicated(datafinal_GS$ID),c("result","true_result","GS_vis4","GS_vis4_mis","BGid", "strat", "bghhsub_s2","xstar","xhat","z1","z2")]
#Now create final dataset for standard IC analysis, omitting anyone with missing GS
  IC_GS_datafinal<-IC_data[complete.cases(IC_data$GS_vis4_mis),]

#Create starting values for our survival parameters: first, to avoid maximization problems due to the ordered constraint of the survival parameters,
#we re-parameterize using a log-log transformation of survival function for S_2, and a change in log-log of the survival function for all others.
#We consider initial values of 0.1 for our survival parameters, then transform these based on this re-parameterization.
  initsurv <- 0.1
  lami <- log(-log(seq(1, initsurv, length.out = J + 1)[-1]))
  lami <- c(lami[1], diff(lami))
  tosurv <- function(x) exp(-exp(cumsum(x)))

#Define a lower bound of infinity for the first re-parameterized survival function and 0 for the subsequent J-1 terms.
  lowlam <- c(-Inf, rep(0, J - 1))
  lowerLBFGS <- c(rep(-Inf, nbeta),lowlam)
  upperLBFGS <- c(rep(Inf, nbeta+J))

#Create vector parmi consisting of a starting value for beta, and starting values for re-parameterized survival paramters
parmi <- c(rep(0.5,nbeta),lami)

#Create survey designs using survey package for the two models
  samp_design_reg = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=IC_GS_datafinal)
  samp_design_reg_complete = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=IC_data)

#Fit proposed method and standard, IC approach for weighted estimators
    proposed_fit_data_weight<-optim(par=parmi, fn=log_like_proposed,gr=gradient_proposed,lower = lowerLBFGS,upper=upperLBFGS,method = "L-BFGS-B",nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights=weights,
                    purpose="SUM",hessian=TRUE)
    inverted_hessian<-solve(proposed_fit_data_weight$hessian)

    proposed_fit_data_weight_GS<-svyglm(GS_vis4~xhat+z1+z2,family=quasibinomial(link="cloglog"),data=IC_GS_datafinal,design=samp_design_reg)
    proposed_fit_data_weight_GS$coefficients
#Calculate sandwich variance to incorporate stratification and clustering
    U_prop1<-gradient_proposed(proposed_fit_data_weight$par,nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights=noweights,
                           purpose="INDIVIDUAL")
    infl1<- U_prop1%*%inverted_hessian
    mySandVar<- vcov(svytotal(infl1,samp_design_reg_complete,influence=TRUE))
    
#Save estimated parameters
    beta1est_aux_w<-  proposed_fit_data_weight$par[1]
    sandVar<-  sqrt(diag(mySandVar))[1]
  
    fitsum_truth4Year<-summary(proposed_fit_data_weight_GS)
    beta1est_truth4Year<-fitsum_truth4Year$coefficients["xhat",1]
    beta1se_truth4Year<-fitsum_truth4Year$coefficients["xhat",2]

#We now begin multiple imputation. Can either run this code (takes some time), or skip to line 241 where we load output from it 
#Number of multiple imputation runs
    M=25 
    
#Set seed for multiple imputation
    set.seed(2328)
    
#Create empty lists for saving results
    list_mi = list()

# Multiple Imputation starts now
for(j in 1:M){cat("Imputation: ",j,".\n")
  
  list_mi[[j]] = list()
  # Sampling new variance-covariance matrix from inverse chi-squared distribution
  samp.solnas.boot <- samp.solnas[sample(1:nrow(samp.solnas),nrow(samp.solnas),replace = T),]
  vcov.star <- (((sigma(lm.lsodi)^2)*(lm.lsodi$df.residual))/rchisq(n=1,df=lm.lsodi$df.residual))*ixtx.lsodi
  
  IC_data$xhat_mi <- with(IC_data,cbind(1,xstar,z1,z2))%*%mvrnorm(n=1,mu=as.numeric(lm.lsodi$coefficients),Sigma=vcov.star)
  IC_GS_datafinal<-IC_data[complete.cases(IC_data$GS_vis4_mis),]
  
  #update sampling design
  samp_design_reg <- update(samp_design_reg,xhat_mi = IC_GS_datafinal$xhat_mi)
  samp_design_reg_complete <- update(samp_design_reg_complete,xhat_mi = IC_data$xhat_mi)
  
  #outcome model formula again
  formula_mi=result~xhat_mi+z1+z2
  
  #Create Xmatrix, which is our model matrix with one column per covariate and currently is in long form (one row per visit)
  Xmat_mi <- model.matrix(formula_mi, data = IC_data)[, -1, drop = F]
  
  
  #Fit proposed method and standard, IC approach for weighted estimators
  proposed_fit_data_weight_mi<-optim(par=parmi, fn=log_like_proposed,gr=gradient_proposed,lower = lowerLBFGS,upper=upperLBFGS,method = "L-BFGS-B",nsub,J,nbeta,Dm,Cm,Xmat=Xmat_mi,GSdelta,GSVis,weights=weights,
                                     purpose="SUM",hessian=TRUE)
  inverted_hessian_mi<-solve(proposed_fit_data_weight_mi$hessian)
  
  proposed_fit_data_weight_GS_mi<-svyglm(GS_vis4~xhat_mi+z1+z2,family=quasibinomial(link="cloglog"),data=IC_GS_datafinal,design=samp_design_reg)
  
  #Calculate sandwich variance to incorporate stratification and clustering
  U_prop1_mi<-gradient_proposed(proposed_fit_data_weight_mi$par,nsub,J,nbeta,Dm,Cm,Xmat=Xmat_mi,GSdelta,GSVis,weights=noweights,
                                purpose="INDIVIDUAL")
  infl1_mi<- U_prop1_mi%*%inverted_hessian_mi
  mySandVar_mi<- vcov(svytotal(infl1_mi,samp_design_reg_complete,influence=TRUE))
  
  #Save estimated parameters
  beta1est_aux_w_mi<-  proposed_fit_data_weight_mi$par[1]
  sandVar_mi<-  sqrt(diag(mySandVar_mi))[1]
  
  fitsum_truth4Year_mi<-summary(proposed_fit_data_weight_GS_mi)
  beta1est_truth4Year_mi<-fitsum_truth4Year_mi$coefficients["xhat_mi",1]
  beta1se_truth4Year_mi<-fitsum_truth4Year_mi$coefficients["xhat_mi",2]
  
  # Saving parameters
  list_mi[[j]] <-data.frame(Imp=j,beta1est_aux_w_mi,sandVar_mi,beta1est_truth4Year_mi,beta1se_truth4Year_mi)
  
  
  list_mi[[j]] <- do.call(rbind,list_mi[[j]])
  
  print(j)
}

save(list_mi,file=paste0('U:/Paper 2/SampleAnalysis_MIVarianceResults3.RData'))

load(file=paste0('U:/Paper 2/SampleAnalysis_MIVarianceResults3.RData'))

output_mi<-as.data.frame(matrix(unlist(list_mi), ncol=5, byrow=T))
colnames(output_mi)<-c("Imp","beta_proposed_mi","se_proposed_mi","beta_standard_mi","se_standard_mi")

MI_Var<-function(V,betas){
  var<-(median(V)+(mad(betas)^2))
  return(sqrt(var))
}


se_proposed_final<-MI_Var((output_mi$se_proposed_mi)^2,output_mi$beta_proposed_mi)
se_standard_final<-MI_Var((output_mi$se_standard_mi)^2,output_mi$beta_standard_mi)



myfinaltable<-cbind(paste(round(exp(beta1est_aux_w*log(1.2)),2),"(",round(exp(beta1est_aux_w-1.96*se_proposed_final)^log(1.2),2),",",
                          round(exp(beta1est_aux_w+1.96*se_proposed_final)^log(1.2),2),")"),
                    paste(round(exp(beta1est_truth4Year*log(1.2)),2),"(",round(exp(beta1est_truth4Year-1.96*se_standard_final)^log(1.2),2),",",
                          round(exp(beta1est_truth4Year+1.96*se_standard_final)^log(1.2),2),")"),round((se_standard_final^2)/(se_proposed_final^2),2))


library(xtable)
xtable(myfinaltable)
