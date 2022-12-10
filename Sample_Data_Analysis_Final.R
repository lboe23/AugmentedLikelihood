#############################################################################################################################################################################################################
#Code for running a sample data analysis from: "An Augmented Likelihood Approach for the Cox Proportional Hazards Model with Interval-Censored Auxiliary and Validated Outcome Data – with Application to the HCHS/SOL Study"
#############################################################################################################################################################################################################

#We need to load in our own functions for the log-likelihood/gradient and some functions from Gu et al. 2015 that help our calculations,
# as well as our R functions for fitting proposed model.
#Note that the file and directory below will need to be changed in order for another user to get the code to run.
Rcpp::sourceCpp('RcppFunc.cpp')
source('PROPOSED_AUGMENTEDLIKELIHOOD_FUNCTIONS.R')


library(survey)
library(dplyr)
library(MASS)
library(stringr)
library(data.table)
library(xtable)
library(numDeriv)

#Set proportion of missingness in gold standard
prop_m<-0.20

#Load in sample simulated data
load(file=paste0('G:/Old Papers/SMMR/SampleData.RData'))

#Number of people in data set 
N<-dim(samp)[1]

#Now estimate sensitivity and specificity based on visit 4 measures, comparing auxiliary to gold standard
sens_spec_tbl<-table(True=samp$true_result_4,Self=samp$result_4)
Yi_sens<-sens_spec_tbl[4]
Yi_spec<-sens_spec_tbl[1]
n_sens<-sens_spec_tbl[2]+sens_spec_tbl[4]
n_spec<-sens_spec_tbl[3]+sens_spec_tbl[1]
sensitivity<-Yi_sens/n_sens
specificity<-Yi_spec/n_spec

#Create variables indicating whether or not someone contributes to sensitivity and specificity estimating equations'
samp$Se_EstEq<-ifelse(samp$true_result_4==1,1,0)
samp$Sp_EstEq<-ifelse(samp$true_result_4==0,1,0)



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

#Create interval censored data - unique ID's only so data has N rows
IC_data<-datafinal_GS[!duplicated(datafinal_GS$ID),c("ID","result","true_result","GS_vis4","GS_vis4_mis","BGid",
                                                     "strat", "bghhsub_s2","xstar","z1","z2","xstarstar",
                                                     "Se_EstEq","Sp_EstEq","solnas")]

#Strata- cross-classify original strata and subset indicator
IC_data$newstrata<-interaction(IC_data$strat,IC_data$solnas)

#Create survey designs using survey package for the two models
samp_design_reg_complete = svydesign(id=~BGid, strata=~newstrata, weights=~bghhsub_s2, data=IC_data,nest=TRUE)

#Fitting calibration model
#Create survey design - subset to solnas design
samp.solnas.design <- subset(samp_design_reg_complete,solnas==T)
lm.calib <- svyglm(xstarstar ~  xstar+z1+z2,design=samp.solnas.design) 

#Predict xhat in data and in survey design
datafinal$xhat<-predict(lm.calib, newdata=datafinal)
samp_design_reg_complete <- update(samp_design_reg_complete,xhat =predict(lm.calib,
                                              newdata=samp_design_reg_complete$variables) )

#Now that we have xhat, further subset designs for IC, GS analysis as well as sandwich for GS analysis
samp_design_reg = subset(samp_design_reg_complete,GS_vis4_mis==1 | GS_vis4_mis==0)
samp_design_reg_sandwich = subset(samp_design_reg_complete,complete.cases(GS_vis4_mis) | solnas==T)

#Write down formula for survival model: outcome is error-prone result Y* (result), covariate is x (X1_cov)
formula=result~xhat+z1+z2
formula_naive=result~xstar+z1+z2


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
#Xmat for proposed apporach with estimated exposure AND for naive
Xmat <- model.matrix(formula, data = datafinal)[, -1, drop = F]
Xmat_naive <- model.matrix(formula_naive, data = datafinal)[, drop = F]

beta.nm <- colnames(Xmat)
nbeta <- ncol(Xmat)
uid <- getrids(id, nsub)

#Create unique Xmat with only one row per person
Xmat <- Xmat[uid, , drop = F]
Xmat_naive <- Xmat_naive[uid, , drop = F]

#Create unique vector with only one row per person indicating whether each person had GS or not
#This will be used to calculate likelihood for proposed estimator...if else based on if GSdelta=NA, 0 or 1.
GSdelta <- datafinal[uid,"GS_vis4_mis"]
GSVis<-rep(4,nsub)
noweights<-rep(1,nsub)

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


#Fit proposed method 
proposed_fit_data_weight<-optim(par=parmi, fn=log_like_proposed,gr=gradient_proposed,lower = lowerLBFGS,upper=upperLBFGS,method = "L-BFGS-B",nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights=weights,
                                purpose="SUM",hessian=TRUE)
inverted_hessian<-solve(proposed_fit_data_weight$hessian)

#Fit standard, IC approach 
proposed_fit_data_weight_GS<-svyglm(GS_vis4~xhat+z1+z2,family=quasibinomial(link="cloglog"),design=samp_design_reg)
#Fit standard, IC approach with error-prone exposure
proposed_fit_data_weight_GS_naive<-svyglm(GS_vis4~xstar+z1+z2,family=quasibinomial(link="cloglog"),design=samp_design_reg)

#Calculate sandwich variance to incorporate stratification and clustering and extra uncertainty added
#   by estimated exposure, Se, and Sp for Proposed Approach
########################################################################################################
  params=proposed_fit_data_weight$par
  k_dim<-length(proposed_fit_data_weight$par)
  j_dim<-length(coef(lm.calib))
  
  #Now create jacobians for these functions
  JacobianAll.alphas<-numDeriv::jacobian(func=gradient_proposed_alphas, x=coef(lm.calib))
  JacobianAll.sens<-numDeriv::jacobian(func=gradient_proposed_sens, x=sensitivity)
  JacobianAll.spec<-numDeriv::jacobian(func=gradient_proposed_spec, x=specificity)
  
  #Split these matrices into a list
  JacobianList.alphas<-lapply(split(JacobianAll.alphas,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  JacobianList.sens<-lapply(split(JacobianAll.sens,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  JacobianList.spec<-lapply(split(JacobianAll.spec,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  
  #Create components of A matrices
  A.U1.U1<- numDeriv::jacobian(func=esteqsens, x=sensitivity) 
  A.U1.U2<- matrix(0,nrow=1,ncol=1)
  A.U1.U3<- matrix(0,nrow=1,ncol=j_dim)
  A.U1.U4<- matrix(0,nrow=1,ncol=k_dim)
  A.U2.U1<- matrix(0,nrow=1,ncol=1)
  A.U2.U2<- numDeriv::jacobian(func=esteqspec, x=specificity) 
  A.U2.U3<- matrix(0,nrow=1,ncol=j_dim)
  A.U2.U4<- matrix(0,nrow=1,ncol=k_dim)
  A.U3.U1<- matrix(0,nrow=j_dim,ncol=1)
  A.U3.U2<- matrix(0,nrow=j_dim,ncol=1)
  A.U3.U3<- -solve(lm.calib$naive.cov/mean(lm.calib$prior.weights))/N
  A.U3.U4<- matrix(0,nrow=j_dim,ncol=k_dim)
  A.U4.U1<-add(JacobianList.sens)*mean(weights)/N
  A.U4.U2<-add(JacobianList.spec)*mean(weights)/N
  A.U4.U3<-add(JacobianList.alphas)*mean(weights)/N
  A.U4.U4<- proposed_fit_data_weight$hessian*mean(weights)/N
  
  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll<-rbind(cbind(A.U1.U1,A.U1.U2,A.U1.U3,A.U1.U4),cbind(A.U2.U1,A.U2.U2,A.U2.U3,A.U2.U4),
              cbind(A.U3.U1,A.U3.U2,A.U3.U3,A.U3.U4),cbind(A.U4.U1,A.U4.U2,A.U4.U3,A.U4.U4))
  
  #Invert this matrix; needed for sandwich
  A.inv<-solve(AAll)
  
  #Now estfuns for all 4 EE - create empty ones that we will fill in
  estfun.Sens<-matrix(0,nrow=N,ncol=1)
  estfun.Spec<-matrix(0,nrow=N,ncol=1)
  estfun.calib<-matrix(0,nrow=N,ncol=length(coef(lm.calib)))

  #Now fill in pieces of estimating equations depending on who contributes to each
  estfun.Sens[samp_design_reg_complete$variables$Se_EstEq==1,]<-esteqsens(sensitivity)
  estfun.Spec[samp_design_reg_complete$variables$Sp_EstEq==1,]<-esteqspec(specificity)
  estfun.calib[samp_design_reg_complete$variables$solnas==1,]<-estfuns.lmNEW(lm.calib)
  estfun.outcome<-gradient_proposed(proposed_fit_data_weight$par,nsub=nsub,J=J,nbeta=nbeta,Dm=Dm,
                                    Cm=Cm,Xmat=Xmat,GSdelta=GSdelta,GSVis=GSVis,weights=noweights,
                                    purpose="INDIVIDUAL")
 
  #Combine estfuns for Se, Sp, calibration and outcome models
  estfun.all<-cbind(estfun.Sens,estfun.Spec, estfun.calib,estfun.outcome)
  
  #Compute influence functions
  infl<- as.matrix(estfun.all)%*%t(A.inv)/N
  
  #Compute the sandwich using convenient survey package functions, applying vcov(svytotal())
  sand1<-vcov(svytotal(infl,samp_design_reg_complete))

#Calculate sandwich variance to incorporate stratification and clustering and extra uncertainty added
#   by estimated exposure for Gold-Standard, Interval-Censored Approach
########################################################################################################
  
  #Save estimated sample size for Gold Standard Approach (smaller due to missingness)
  N_GoldStd<-dim(samp_design_reg)[1]
  N_GoldStd_All<-dim(samp_design_reg_sandwich)[1]
  
  #A matrix for IC approach is just four quadrants-- start with upper right, all zero
  A.upperright.IC<- matrix(0,nrow=length(lm.calib$coefficients),ncol=length(proposed_fit_data_weight_GS$coefficients))
  
  #Here is the upper left computed using the Hessian from calibration model
  A.upperleft.IC<- -solve(lm.calib$naive.cov/mean(lm.calib$prior.weights))/N_GoldStd_All
  
  #Bottom-right quadrant is just Hessian from outcome model
  A.bottomright.IC<- -solve(proposed_fit_data_weight_GS$naive.cov/mean(proposed_fit_data_weight_GS$prior.weights))/N_GoldStd
  
  #Now, use Jacobian function from numderiv to obtain derivatives of outcome model estimating equation wrt alpha
  #Then split it up, because it computes large matrix with each person's contributions for each parameter
  JacobianAll.outcome.IC<-numDeriv::jacobian(func=outcome.alphas.estfuns.glm, x=coef(lm.calib))
  JacobianList.outcome.IC<-lapply(split(JacobianAll.outcome.IC,rep(c(1:N_GoldStd),each=1)),
                                  matrix,nrow=length(proposed_fit_data_weight_GS$coefficients))
  
  #Add these matrices for subjects 1...N_GoldStd then divide by N_GoldStd
  A.bottomleft.IC<-add(JacobianList.outcome.IC)*mean(proposed_fit_data_weight_GS$prior.weights)/N_GoldStd
  
  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll.IC<-rbind(cbind(A.upperleft.IC,A.upperright.IC),cbind(A.bottomleft.IC,A.bottomright.IC))
  
  #Invert this matrix; needed for sandwich
  A.inv.IC<-solve(AAll.IC)
  
  #Calculate estimating equation contributions for calibration and outcome models
  estfun.calib.IC<-matrix(0,nrow= N_GoldStd_All,ncol=length(coef(lm.calib)))
  estfun.outcome.IC<-matrix(0,nrow= N_GoldStd_All,ncol=length(coef(proposed_fit_data_weight_GS)))
  
  #Now create estimating equation contributions for calibration model (0 if not in model)
  estfun.calib.IC[samp_design_reg_sandwich$variables$solnas==1 ,] <- estfuns.lmNEW(lm.calib)
  estfun.outcome.IC[samp_design_reg_sandwich$variables$ID %in% samp_design_reg$variables$ID,] <- 
    estfuns.glmNEW(proposed_fit_data_weight_GS)/(proposed_fit_data_weight_GS$prior.weights)
  
  #Combine estfuns for both models
  estfun.all.IC<-cbind(estfun.calib.IC,estfun.outcome.IC)
  
  #Compute influence functions
  infl.IC<- as.matrix(estfun.all.IC)%*%t(A.inv.IC)/N_GoldStd_All
  
  #Compute the sandwich using vcov(svytotal())
  sand2<-vcov(svytotal(infl.IC,samp_design_reg_sandwich))


#Now save final results
beta_proposed<-proposed_fit_data_weight$par[1]
beta_standard_IC<-proposed_fit_data_weight_GS$coefficients[2]
se_proposed_final<-sqrt(diag(sand1))[7]
se_standard_final<-sqrt(diag(sand2))[6]


#Combine results into a final table
myfinaltable<-cbind(paste(round(exp(beta_proposed*log(1.2)),2),"(",round(exp(beta_proposed-1.96*se_proposed_final)^log(1.2),2),",",
                            round(exp(beta_proposed+1.96*se_proposed_final)^log(1.2),2),")"),
                      paste(round(exp(beta_standard_IC*log(1.2)),2),"(",round(exp(beta_standard_IC-1.96*se_standard_final)^log(1.2),2),",",
                            round(exp(beta_standard_IC+1.96*se_standard_final)^log(1.2),2),")"),round((se_standard_final^2)/(se_proposed_final^2),2))
  
  
xtable(myfinaltable)
