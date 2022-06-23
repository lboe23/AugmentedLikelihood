#(1) We begin by generating a list of individuals using a three-stage design with stratification, block groups, and households.
lapply(c('data.table','mvtnorm','magrittr','survey','dplyr','tidyr','ggplot2','ggpubr','kableExtra','stringr'),library,character.only = T)

#Note that we need to load in multiNorm.RData, logiReg.RData, calibrCoeff.RData, and outPar.Rdata from
#https://github.com/plbaldoni/HCHSsim/tree/master/data

#Now begin by assigning parameters and choose covariate distribution

#Settings to vary: average sample size, censoring rate for latent true event time
#meanN<-1000
meanN<-10000
#myCR<-0.9
#myCR<-0.7
myCR<-0.5


#Pick covariate distribution
#covardist<-"Normal"
covardist<-"Gamma"

### Assign population name
popname = paste0('TargetPopulationData',covardist,'_CR',myCR)


####### generate HCHS population #######
# Set the seed
set.seed(123)

# Number of BG (block groups) in the population
Nbg=376

# Cumulative sum of BG in the population across all 4 strata
bgstrat = cumsum(c(58,21,130,167))

# Number of HH per BG, at least 1 HH per BG, mean number of HHs/BG=451 (this is the number of ALL HHs in each BG including non-eligible HHs)
bg.size.temp=1+round(rexp(Nbg,1/450))

# Number of eligible HHs with Hispanic surname in each BG
# The proportion of eligible HHs with Hispanic surname per stratum is 0.2, 0.2125, 0.2975, and 0.2975
num.hisp.strat=round(c(bg.size.temp[1:bgstrat[1]]*.2,bg.size.temp[(bgstrat[1]+1):bgstrat[2]]*.2125,bg.size.temp[(bgstrat[2]+1):bgstrat[3]]*.2975,bg.size.temp[(bgstrat[3]+1):bgstrat[4]]*.2975))

# Number of eligible HHs with Other surname in each BG
# The proportion of eligible HHs with Hispanic surname per stratum is 0.15, 0.225, 0.26, and 0.2925
num.other.strat=round(c(bg.size.temp[1:bgstrat[1]]*.15,bg.size.temp[(bgstrat[1]+1):bgstrat[2]]*.225,bg.size.temp[(bgstrat[2]+1):bgstrat[3]]*.26,bg.size.temp[(bgstrat[3]+1):bgstrat[4]]*.2925))

# Number of eligible HHs in each BG
# Note that bg.size < bg.size.temp, as the latter contains non-eligible HHs
bg.size=num.hisp.strat+num.other.strat

# Number of HHs in target population
Nhh=sum(bg.size)

# HH size, at least 1 subject per HH, mean number of subjects/HH=2
hh.size=1+rpois(Nhh,1)

# Number of subjects in target population
N=sum(hh.size)

# all ID's unique (e.g., subid=k only for one subject within one HH)
BGid=rep(rep(rep(1:Nbg, times=bg.size), times=hh.size), times=2)
hhid=rep(rep(1:Nhh, times=hh.size), times=2)
subid=rep(1:N, times=2)
v.num=rep(c(1,2),each=N)
strat=rep(NA,times=Nbg); strat[BGid<=bgstrat[1]]=1; strat[BGid>bgstrat[1] & BGid<=bgstrat[2]]=2; strat[BGid>bgstrat[2] & BGid<=bgstrat[3]]=3; strat[BGid>bgstrat[3] & BGid<=bgstrat[4]]=4


# Now creating hips.strat variable
# It creates a matrix with Nbg columns (BGs) and it fills up with TRUE (if hispanic and in target population),
# FALSE (if other surname and in target population), and NA (if otherwise)
A=matrix(NA,nrow=max(bg.size),ncol=length(bg.size))
for (i in 1:length(bg.size)){
  A[,i]=c(rep(TRUE,times=num.hisp.strat[i]),rep(FALSE,times=num.other.strat[i]),rep(NA,times=max(bg.size)-bg.size[i]))    # create matrix A with each column corresponding to a BG (containing 1's for each Hispanic HH followed by 0's for each other HH)
}
# Indicator for Hispanic surname (na.omit(c(A)) gives the household-level hispanic surname)
hisp.strat=rep(rep(c(na.omit(c(A))),times=hh.size),times=2)

pop = data.table(strat,BGid,hhid,hisp.strat,subid,v.num)
pop = pop[order(strat,BGid,hhid,hisp.strat,subid,v.num),]

# Sanity check: Check if there is any household with multiple distinct hispanic surname (should be 0 or 1)
pop[,.(Check = mean(hisp.strat)%in%c(0,1)),by='hhid'][,all(Check)]
# Check the proportion of hispanic surname per Strat (should be somewhat close to 0.57, 0.48, 0.53, 0.50)
pop[,mean(hisp.strat),by='strat']
# Check how many BGid per stratum we have (it should be 58, 21, 130, and 167)
pop[,unique(BGid),by='strat'][,.N,by='strat']

# Saving output
if(!dir.exists('./RData')){system('mkdir RData')}
save(pop,file='U:/TargetPopulation.RData',compress = 'xz')


#(2) Next, populate the target population data sets with several variables to be used in this simulation study.

after_first_pos <- function(x){
  npos<-cumsum(x==1)
  (npos==0) | (npos==1 & x==1)
}

mybeta<- log(1.5)

#test times is vector of pre-scheduled test times
sensitivity<-0.80
specificity<-0.90
testtimes<-c(1,2,3,4)
ntest <- length(testtimes)

#Exponential Baseline Hazard
blambda<-ifelse(myCR==0.5,0.17,ifelse(myCR==0.7,0.08,ifelse(myCR==0.9,0.023,NA)))

### Expit function
expit = function(xb){return(exp(xb)/(1+exp(xb)))}

### Setting seed
set.seed(04142019)
### Loading data for Multivariate Gaussian ###
load('U://multiNorm.RData')

### Loading parameter estimates to simulate binary variables
load('U://logiReg.RData')

### Loading Target Population Data ###
load(paste0('U://TargetPopulation.RData'))
pop = pop[v.num==1,]

### sim_sexbkg: simulate multinomial distribution of sex*background ###
sim_sexbkg = function(dat,p,sex=c('F','M'),bkg=c('D','PR','O')){
  tb =  expand.grid(sex,bkg);colnames(tb) = c('sex','bkg')
  tb$idx = 1:nrow(tb)
  tb$prob = p
  
  dat$idx = apply(rmultinom(nrow(dat),size=1,prob=tb$prob),2,which.max)
  dat = merge(dat,subset(tb,select=-prob),by='idx',all.x=T)
  dat = dat[order(dat$strat,dat$BGid,dat$hhid,dat$subid),]
  return(dat)
}

### sim_covar: simulate covariates ###
sim_covar = function(dat,ranges,means,covs,dat.par,count.seed=1,
                     covars=c('AGE','BMI','LN_NA_AVG','LN_K_AVG','LN_KCAL_AVG','LN_PROTEIN_AVG')){
  dat <- as.data.frame(dat)
  ncovars = length(covars)
  df = data.frame()
  groups = names(table(dat$idx))
  namevar = c('age','bmi','ln_na_avg','ln_k_avg','ln_kcal_avg','ln_protein_avg')
  
  m = list()
  v = list()
  
  # F=Female, M=Male
  # D= Dominican, PR=Puerto Rican, O=other
  m[[1]] = as.numeric(means[1,]) #F,D
  m[[2]] = as.numeric(means[2,]) #M,D
  m[[3]] = as.numeric(means[3,]) #F,PR
  m[[4]] = as.numeric(means[4,]) #M,PR
  m[[5]] = as.numeric(means[5,]) #F,O
  m[[6]] = as.numeric(means[6,]) #M,O
  
  v[[1]] = matrix(as.numeric(covs[1,]),ncovars,ncovars,byrow = T) #F,D
  v[[2]] = matrix(as.numeric(covs[2,]),ncovars,ncovars,byrow = T) #M,D
  v[[3]] = matrix(as.numeric(covs[3,]),ncovars,ncovars,byrow = T) #F,PR
  v[[4]] = matrix(as.numeric(covs[4,]),ncovars,ncovars,byrow = T) #M,PR
  v[[5]] = matrix(as.numeric(covs[5,]),ncovars,ncovars,byrow = T) #F,O
  v[[6]] = matrix(as.numeric(covs[6,]),ncovars,ncovars,byrow = T) #M,O
  
  for(i in groups){cat('Group: ',i,".\n")
    subdat = subset(dat,idx==i)
    
    m.topass = m[[as.numeric(i)]]
    v.topass = v[[as.numeric(i)]]
    
    ### Simulating Age, BMI, and log-transformed nutrients
    subdat = cbind(subdat,rmvnorm(nrow(subdat),mean=m.topass,sigma=v.topass))
    
    #Simulate a single covariate for the proportional hazards model
    #Consider gamma or normal covariate for complex survey design
    X1 = rgamma(nrow(subdat), shape=1/5, scale = 1)
    
    # X1 = rnorm(nrow(subdat), 1/5, 1)
    
    #Combine my covariate with rest of data 
    subdat = cbind(subdat,X1)
    
    colnames(subdat) = c(colnames(dat),namevar)
    
    idx.age = with(subdat,which(!between(age,ranges[1,1],ranges[1,2]) |
                                  !between(bmi,ranges[2,1],ranges[2,2]) |
                                  !between(ln_na_avg,ranges[3,1],ranges[3,2]) |
                                  !between(ln_k_avg,ranges[4,1],ranges[4,2]) |
                                  !between(ln_kcal_avg,ranges[5,1],ranges[5,2]) |
                                  !between(ln_protein_avg,ranges[6,1],ranges[6,2])))
    ### Checking for values out of the range in the population data
    for(j in idx.age){cat(count.seed,"\r")
      flag = T
      while(flag){
        subdat[j,namevar] <- rmvnorm(1,mean=m.topass,sigma=v.topass)
        count.seed = count.seed + 1
        flag = with(subdat,!between(age[j],ranges[1,1],ranges[1,2]) |
                      !between(bmi[j],ranges[2,1],ranges[2,2]) |
                      !between(ln_na_avg[j],ranges[3,1],ranges[3,2]) |
                      !between(ln_k_avg[j],ranges[4,1],ranges[4,2]) |
                      !between(ln_kcal_avg[j],ranges[5,1],ranges[5,2]) |
                      !between(ln_protein_avg[j],ranges[6,1],ranges[6,2]))
      }
    }
    ### Now, simulating binary variables ###
    subdat.par = subset(dat.par,idx==i & MODEL =='USBORN',select=c(idx,Variable,Estimate))
    subdat.par = subdat.par[match(c("Intercept","AGE","BMI","LN_NA_AVG","LN_K_AVG","LN_KCAL_AVG","LN_PROTEIN_AVG"),subdat.par$Variable),]
    p.par = expit(as.matrix(cbind('Intercept'=1,subdat[,namevar]))%*%subdat.par$Estimate)
    subdat$usborn = as.numeric(1*(runif(nrow(p.par))<=p.par))
    
    subdat.par = subset(dat.par,idx==i & MODEL =='CHOL',select=c(idx,Variable,Estimate))
    subdat.par = subdat.par[match(c("Intercept","AGE","BMI","LN_NA_AVG","LN_K_AVG","LN_KCAL_AVG","LN_PROTEIN_AVG","US_BORN"),subdat.par$Variable),]
    p.par = expit(as.matrix(cbind('Intercept'=1,subdat[,c(namevar,'usborn')]))%*%subdat.par$Estimate)
    subdat$high_chol = as.numeric(1*(runif(nrow(p.par))<=p.par))
    
    ### Saving data
    df = rbind(df,subdat)
  }
  df1 = df[order(df$strat,df$BGid,df$hhid,df$subid),];rm(df)
  df1$age.strat = (df1$age>=45)
  df1$v.num=1
  
  #df2 = df1
  #df2$v.num = 2
  
  #df = rbind(df1,df2)
  df = df1
  df = df[order(df$strat,df$BGid,df$hhid,df$subid,df$v.num),]
  
  cat('Done!',"\n")
  df = as.data.table(df)
  
  return(list('Data'=df,'HCHS.Mean'=m,'HCHS.Cov'=v))
}

ls.pop = sim_sexbkg(dat=pop,p=c(0.1953,0.1300,0.2065,0.1984,0.1421,0.1277))
rls.pop = sim_covar(dat=ls.pop,ranges = ranges,means = means,covs = covs,
                    dat.par = dat.par)
ALL_N<-dim(ls.pop)[1]

#Now let's try adding gamma covariates where distributions differ by strata and bock group-
#Start by assigning base alpha and beta to vary by 4 strata
#Then use a uniform to add 15% change to each 
nbg1<-length(table(STRAT=rls.pop$Data[rls.pop$Data$strat==1,]$strat,BG=rls.pop$Data[rls.pop$Data$strat==1,]$BGid))
nbg2<-length(table(STRAT=rls.pop$Data[rls.pop$Data$strat==2,]$strat,BG=rls.pop$Data[rls.pop$Data$strat==2,]$BGid))
nbg3<-length(table(STRAT=rls.pop$Data[rls.pop$Data$strat==3,]$strat,BG=rls.pop$Data[rls.pop$Data$strat==3,]$BGid))
nbg4<-length(table(STRAT=rls.pop$Data[rls.pop$Data$strat==4,]$strat,BG=rls.pop$Data[rls.pop$Data$strat==4,]$BGid))

tablestrat1<-table(BG=rls.pop$Data[rls.pop$Data$strat==1,]$BGid,STRAT=rls.pop$Data[rls.pop$Data$strat==1,]$strat)
tablestrat2<-table(BG=rls.pop$Data[rls.pop$Data$strat==2,]$BGid,STRAT=rls.pop$Data[rls.pop$Data$strat==2,]$strat)
tablestrat3<-table(BG=rls.pop$Data[rls.pop$Data$strat==3,]$BGid,STRAT=rls.pop$Data[rls.pop$Data$strat==3,]$strat)
tablestrat4<-table(BG=rls.pop$Data[rls.pop$Data$strat==4,]$BGid,STRAT=rls.pop$Data[rls.pop$Data$strat==4,]$strat)

#Assign different alphas and betas for each stratum 
alpha1<-.25
beta1<-1.25
alpha2<-.15
beta2<-.75
alpha3<-.3
beta3<-1.5
alpha4<-.1
beta4<-.5

#Assign epsilons - slight shift for each BG
epsilon_alpha1<-runif(nbg1,-.15*alpha1,.15*alpha1)
epsilon_alpha2<-runif(nbg2,-.15*alpha2,.15*alpha2)
epsilon_alpha3<-runif(nbg3,-.15*alpha3,.15*alpha3)
epsilon_alpha4<-runif(nbg4,-.15*alpha4,.15*alpha4)

epsilon_beta1<-runif(nbg1,-.15*beta1,.15*beta1)
epsilon_beta2<-runif(nbg2,-.15*beta2,.15*beta2)
epsilon_beta3<-runif(nbg3,-.15*beta3,.15*beta3)
epsilon_beta4<-runif(nbg4,-.15*beta4,.15*beta4)

#Save as data frame instead of data table
rls.pop$Data<-as.data.frame(rls.pop$Data)

#Create empty myX which will be filled in by block group below
rls.pop$Data$myX<-NA

#Strata 1
for(l in 1:nbg1){
  assign(paste("alpha1",l, sep = "_"),alpha1+epsilon_alpha1[l])
  assign(paste("beta1",l, sep = "_"),beta1+epsilon_beta1[l])
  
  if (covardist=="Gamma") {
    rls.pop$Data[rls.pop$Data$strat==1 & rls.pop$Data$BGid==l,]$myX = rgamma(tablestrat1[l,], 
                                                                             shape=get(paste("alpha1",l, sep = "_")), 
                                                                             scale = get(paste("beta1",l, sep = "_")))
  }else if (covardist=="Normal") {
    rls.pop$Data[rls.pop$Data$strat==1 & rls.pop$Data$BGid==l,]$myX = rnorm(tablestrat1[l,], 
                                  get(paste("alpha1",l, sep = "_")), 
                                  sqrt(get(paste("beta1",l, sep = "_"))))
    
  }else (NA)
  

  
}

#Strata 2
for(l in 1:nbg2){
  assign(paste("alpha2",l, sep = "_"),alpha2+epsilon_alpha2[l])
  assign(paste("beta2",l, sep = "_"),beta2+epsilon_beta2[l])
  
  BGidnum<-nbg1+l

  if (covardist=="Gamma") {
    rls.pop$Data[rls.pop$Data$strat==2 & rls.pop$Data$BGid==BGidnum,]$myX = rgamma(tablestrat2[l,], 
                                                                                   shape=get(paste("alpha2",l, sep = "_")), 
                                                                                   scale = get(paste("beta2",l, sep = "_")))
  }else if (covardist=="Normal") {
    rls.pop$Data[rls.pop$Data$strat==2 & rls.pop$Data$BGid==BGidnum,]$myX = rnorm(tablestrat2[l,], 
                                                                get(paste("alpha2",l, sep = "_")), 
                                                                sqrt(get(paste("beta2",l, sep = "_"))))
    
  }else (NA)
  
  
  
}

#Strata 3
for(l in 1:nbg3){
  assign(paste("alpha3",l, sep = "_"),alpha3+epsilon_alpha3[l])
  assign(paste("beta3",l, sep = "_"),beta3+epsilon_beta3[l])
  
  BGidnum<-nbg1+nbg2+l
  
  if (covardist=="Gamma") {
    rls.pop$Data[rls.pop$Data$strat==3 & rls.pop$Data$BGid==BGidnum,]$myX = rgamma(tablestrat3[l,], 
                                                                                   shape=get(paste("alpha3",l, sep = "_")), 
                                                                                   scale = get(paste("beta3",l, sep = "_")))
    
  }else if (covardist=="Normal") {
    rls.pop$Data[rls.pop$Data$strat==3 & rls.pop$Data$BGid==BGidnum,]$myX = rnorm(tablestrat3[l,], 
                                                                                   get(paste("alpha3",l, sep = "_")), 
                                                                                   sqrt(get(paste("beta3",l, sep = "_"))))
    
  }else (NA)
  
  
}

#Strata 4
for(l in 1:nbg4){
  assign(paste("alpha4",l, sep = "_"),alpha4+epsilon_alpha4[l])
  assign(paste("beta4",l, sep = "_"),beta4+epsilon_beta4[l])
  
  BGidnum<-nbg1+nbg2+nbg3+l
  
  if (covardist=="Gamma") {
    rls.pop$Data[rls.pop$Data$strat==4 & rls.pop$Data$BGid==BGidnum,]$myX = rgamma(tablestrat4[l,], 
                                                                                   shape=get(paste("alpha4",l, sep = "_")), 
                                                                                   scale = get(paste("beta4",l, sep = "_")))
    
  }else if (covardist=="Normal") {
    rls.pop$Data[rls.pop$Data$strat==4 & rls.pop$Data$BGid==BGidnum,]$myX = rnorm(tablestrat4[l,], 
                                                                                  get(paste("alpha4",l, sep = "_")), 
                                                                                  sqrt(get(paste("beta4",l, sep = "_"))))
    
  }else (NA)
  
}


#Summarize and plot to see differences 
summary(rls.pop$Data[rls.pop$Data$strat==1,]$myX)
summary(rls.pop$Data[rls.pop$Data$strat==2,]$myX)
summary(rls.pop$Data[rls.pop$Data$strat==3,]$myX)
summary(rls.pop$Data[rls.pop$Data$strat==4,]$myX)

ggplot(rls.pop$Data,aes(x=myX,group=strat,fill=strat))+
  geom_histogram(position="dodge",bins=10)+theme_bw()


### Centering continous naive nutrients
rls.pop$Data$c_ln_na_avg = as.numeric(scale(rls.pop$Data$ln_na_avg,scale = F))
rls.pop$Data$c_ln_k_avg = as.numeric(scale(rls.pop$Data$ln_k_avg,scale = F))
rls.pop$Data$c_ln_kcal_avg = as.numeric(scale(rls.pop$Data$ln_kcal_avg,scale = F))
rls.pop$Data$c_ln_protein_avg = as.numeric(scale(rls.pop$Data$ln_protein_avg,scale = F))

rls.pop$Data$c_age = as.numeric(scale(rls.pop$Data$age,scale = F))
rls.pop$Data$c_bmi = as.numeric(scale(rls.pop$Data$bmi,scale = F))
rls.pop$Data$female = 1*(rls.pop$Data$sex=='F')
rls.pop$Data$bkg_pr = 1*(rls.pop$Data$bkg=='PR')
rls.pop$Data$bkg_o = 1*(rls.pop$Data$bkg=='O')

### Saving the data
dt.pop <- rls.pop$Data

#Now, simulate outcome. Consider just age and BMI for now.
#Data is going to be in long form: creates 4 IDs per subject 1-1000
ID <- rep(1:ALL_N, each = ntest)

time <- rep(testtimes, times = ALL_N)

#Formula for true lambda
lambda<-blambda*exp((dt.pop$myX*mybeta))

#Uses exponential distribution to generate time to event of interest (event time = ET)
ET <- rexp(ALL_N, lambda)

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

data1<-data.frame(ID,time,result,true_result)
keep<-unlist(tapply(data1$result,data1$ID,after_first_pos))
#unlist - unlist a list of vectors into a single vector
#tapply applies my "after_first_pos" function to the "result" of data vector
datafinal_1<-data1[keep,]

data.table::setDT(data1)
data_wide_ALL<-dcast(melt(data1, id.vars=c("ID", "time")), ID~variable+time)
data_wide_ALL$GS_vis4<-data_wide_ALL$true_result_4
prop.table(table(data_wide_ALL$GS_vis4))

dt.pop<-cbind(dt.pop,data_wide_ALL)

### Now loading model coefficients and variance parameter for simulation of true and biomarker intake ###
load('U://calibrCoeff.RData')

dt.pop$ln_na_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_na_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%sodicoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[1])))
dt.pop$ln_k_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_k_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%potacoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))
dt.pop$ln_kcal_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_kcal_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%kcalcoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))
dt.pop$ln_protein_true = as.numeric(as.matrix(cbind('Int'=1,subset(dt.pop,select=c('c_ln_protein_avg','c_age','c_bmi','female','usborn','high_chol','bkg_pr','bkg_o'))))%*%protcoeff$Estimate+rnorm(nrow(dt.pop),0,sqrt(var.df$Var.True[2])))

# Centering the true nutrients

dt.pop$c_ln_na_true = as.numeric(scale(dt.pop$ln_na_true,scale=F))
dt.pop$c_ln_k_true = as.numeric(scale(dt.pop$ln_k_true,scale=F))
dt.pop$c_ln_kcal_true = as.numeric(scale(dt.pop$ln_kcal_true,scale=F))
dt.pop$c_ln_protein_true = as.numeric(scale(dt.pop$ln_protein_true,scale=F))

#Checking correlations:
with(dt.pop,cor(c_ln_na_true,c_ln_na_avg))
with(dt.pop,cor(c_ln_k_true,c_ln_k_avg))
with(dt.pop,cor(c_ln_kcal_true,c_ln_kcal_avg))
with(dt.pop,cor(c_ln_protein_true,c_ln_protein_avg))

### Now, simulating biomarker
dt.pop$ln_na_bio1 = dt.pop$ln_na_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[1]))
dt.pop$ln_na_bio2 = dt.pop$ln_na_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[1]))

dt.pop$ln_k_bio1 = dt.pop$ln_k_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[2]))
dt.pop$ln_k_bio2 = dt.pop$ln_k_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[2]))

dt.pop$ln_kcal_bio1 = dt.pop$ln_kcal_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[3]))
dt.pop$ln_kcal_bio2 = dt.pop$ln_kcal_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[3]))

dt.pop$ln_protein_bio1 = dt.pop$ln_protein_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[4]))
dt.pop$ln_protein_bio2 = dt.pop$ln_protein_true + rnorm(nrow(dt.pop),0,sqrt(var.df$Biomarker[4]))

### Centering the biomarkers
dt.pop$c_ln_na_bio1 = as.numeric(scale(dt.pop$ln_na_bio1,scale = F))
dt.pop$c_ln_na_bio2 = as.numeric(scale(dt.pop$ln_na_bio2,scale = F))

dt.pop$c_ln_k_bio1 =  as.numeric(scale(dt.pop$ln_k_bio1,scale = F))
dt.pop$c_ln_k_bio2 = as.numeric(scale(dt.pop$ln_k_bio2,scale = F))

dt.pop$c_ln_kcal_bio1 = as.numeric(scale(dt.pop$ln_kcal_bio1,scale = F))
dt.pop$c_ln_kcal_bio2 = as.numeric(scale(dt.pop$ln_kcal_bio2,scale = F))

dt.pop$c_ln_protein_bio1 = as.numeric(scale(dt.pop$ln_protein_bio1,scale = F))
dt.pop$c_ln_protein_bio2 = as.numeric(scale(dt.pop$ln_protein_bio2,scale = F))

### Now, I will simulate two outcomes (continuous sbp and hypertension2_indicator)
load('U://outPar.RData')

### Now, simulating hypertension outcome ###
suboutcome.par = data.frame(Variable = c('Intercept',all.vars(m1.formula[-2])), Estimate = m1.coefficients)
# Changing the main effect size: 1.3 OR for a 20% increase (1.2) in Sodium. From Prentice (2017)
suboutcome.par$Estimate[suboutcome.par$Variable=='C_LN_NA_CALIBR'] = log(1.3)/log(1.2)
# Adjusting the intercept for an overall prevalence of 0.08 (from real data)
suboutcome.par$Estimate[suboutcome.par$Variable=='Intercept'] = log(0.08/(1-0.08))-sum(suboutcome.par$Estimate[-1]*apply(subset(dt.pop,select=c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')),2,mean))

suboutcome.par = suboutcome.par[match(c("Intercept","C_AGE","C_BMI","C_LN_NA_CALIBR",'HIGH_TOTAL_CHOL','US_BORN','FEMALE','BKGRD_PR','BKGRD_OTHER'),suboutcome.par$Variable),]
p.par = expit(as.matrix(cbind('Intercept'=1,dt.pop[,c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')]))%*%suboutcome.par$Estimate)
dt.pop$hypertension = as.numeric(1*(runif(nrow(p.par))<=p.par))

htn_parameters <- suboutcome.par
if(!dir.exists('./Output')){system('mkdir Output')}
save(htn_parameters,file = 'U://htn_parameters.RData')
rm(htn_parameters)

### Now, simulating SBP outcome ###
# Based on real data, 1 unit increase in log_NA (or exp(1)~2.7 times increase in sodium)
# leads to an increment of 1.4058215 units in SBP. Very small effect.
# I will make the effect of sodium larger.
# I will do:
# leads to an increment of 5 units in SBP. Large effect

suboutcome.par = data.frame(Variable = c('Intercept',all.vars(m2.formula[-2])), Estimate = m2.coefficients)
# Increasing the main effect size: 0.182 unit increase in log_NA (or exp(0.182)~1.2 times increase in sodium)
suboutcome.par$Estimate[suboutcome.par$Variable=='C_LN_NA_CALIBR'] = 5/0.182
# Adjusting the intercept for an overall average of 120 SBP
suboutcome.par$Estimate[suboutcome.par$Variable=='Intercept'] = 120-sum(suboutcome.par$Estimate[-1]*apply(subset(dt.pop,select=c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')),2,mean))

suboutcome.par = suboutcome.par[match(c("Intercept","C_AGE","C_BMI","C_LN_NA_CALIBR",'HIGH_TOTAL_CHOL','US_BORN','FEMALE','BKGRD_PR','BKGRD_OTHER'),suboutcome.par$Variable),]
xb.par = as.matrix(cbind('Intercept'=1,dt.pop[,c('c_age','c_bmi','c_ln_na_true','high_chol','usborn','female','bkg_pr','bkg_o')]))%*%suboutcome.par$Estimate
dt.pop$sbp = as.numeric((xb.par+rnorm(nrow(xb.par),0,m2.sigma)))

sbp_parameters <- list()
sbp_parameters[['coefficients']] <- suboutcome.par
sbp_parameters[['variance']] <- m2.sigma^2
if(!dir.exists('./Output')){system('mkdir Output')}
save(sbp_parameters,file = 'U://sbp_parameters.RData')
rm(sbp_parameters)

### Organizing, Centering and Saving data
dt.pop$bkg %<>% as.character()

pop1 = dt.pop
pop1$v.num=1

pop2 = dt.pop
pop2$v.num=2

pop = rbind(pop1,pop2)



pop = pop[order(pop$strat,pop$BGid,pop$hhid,pop$subid,pop$v.num),]
save(pop,file=paste0('U://TargetPopulationData',covardist,'_CR',myCR,'.RData'),compress = 'xz')

# (3) Finally, we will generate 1000 samples from the target population using a stratified three-stage sampling scheme. 
# Sampling probability weights and design variables are included in the final data sets.

S=1000
set.seed(20190414)

### Generate Samples for Simulation Study ###
sampname = 'SampleData'

####### read-in population datasets #######
load(paste0('U:/',popname,".RData"))

####### function to generate sample #######
samp.hh.bystrat=function(s,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other){
  hh.select.s=rep(FALSE,dim(pop.unq[pop.unq$strat==s,])[1])         #initially set hh.select.s=FALSE for all subjects in stratum s, so that unselected BG's will be set to FALSE
  for (j in bg.select[bg.select.s==s]){
    hh.list.hisp=hh.list$hhid[hh.list$BGid==j & hisp.strat==TRUE]   #list of unique HHs in BGid j, w/ Hisp surname
    hh.list.other=hh.list$hhid[hh.list$BGid==j & hisp.strat==FALSE] #list of unique HHs in BGid j, w/ other surname
    hh.select.hisp=sample(hh.list.hisp,round(prob.hh.hisp[s]*length(hh.list.hisp)))      #list of randomly sampled HHs in BGid w/ Hispanic surname
    hh.select.other=sample(hh.list.other,round(prob.hh.other[s]*length(hh.list.other)))  #list of randomly sampled HHs in BGid w/ other surname
    hh.select.s[pop.unq$hhid[pop.unq$strat==s] %in% c(hh.select.hisp,hh.select.other)]=TRUE  #select subjects from randomly sampled HHs in BGid j (one entry per subject in BGid j)
  }
  return(hh.select.s)
}

samp.gen = function(pop,prob.bg,num.bg,hisp.prop,other.prop,prob.hh.hisp,prob.hh.other,prob.age){
  bgstrat = cumsum(num.bg)
  pop.unq=pop[pop$v.num==1,]   #this is the population only including visit 1 records (to simplify sampling)
  
  ### re-create hh.size ###
  hh.list=unique(pop.unq[,c("BGid","hhid")])  # list of BGid & hhid for each unique hhid
  bg.size=table(hh.list$BGid)   # number of HHs per BG,
  hh.size=table(pop.unq$hhid)   # number of subjects per HH
  num.hisp.strat=round(c(bg.size[as.numeric(names(bg.size))<=bgstrat[1]]*hisp.prop[1]/(hisp.prop[1]+other.prop[1]),bg.size[as.numeric(names(bg.size))>=(bgstrat[1]+1) & as.numeric(names(bg.size))<=bgstrat[2]]*hisp.prop[2]/(hisp.prop[2]+other.prop[2]),bg.size[as.numeric(names(bg.size))>=(bgstrat[2]+1) & as.numeric(names(bg.size))<=bgstrat[3]]*hisp.prop[3]/(hisp.prop[3]+other.prop[3]),bg.size[as.numeric(names(bg.size))>=(bgstrat[3]+1)]*hisp.prop[4]/(hisp.prop[4]+other.prop[4])))
  num.other.strat=bg.size-num.hisp.strat
  A=matrix(rep(NA,times=max(bg.size)*length(bg.size)),nrow=max(bg.size),ncol=length(bg.size))
  for (i in 1:length(bg.size)){
    A[,i]=c(rep(TRUE,times=num.hisp.strat[i]),rep(FALSE,times=num.other.strat[i]),rep(NA,times=max(bg.size)-bg.size[i]))    # create matrix A with each column corresponding to a BG (containing 1's for each Hispanic HH followed by 0's for each other HH)
  }
  hisp.strat=na.omit(c(A))   #indicator for Hispanic surname (one entry per HH)
  
  ### generate raw weights (these do not depend on the sample (only depend on sampling probabilities)) ###
  pop$W.bg=1/prob.bg[pop$strat]   # BG stage raw weight (based on BG sampling fraction; 1/sampling fraction for that stratum)
  pop$W.hh=ifelse(pop$hhid %in% hh.list$hhid[hisp.strat],1/prob.hh.hisp[pop$strat],1/prob.hh.other[pop$strat])   #HH stage raw weight
  pop$W.sub=ifelse(pop$age.strat,1/prob.age[2],1/prob.age[1])    # subject stage raw weight
  
  ### select random sample from population ###
  #select stratified random sample of BGs from pop & save list of BGs
  bg.select=c(sample(1:bgstrat[1],round(prob.bg[1]*num.bg[1])),sample((bgstrat[1]+1):bgstrat[2],round(prob.bg[2]*num.bg[2])),sample((bgstrat[2]+1):bgstrat[3],round(prob.bg[3]*num.bg[3])),sample((bgstrat[3]+1):bgstrat[4],round(prob.bg[4]*num.bg[4])))
  bg.select.s=1*(bg.select<=bgstrat[1])+2*(bg.select>=(bgstrat[1]+1) & bg.select<=bgstrat[2])+3*(bg.select>=(bgstrat[2]+1) & bg.select<=bgstrat[3])+4*(bg.select>=(bgstrat[3]+1) & bg.select<=bgstrat[4])  #stratum for each BG
  
  #select stratified random sample of HHs from each selected BG & save indicator of HH selection for all subjects within selected HHs
  hh.select=c(samp.hh.bystrat(1,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),samp.hh.bystrat(2,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),
              samp.hh.bystrat(3,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other),samp.hh.bystrat(4,pop.unq,bg.select,bg.select.s,hh.list,hisp.strat,prob.hh.hisp,prob.hh.other))  #indicator of HH selection
  
  #select random sample of subjects from each selected HH & save indicator of subject selection
  sub.select=rep(FALSE,dim(pop.unq)[1])           # initially set sub.select=FALSE for all subjects, so that unselected HH's will be set to FALSE
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[!pop.unq$age.strat & hh.select],round(prob.age[1]*dim(pop.unq[!pop.unq$age.strat & hh.select,])[1]))]=TRUE   # randomly sample younger subjects among sampled HH's
  sub.select[pop.unq$subid %in% sample(pop.unq$subid[pop.unq$age.strat & hh.select],round(prob.age[2]*dim(pop.unq[pop.unq$age.strat & hh.select,])[1]))]=TRUE     # randomly sample older subjects among sampled HH's
  
  samp=pop[rep(sub.select,each=2),]   # create sample by restricting pop to selected subjects (replicate sub.select vector to select both visits from selected subjects)
  
  # generate normalized weight
  samp$W.bghhsub=samp$W.bg*samp$W.hh*samp$W.sub     # raw combined weights
  samp$bghhsub_s2=samp$W.bghhsub/mean(samp$W.bghhsub[samp$v.num==1])    # normalized weight (raw combined weight/mean combined weight)
  
  #samp=samp[, colnames(samp)!="hhid"]
  
  return(samp)
}

my_prop<-meanN/3944
myNAns<-NULL

####### generate S samples from same population for each scenario #######
for(i in 1:S){
  
  samp=samp.gen(pop,
                prob.bg=c(.25,.25,.6,.6),
                num.bg=c(58,21,130,167),
                hisp.prop=c(.2,.2125,.2975,.2975),
                other.prop=c(.15,.225,.26,.2925),
                prob.hh.hisp=c(.18,.225,.14,.14)*my_prop,
                prob.hh.other=c(.025,.0175,.035,.04)*my_prop,
                prob.age=c(.3575,.55))
  
  
  samp=subset(samp,v.num==1)
  samp$dat.num=i
  samp$solnas = F; samp$solnas[sample(x=nrow(samp),size=450,replace=F)] = T
  
  ### Sorting
  samp = samp[order(samp$strat,samp$BGid,samp$hhid,samp$subid),]
  
  #Save Sample Size
  
  saveN<-dim(samp)[1]
  myNAns<-rbind(myNAns,saveN)
  
  save(samp,file=paste0('U:/',popname,'_N',meanN,'/',sampname,"_",str_pad(i,nchar(S),pad=0),".RData"),compress = 'xz')
  
  print(i)
}



