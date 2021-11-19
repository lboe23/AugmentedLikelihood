#LOG LIKELIHOOD
log_like_proposed<-function(params,nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights,purpose)
{
  beta_func=params[1:nbeta]
  
  for (i in 1:J){
    assign(paste("gam",i, sep = ""), params[i+nbeta])}
  
  #Create a vector of all of the gam 1-J
  allgam<-NULL
  for (i in 1:J){
    allgam<-(cbind(allgam,get(paste("gam",i, sep = ""))))}
  
  #Create a cumulative sum of all the gams vector
  sumgams<-cumsum(allgam)
  
  #now, S1=1, but all other S's are the cumulative sums to the exp(-exp())
  S1<-1
  
  #Create function for exp(-exp())
  EnegE<-function(x){
    exp(-exp(x))
  }
  
  #Assign all other S's
  for (i in 2:(J+1)){
    assign(paste("S",i, sep = ""), EnegE(sumgams[i-1]))
  }
  
  u<-(Xmat)%*%beta_func
  
  l<-vector(length=nsub)
  for(i in c(1:nsub)) {
    xibeta<-u[i]
    wi<-weights[i]
    Di<-Dm[i,]
    Ci<-Cm[i,]
    
    #If missing, use method from Gu et al.
    if(is.na(GSdelta[i])){
      #Now, assign Svec for missing GS
      Svec<-NULL
      for (k in 1:J){
        Svec<-(cbind(Svec,(get(paste("S",k, sep = ""))^exp(xibeta)-get(paste("S",k+1, sep = ""))^exp(xibeta))))}
      
      Svec<-cbind(Svec,get(paste("S",J+1, sep = ""))^exp(xibeta))
      Svec<-as.vector(Svec)
  
    }else if (GSdelta[i]==0){
      #Now, assign Svec for GS = 0

      SvecTEST<-NULL
      for (k in 1:J){
        SvecTEST<-(cbind(SvecTEST,(get(paste("S",k, sep = ""))^exp(xibeta)-get(paste("S",k+1, sep = ""))^exp(xibeta))))}
      
      SvecTEST<-cbind(SvecTEST,get(paste("S",J+1, sep = ""))^exp(xibeta))
      SvecTEST<-as.vector(SvecTEST)
      
   
      Svec<-c(rep(0,GSVis[i]),SvecTEST[(GSVis[i]+1):length(SvecTEST)])
      
      
      
    }else if (GSdelta[i]==1){
      #Now, assign Svec for GS = 1
      SvecTEST<-NULL
      for (k in 1:J){
        SvecTEST<-(cbind(SvecTEST,(get(paste("S",k, sep = ""))^exp(xibeta)-get(paste("S",k+1, sep = ""))^exp(xibeta))))}
      
      SvecTEST<-cbind(SvecTEST,get(paste("S",J+1, sep = ""))^exp(xibeta))
      SvecTEST<-as.vector(SvecTEST)
      
      Svec<-c(SvecTEST[1:GSVis[i]],rep(0,(length(SvecTEST)-GSVis[i])))
      
    } else{
      Svec<-0
    }
    Li<-t(Svec)%*%Ci
    l[i]<-wi*log(Li)
    
  }
  
  if (purpose=="SUM"){
    return(-(sum(l)))
    
  }else if (purpose=="INDIVIDUAL"){
    return(-(l))
  }
  else{
    return("STATE PURPOSE OF FUNCTION: SUM or INDIVIDUAL")
  }
}

#GRADIENT/ESTIMATING FUNCTION
gradient_proposed<-function(params,nsub,J,nbeta,Dm,Cm,Xmat,GSdelta,GSVis,weights,purpose)
{
  beta_func=params[1:nbeta]
  
  for (i in 1:J){
    assign(paste("gam",i, sep = ""), params[i+nbeta])}
  
  #Create a vector of all of the gam 1-J
  allgam<-NULL
  for (i in 1:J){
    allgam<-(cbind(allgam,get(paste("gam",i, sep = ""))))}
  
  #Create a cumulative sum of all the gams vector
  sumgams<-cumsum(allgam)
  
  #now, S1=1, but all other S's are the cumulative sums to the exp(-exp())
  S1<-1
  
  #Create function for exp(-exp())
  EnegE<-function(x){
    exp(-exp(x))
  }
  
  #Assign all other S's
  for (i in 2:(J+1)){
    assign(paste("S",i, sep = ""), EnegE(sumgams[i-1]))
  }

  u<-(Xmat)%*%beta_func

  for (l in 1:nbeta){
    assign(paste("g",l, sep = ""), vector(length=nsub))}
  
  for (i in 2:(J+1)){
    assign(paste("S",i,"grad", sep = ""), vector(length=nsub))
  }
  
  for(i in c(1:nsub)) {
    xibeta<-u[i]
    wi<-weights[i]
    
    for (h in 1:nbeta){
      assign(paste("x",h,"i", sep = ""), Xmat[i,h])}
    
    Di<-Dm[i,]
    Ci<-Cm[i,]
    
    GAM.x <- sumgams+xibeta
  
    #Now, assign Svec for missing GS
    SvecAll<-NULL
    for (k in 1:J){
      SvecAll<-(cbind(SvecAll,(get(paste("S",k, sep = ""))^exp(xibeta))*log(get(paste("S",k, sep = "")))
                   -(get(paste("S",k+1, sep = ""))^exp(xibeta))*log(get(paste("S",k+1, sep = "")))))}
    SvecAll<-cbind(SvecAll,(get(paste("S",J+1, sep = ""))^exp(xibeta))*log(get(paste("S",J+1, sep = ""))))
    SvecAll<-as.vector(SvecAll)
    
    #Now, assign Svec for missing GS
    SdenomAll<-NULL
    for (k in 1:J){
      SdenomAll<-(cbind(SdenomAll,(get(paste("S",k, sep = ""))^exp(xibeta)-get(paste("S",k+1, sep = ""))^exp(xibeta))))}
    
    SdenomAll<-cbind(SdenomAll,get(paste("S",J+1, sep = ""))^exp(xibeta))
    SdenomAll<-as.vector(SdenomAll)
    
    #Helper function for SGradNum
    Sadditionfunc<-function(h){
      S2GradNum_vec<-NULL
      
      if (h<J){
        for (k in (h+1):J){
          S2GradNum_vec<-cbind(S2GradNum_vec,Ci[k]*(exp(-exp(GAM.x[k])+GAM.x[k])-exp(-exp(GAM.x[k-1])+GAM.x[k-1])))
        }
        S2GradNum_vec2<-c(Ci[h]*(exp(-exp(GAM.x[h])+GAM.x[h])),S2GradNum_vec,-Ci[J+1]*(exp(-exp(GAM.x[J])+GAM.x[J])))
      }else if (h==J) {

        S2GradNum_vec2<-c(Ci[h]*(exp(-exp(GAM.x[h])+GAM.x[h])),S2GradNum_vec,-Ci[J+1]*(exp(-exp(GAM.x[J])+GAM.x[J])))
        
      }
      return((S2GradNum_vec2))
      
    }
    
     
    for (l in 2:(J+1)){
      assign(paste("S",l,"GradNumVecAll",sep = ""), Sadditionfunc(l-1))
    }
    
    
   
    S2GradDen_vec<-NULL
    for (h in 2:J){
      S2GradDen_vec<-cbind(S2GradDen_vec,Ci[h]*(exp(-exp(GAM.x[h-1]))-exp(-exp(GAM.x[h]))))
    }
    S2GradDen_vecAll<-(cbind(Ci[1]*(S1^exp(xibeta)-exp(-exp(GAM.x[1]))),S2GradDen_vec,Ci[J+1]*(exp(-exp(GAM.x[J])))))
    
#GRADIENT FOR GSdelta NA    
    if(is.na(GSdelta[i])){
      Svec<-SvecAll
      Sdenom<-SdenomAll
      S2GradDen<-sum(S2GradDen_vecAll)
      for (l in 2:(J+1)){
        assign(paste("S",l,"GradNum",sep = ""),  sum(get(paste("S",l,"GradNumVecAll",sep = ""))))
      }

#GRADIENT FOR GSdelta =0   
    }else if (GSdelta[i]==0){
      
      Svec<-c(rep(0,GSVis[i]),SvecAll[(GSVis[i]+1):length(SvecAll)])
      Sdenom<-c(rep(0,GSVis[i]),SdenomAll[(GSVis[i]+1):length(SdenomAll)])
      S2GradDen<-sum(c(rep(0,GSVis[i]),S2GradDen_vecAll[(GSVis[i]+1):length(S2GradDen_vecAll)]))
       
      for (l in 2:(J+1)){
          if ((GSVis[i]+1)>=l){
          assign(paste("S",l,"GradNum",sep = ""),  sum(get(paste("S",GSVis[i]+1,"GradNumVecAll",sep = ""))[-1]) )
          }else if((GSVis[i]+1)<l){
            assign(paste("S",l,"GradNum",sep = ""),  sum(get(paste("S",l,"GradNumVecAll",sep = ""))) )
            
            }
        }

#GRADIENT FOR GSdelta =1
    }else if (GSdelta[i]==1){
      Svec<-c(SvecAll[1:GSVis[i]],rep(0,(length(SvecAll)-GSVis[i])))
      Sdenom<-c(SdenomAll[1:GSVis[i]],rep(0,(length(SdenomAll)-GSVis[i])))
      S2GradDen<-sum(c(S2GradDen_vecAll[1:GSVis[i]],rep(0,(length(S2GradDen_vecAll)-GSVis[i]))))
      
      for (l in 2:(J+1)){
        assign(paste("S",l,"GradNumVecAllLength",sep = ""), c(rep(0,J+1-length( get(paste("S",l,"GradNumVecAll",sep = "")))), get(paste("S",l,"GradNumVecAll",sep = ""))))
      }
      
      for (l in 2:(J+1)){
        assign(paste("S",l,"GradNum",sep = ""), sum(get(paste("S",l,"GradNumVecAllLength",sep = ""))[1:GSVis[i]]))
      }
      
         
    } else{
      Svec<-0
      Sdenom<-0
      S2GradNum<-0
    }
    Gnum<-t(Svec)%*%Ci
    Gdenom<-t(Sdenom)%*%Ci

    #Assign g and Sgrad pieces
    for (k in 1:nbeta){
      assign(paste0("g", k), `[<-` (get(paste0("g", k)), i, (wi)*(Gnum*get(paste0("x",k,"i"))*exp(xibeta))/Gdenom ))}
    
    for (k in 2:(J+1)){
      assign(paste0("S", k,"grad"), `[<-` (get(paste0("S", k,"grad")), i, (wi)*(get(paste0("S", k,"GradNum")))/S2GradDen ))}
    
   
    
    
  }
  
  #Now, outside loop, sum up gs for gradient
  
  for (k in 1:nbeta){
    assign(paste0("sum_g", k), -sum(get(paste0("g", k))) )
  }
  
  sums<-NULL
  for (k in 1:nbeta){
    sums<-cbind(sums,get(paste0("sum_g", k)))}
  
  #Sum up S's for gradient
  for (k in 2:(J+1)){
    assign(paste0("sum_S", k), -sum(get(paste0("S", k,"grad"))) )
  }
  
  sums_S<-NULL
  for (k in 2:(J+1)){
    sums_S<-cbind(sums_S,get(paste0("sum_S", k)))}
  
  #If returning each person's estimating function contribution, need to create matrix
  all_g<-NULL
  for (k in 1:nbeta){
    all_g<-cbind(all_g,get(paste0("g", k)))}
  
  all_S<-NULL
  for (k in 2:(J+1)){
    all_S<-cbind(all_S,get(paste0("S", k,"grad")))}
  
  ncols=J+nbeta
  
  if (purpose=="SUM"){
    return(c(sums,sums_S))
    
  }else if (purpose=="INDIVIDUAL"){
    return(matrix(c(all_g,all_S),ncol=ncols,byrow=FALSE))
  }
  else{
    return("STATE PURPOSE OF FUNCTION: SUM or INDIVIDUAL")
  }
  
 
}




