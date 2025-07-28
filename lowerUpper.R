lowerUpper<-function(obj,data,parameters) {
    nl<-names(obj$par)
  
  spLimit<-function(nl,key,param,l=20) {
    a<-array2DF(key) %>% mutate(Var2=NULL) %>% filter(Value>=1) %>%unique() %>% arrange(Value)
    found<-nl==param
    nl[found]<-paste(nl[found],substr(a$Var1,1,l),sep='_')
    nl
  }
  
  
  nlSp<-nl
  #sort(unique(nl))
 # nlSp<-spLimit(nlSp,key=data$keyVarObsCatch,param="logSdLogObsCatch") 
  nlSp<-spLimit(nlSp,key=data$keyCatchability,param="logCatchability",l=3) 
  nlSp<-spLimit(nlSp,key=data$keyVarObsSurvey,param="logSdLogObsSurvey",l=3) 
  nlSp<-spLimit(nlSp,key=data$keyLogFstaSd,param="logSdLogFsta") 
  nlSp<-spLimit(nlSp,key=data$keyVarLogN,param="logSdLogN") 
  
  nlSp  
  
  found<-nlSp=='rho';      nlSp[found]<- paste(nlSp[found],data$spNames[data$useRho],sep='_')
  useRec<-!(data$inclSsbR==0  & !data$doProcessN_any)
  found<-nlSp=='rec_loga'; if (sum(found)>0) nlSp[found]<- paste(nlSp[found],data$spNames[data$stockRecruitmentModelCode>0][useRec],sep='_')
  found<-nlSp=='rec_logb'; if (sum(found)>0) nlSp[found]<- paste(nlSp[found],data$spNames[data$stockRecruitmentModelCode>0 & (!data$stockRecruitmentModelCode %in% c(3,6))],sep='_')
  
  found<-nlSp=='logSsbRsd';      nlSp[found]<- paste(nlSp[found],data$spNames[data$inclSsbR>0],sep='_')

  lower <- obj$par*0-Inf
  upper <- obj$par*0+Inf
  lower[nl=='rho'] <- 0.01
  upper[nl=='rho'] <- 0.99
  
  lower[nl=="logSdLogObsSurvey"]<-rep(log(data$minVarObsSurvey),length(parameters$logSdLogObsSurvey))
  upper[nl=="logSdLogObsSurvey"]<-rep(log(2.0),length(parameters$logSdLogObsSurvey))
  
  lower[nl=="logSdLogObsCatch"]<-rep(log(data$minVarObsCatch),length(parameters$logSdLogObsCatch))
  upper[nl=="logSdLogObsCatch"]<-rep(log(2.0),length(parameters$logSdLogObsCatch))
  
  lower[nl=='logSsbRsd'] <- rep(log(0.1),length(parameters$logSsbRsd))

  lower[nl=='logSdLogFsta'] <- log(0.05)
  
  if (data$sms.mode>0) if (any ( found<-nlSp=='logStomObsVar')) {
    lower[nl=="logStomObsVar"]<-log(0.1^2)
  }
  
 # process N
  if (FALSE) {
    firstGroup<-data$keyVarLogN[,1]  #recruits
    secondGroup<-unique(as.vector(data$keyVarLogN[,-1])); 
    secondGroup<-secondGroup[secondGroup>0]
    names(secondGroup)<-data$spNames
    upper[nl=="logSdLogN"][secondGroup['SSA']]<-log(0.05)
  }
  
  
  return(list(lower=lower,upper=upper))
}
