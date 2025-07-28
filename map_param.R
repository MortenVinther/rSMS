
map_param<-function(data,parameters) {
  
  my.map<-lapply(parameters,function(x) {if (length(x)>0) a<-factor(1L:length(x)) else a<-factor(numeric(0)); a})
  
  # adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
  if (data$zeroCatchYearExists==1 | length(parameters$logSeparAgeF)>0 ) {
    if (any(data$nlogF>0)) {
      UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
      for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
        zy<-data$zeroCatchYear[[s]]
        fromTo<-data$nlogFfromTo[s,]
        UfMap[fromTo[1]:fromTo[2],zy]<-NA
        parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
      }
      if (length(parameters$logSeparAgeF)>0) for (s in 1:data$nSpecies){
        fSepar<-data$info[s,'fSepar']
        if (fSepar<90) {
          ii<- max(data$nlogFfromTo[s,]) # last random walk group for the species is used as a year effect in separable F model
          #cat('ii:',ii,'\n')
          parameters$Uf[ii,1]<-0.0  
          UfMap[ii,1]<-NA
        }
      }
      UfMap<-factor(UfMap)
    }
  }
 
 
  if ((data$zeroCatchYearExists==1 | length(parameters$logSeparAgeF)>0) & (any(data$nlogF>0)))  my.map<-purrr::assign_in(my.map,'Uf',UfMap)

  
  # my.map<-list(Uf=UfMap) else my.map=list()
 
  aMap<-1L:length(parameters$rec_loga)
  bMap<-1L:length(parameters$rec_logb)
  aMap[data$stockRecruitmentModelCode==0]<-NA
  bMap[data$stockRecruitmentModelCode==0]<-NA
  bMap[data$stockRecruitmentModelCode==3]<-NA
  bMap[data$stockRecruitmentModelCode==6]<-NA
  noRec<-data$inclSsbR==0  & !data$doProcessN_any
  aMap[noRec]<-NA
  bMap[noRec]<-NA
  aMap<-factor(aMap)
  bMap<-factor(bMap)
  #my.map<-c(my.map,list(rec_loga=aMap),list(rec_logb=bMap))
  my.map<-purrr::assign_in(my.map,'rec_loga',aMap)
  my.map<-purrr::assign_in(my.map,'rec_logb',bMap)

  
  if (any(data$useRho<=0)) {
    rhoMap<-1L:length(data$useRho)
    rhoMap[data$useRho==0]<-NA
    rhoMap<-factor(rhoMap)
    my.map<-purrr::assign_in(my.map,'rho',rhoMap)
  }
  
  #FfromSeparableModel
  if (any(data$FfromSeparableModel)) {
    # set the last quarter (used) to NA (to fix this parameter to 1.0)
    key<-data$keylogSeasonF
    if (data$nSpecies==1) spl<-1 else spl<-c(1,2)
    all0<-apply(data$zeroCatchSeasonAge,MARGIN=spl,function(x) all(x==0))
    if (length(spl)==1) all0<-matrix(all0,c(data$nSpecies,length(all0)))
    # ftable(all0)
    
    for (s in (1:data$nSpecies)) if (nrow(key[[s]])>0) {
      for (q in (1:dim(all0)[[2]])) if (all0[s,q] ) key[[s]][,q,]<-0L
      
      maxQ<-unique(apply(key[[s]],c(1,3),max))
      key[[s]][key[[s]] %in%  maxQ]<-0L
      key[[s]][,1:(data$fqa[1]-1),data$fa+data$off.age] <- 0L
    } 
    incl<-unlist(lapply(key,unique))
    incl<-sort(incl[incl>0])
    seasonMap<-1L:length(parameters$logFSeasonal)
    #incl<-c(3,4,5,6)
    seasonMap[setdiff(seasonMap,incl)]<-NA
    if (data$fixMapseasonalF[1]  >0)   seasonMap[data$fixMapseasonalF]<-NA
    seasonMap<-factor(seasonMap)
     
    my.map<-purrr::assign_in(my.map,'logFSeasonal',seasonMap)
  }
 
  
  if (data$sms.mode>0) {
    if( any(!data$usestomObsVar) ) {
      sOVMap<-1L:length(data$usestomObsVar)
      sOVMap[!data$usestomObsVar]<-NA
      sOVMap<-factor(sOVMap)
      my.map<-purrr::assign_in(my.map,'logStomObsVar',sOVMap)
    }
  
    if (length(my.map$vulnera) >0) {
      vulMap<-1:length(parameters$vulnera)
      vulNo<-data$vulneraIdx[data$otherFoodName,]
      vulNo<-vulNo[vulNo>0]
      vulMap[vulNo]<-NA
      vulMap<-factor(vulMap)
      my.map<-purrr::assign_in(my.map,'vulnera',vulMap)
    }
  }  
  # 
  if (length(parameters$logSeparAgeF)>0) {
     sepMap<-1L:length(parameters$logSeparAgeF)
     for (s in 1:data$nSpecies) {
       if (data$info[s,'fModel']==2) {
         fix<-max(unique(as.vector(data$keylogSeparAgeF[[s]])))
         sepMap[fix]<-NA
       }
     } 
        
     sepMap<-factor(sepMap)
     my.map<-purrr::assign_in(my.map,'logSeparAgeF',sepMap)
     #my.map<-c(my.map,list(logSeparAgeF=sepMap))
   }
  #grep('logSeparAgeF',nl)
  
  ###
  
  if (length(parameters$logNrecruitParam) >0) {
    
    used<-data$logNrecruitParamfromTo[ data$logNrecruitParamfromTo[,1]>0,,drop=FALSE]
    yfMap<-matrix(1L:length(parameters$logNrecruitParam),ncol=data$nYears.model,nrow=nrow(used),byrow=TRUE)
    rownames(yfMap)<-rownames(used)

    # last model year, potential without recruitment 
    yfMap[,data$nYears.model]<-yfMap[,data$nYears.model-data$extendMode]
    yfMap<-factor(t(yfMap))
    my.map<-purrr::assign_in(my.map,'logNrecruitParam',yfMap)
  }
  
  
 if (length(parameters$logYearEffectF) >0) {
   yfMap<-matrix(1L:length(parameters$logYearEffectF),ncol=data$nYears.model,nrow=max(data$idxLogYearEffectF),byrow=TRUE)
   rownames(yfMap)<-rownames(parameters$logYearEffectF)
   yfMap[,1]<-NA
   
   if (data$zeroCatchYearExists==1) {
      for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0  & data$info[s,'fModel']==3) {
       zy<-data$zeroCatchYear[[s]]
       yfMap[data$spNames[s],zy]<-NA
      }} 
   
   sp<-match(rownames(parameters$logYearEffectF),data$spNames)
   for (s in sp) {
     ll<-length(data$catchSepYear[[s]])
     if (ll>1) {
      for (y in (2:ll) ) yfMap[data$spNames[s],data$catchSepYear[[s]][y]+data$off.year ]<-NA
     }
   }
  
   # last model year, potential without catch data 
   yfMap[,data$nYears.model]<-yfMap[,data$nYears.model-data$extendMode]
   yfMap<-factor(yfMap)
   my.map<-purrr::assign_in(my.map,'logYearEffectF',yfMap)
 }
 

   #mapply(function(x,y) length(x)==length(y),parameters,my.map) 

  return(my.map)
}


lock_param<-function(data,parameters,my.map,lockP){
  for (p in lockP) {
    if (p %in% names(my.map)) {
      my.map[[p]][]<-NA
      my.map[[p]]<- factor(my.map[[p]])
    } else {
      m<-1L:length(parameters[[p]])
      m[]<-NA
      my.map[[p]]<-factor(m)
    }
  }
  return(my.map)
}

