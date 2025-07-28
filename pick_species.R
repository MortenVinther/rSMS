pick_species<-function(ps=c(1L), pso=c(0L),inp, smsConf=0L) {


 #test_multi  ps<-my.ps; pso=my.pso; inp=inp_all; smsConf=1L ;ps; pso

  if (smsConf==0) pso<-0L
  ps<-sort(unique(ps))
  nps<-length(ps)
 
  pso<-sort(unique(pso))
  npso<-length(pso)
  
  if (npso==1 & pso[1]==0) npso<-0
  if (smsConf>0 & npso>0) psAll<-c(ps,pso) else psAll=ps
  
  #if (smsConf>0 ) psAll<-c(ps,pso) else psAll=ps

  
  parameters<-inp[["parameters"]]
  data<-inp[["data"]]
  
  p<-parameters
  d<-data
  
  sOld<-data$info[,'s']
  sOld<-data.frame(oldName=rownames(data$info),oldNo=sOld)
  
  
  d$spNames<-  rownames(data$info[data$info[,'s'] %in% ps,,drop=FALSE] )
  d$othspNames<-rownames(data$info[data$info[,'s'] %in% pso,,drop=FALSE] )
  d$allSpNames<-c(d$spNames,d$othspNames)

  d$allSpNamesLong<-data$allSpNamesLong[ match(d$allSpNames,data$allSpNames) ]
  
  d$info<-data$info[psAll,,drop=FALSE] 
  d$info[,'s']<-1L:dim(d$info)[[1]]
  sNew<-d$info[,'s',drop=FALSE]
  sNew<-data.frame(newName=rownames(sNew),newNo=sNew)
  
  d$sOldNew<-left_join(sOld,sNew,by=join_by(oldName==newName)) %>%
    transmute(species=oldName,oldNo,newNo=s)
  
  
  d$nSpecies<-length(ps)
  d$nSpeciesAll<-length(psAll)
  npsAll<- d$nSpeciesAll
  d$preds<-rownames(d$info[d$info[,'predator'] %in% c(1,2),] )
  d$preys<-rownames(d$info[d$info[,'prey'] %in% c(1),] )
  
  d$nAges<-max(d$info[,'last-age'])-data$minAge+1 
  ages<-1:d$nAges

  
  d$fbarRange<-data$fbarRange[ps,,drop=FALSE]     
  d$useRho<-data$useRho[ps]
  d$stockRecruitmentModelCode<-data$stockRecruitmentModelCode[ps]
  d$stockRecruitmentModelCodeAdd<-data$stockRecruitmentModelCodeAdd[ps]
  d$zeroCatchYearExistsSp<-data$zeroCatchYearExistsSp[ps]
  d$zeroCatchYearExists<-any(d$zeroCatchYearExistsSp==1)  
  d$zeroCatchYear<-data$zeroCatchYear[ps]
  d$seasonalCatches<-data$seasonalCatches[ps]
  d$FfromSeparableModel<-data$FfromSeparableModel[ps]
  d$useFrandomWalk<-data$useFrandomWalk[ps]

  d$catchSepYear<-lapply(ps,function(s) data$catchSepYear[[s]])
  
  dm<- dimnames(data$zeroCatchSeasonAge)
  d$zeroCatchSeasonAge<-data$zeroCatchSeasonAge[ps,,]
 # dm[[1]]<-as.character(1:d$nSpecies)
 
  # dimnames(d$zeroCatchSeasonAge)<-dm
  d$doProcessNoise<-data$doProcessNoise[ps]
  d$doProcessN_any<-data$doProcessN_any[ps]
  d$doProcessN_old<-data$doProcessN_old[ps]
  d$doProcessN_none<-data$doProcessN_none[ps]
  
  
  cut_tab<-function(tab,reNumber=FALSE,surv=FALSE) {
   if (surv) tab<-tab[data$keySurvey.overview[,"s"] %in% ps,ages,drop=FALSE]  else tab<-tab[ps,ages,drop=FALSE] 
   if (reNumber) {
     k<-sort(unique(tab[tab>0]))
     kk<-1:length(k)
     names(kk)<-k
     tab[tab>0]<-kk[as.character(tab[tab>0])]
   }
   tab
  }
  
  cutFromTo<-function(tab){
    tab<-tab[ps,,drop=FALSE]
    n<-tab[,2]-tab[,1]+1L
    tab[,2]<-cumsum(n)
    tab[,1]<-tab[,2]-n+1L
    tab
  }
 
  cutFromTo0<-function(tab){
   # tab<-data$logNfirstYparamfromTo
    tab<-tab[ps,,drop=FALSE]
    n<-tab[,2]-tab[,1]+1L
    n[tab[,1]==0]<-0
    tab[,2]<-cumsum(n)
    tab[n==0,2]<-0
    tab[,1]<-tab[,2]-n+1L
    tab[n==0,1]<-0
    tab
  }
  
   
  d$keyLogFsta<-cut_tab(data$keyLogFsta,reNumber=FALSE)   
  d$seasonalCatches<-data$seasonalCatches[ps]
  d$keyLogFstaSd<-cut_tab(data$keyLogFstaSd,reNumber=TRUE)

  d$nlogFfromTo<-cutFromTo0(data$nlogFfromTo)
  d$nlogNfromTo<-cutFromTo0(data$nlogNfromTo)
  d$logNfirstYparamfromTo<-cutFromTo0(data$logNfirstYparamfromTo)
  d$logNrecruitParamfromTo<-cutFromTo0(data$logNrecruitParamfromTo)

  x<-data$inclSsbR[ps]; x[x>0]<-1L; d$inclSsbR<-cumsum(x); d$inclSsbR[x==0]<-0L
  d$SsbRweight<-data$SsbRweight[ps]
  
  #d$nlogF <- d$nlogFfromTo[,2]- d$nlogFfromTo[,1]+ 1L
  d$nlogF <- d$nlogFfromTo[,2]- d$nlogFfromTo[,1]+ ifelse(apply(d$nlogFfromTo,1,sum)>0,1L,0L)
  d$nlogN <- d$nlogNfromTo[,2]- d$nlogNfromTo[,1]+ ifelse(apply(d$nlogNfromTo,1,sum)>0,1L,0L)
  
  cut_tab_3<-function(tab,reNumber=FALSE) {
    tab<-tab[ps] 
    if (reNumber) {
      x<-  do.call(c,lapply(tab,function(x) unique(as.vector(x))))
      k<-x[x>0]  
      kk<-1:length(k)
      names(kk)<-k
      tab<-lapply(tab,function(x)  {
        dims<-dim(x); 
        x[x>0]<-kk[as.character(x[x>0])]
        matrix(x,ncol=dims[[2]],dimnames=dimnames(x))
      })
    }
    
  }
  d$keylogSeparAgeF<-cut_tab_3(tab=data$keylogSeparAgeF,reNumber=TRUE) 
  d$keyVarLogN<-cut_tab(data$keyVarLogN,reNumber=TRUE)
  d$keyVarObsSurvey<-cut_tab(tab=data$keyVarObsSurvey,reNumber=TRUE,surv=TRUE)
  d$keyCatchability<-cut_tab(data$keyCatchability,reNumber=TRUE,surv=TRUE)

  cut_tab_4<-function(tab,reNumber=FALSE) {
    tab<-tab[ps] 
    if (reNumber) {
      k<-sort(unique(c(0,unlist(tab))))
      kk<-0:(length(k)-1L)
      names(kk)<-k
      tab<-map(tab,function(x) array(kk[as.character(x)],dim=dim(x),dimnames=dimnames(x)))
     }
  }
  
  d$keylogSeasonF<-cut_tab_4(tab=data$keylogSeasonF,reNumber=TRUE)
  if (FALSE){
    data$keylogSeasonF[[6]][1,,]  
    data$keylogSeasonF[[6]][25,,] 
    data$keylogSeasonF[[6]][40,,] 
    max(data$keylogSeasonF[[6]])
    
    d$keylogSeasonF[[1]][1,,]  
    d$keylogSeasonF[[1]][25,,] 
    d$keylogSeasonF[[1]][40,,] 
    max(d$keylogSeasonF[[1]])
  }  
  p$logFSeasonal<-rep(log(1/data$nSeasons),max(unlist(lapply(d$keylogSeasonF,function(x) if (nrow(x) >0)  max(x) else 0L))))
  
  
  if (dim(data$bfix)[[1]]>0) d$bfix<-filter(data$bfix,s %in% ps) %>% arrange(qGrp) %>% mutate(qGrp=1:dplyr::n()) else d$bfix<-data$bfix
  
  d$keyVarObsCatch.df<-left_join(filter(data$keyVarObsCatch.df, s %in% ps) %>% 
                                   mutate(varGroup=if_else(vari,1,0), varGroup=cumsum(varGroup)),transmute(d$sOldNew,oldNo,newNo),by = join_by(s==oldNo)) %>% mutate(s=newNo,newNo=NULL)
  
  
  if (any(d$bfix$fix)) {
    fixMap<-d$bfix$qGrp
    d$fixMapseasonalF<-b$qGrp[bfix$fix]
    p$logFSeasonal[fixMapseasonalF]<-0.0
  } else d$fixMapseasonalF<- -9
  
 
  x<-data$idxLogYearEffectF[ps]; x[x>0]<-1L;x<-cumsum(x)
  d$idxLogYearEffectF<- x 
  d$propMat<-data$propMat[ps]; names(d$propMat)<-d$spNames
  d$stockMeanWeight<-data$stockMeanWeight[ps]; names(d$stockMeanWeight)<-d$spNames
  d$catchMeanWeight<-data$catchMeanWeight[ps]; names(d$catchMeanWeight)<-d$spNames
  d$catchNumber<-data$catchNumber[ps]; names(d$catchNumber)<-d$spNames
  
  d$seasFprop<-data$seasFprop[ps]; names(d$seasFprop)<-d$spNames
  d$logSeasFprop<-data$logSeasFprop[ps]; names(d$logSeasFprop)<-d$spNames
  d$natMor<-data$natMor[ps]; names(d$natMor)<-d$spNames
  d$propF<-data$propF[ps]; names(d$propF)<-d$spNames
  d$propM<-data$propM[ps]; names(d$propM)<-d$spNames
  
  #xx<-data.frame(data$keyCatch)
  #tapply(xx$keyVarObsCatch,list(xx$s),range)
  
  d$keyCatch<- data$keyCatch[data$keyCatch[,'s'] %in% ps,]
  
    k<-sort(unique(d$keyCatch[,'keyVarObsCatch']))
    kk<-1:length(k)
    names(kk)<-k
    d$keyCatch[,'keyVarObsCatch']<-kk[as.character(d$keyCatch[,'keyVarObsCatch'])]
    
    k<-ps
    kk<-1:length(ps)
    names(kk)<-k
    d$keyCatch[,'s']<-kk[as.character(d$keyCatch[,'s'])]
    
  d$logCatchObs<-data$logCatchObs[d$keyCatch[,'obs.no']]
  d$keyCatch[,'obs.no']<-1:dim(d$keyCatch)[[1]]
  
 
  # survey
  d$keySurvey.overview<-data$keySurvey.overview[data$keySurvey.overview[,'s'] %in% ps,,drop=FALSE]
  #str( d$keySurvey.overview)
  d$fleetNames<-rownames(d$keySurvey.overview)
  foundTC<-d$keySurvey.overview[,'techCreep']>0
  d$keySurvey.overview[foundTC,'techCreep']<-1:sum(foundTC)
  if (sum(data$keySurvey.overview[,'s'] %in% ps) ==1) {
    d$keySurvey.overview<-matrix(d$keySurvey.overview,nrow=1)
    colnames(d$keySurvey.overview)<- colnames(data$keySurvey.overview)
  }
   
  d$keySurvey<-data$keySurvey[data$keySurvey[,'s'] %in% ps,,drop=FALSE]
  k<-sort(unique(d$keySurvey[,'keyVarObsSurvey']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyVarObsSurvey']<-kk[as.character(d$keySurvey[,'keyVarObsSurvey'])]
  
  k<-sort(unique(d$keySurvey[,'keyCatchability']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyCatchability']<-kk[as.character(d$keySurvey[,'keyCatchability'])]
  d$logSurveyObs<-data$logSurveyObs[d$keySurvey[,'obs.no']]             
 
  k<-d$keySurvey[,'obs.no']
  d$keySurvey[,'obs.no']<-1:length(d$logSurveyObs)
  kk<-d$keySurvey[,'obs.no']
  names(kk)<-k
  d$keySurvey.overview[,'first']<-kk[as.character(d$keySurvey.overview[,'first'])]
  d$keySurvey.overview[,'last']<-kk[as.character(d$keySurvey.overview[,'last'])]
  
   
  k<-ps
  kk<-1:length(ps)
  names(kk)<-k
  d$keySurvey[,'s']<-kk[as.character(d$keySurvey[,'s'])]
  d$keySurvey.overview[,'s']<-kk[as.character(d$keySurvey.overview[,'s'])]
  
  k<-sort(unique(d$keySurvey[,'f']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'f']<-kk[as.character(d$keySurvey[,'f'])]
  d$keySurvey.overview[,'f']<-kk[as.character(d$keySurvey.overview[,'f'])]
  
  d$recruitYears<-data$recruitYears[ps,,drop=FALSE]
 
  if (smsConf==0) {
    d$nOthSpecies<-0L
    d$consum<-NULL
    d$meanL<-NULL
    d$logPropM2<-NULL
    d$natMor1<-NULL
    d$otherN<-NULL
    d$stom<-NULL
    d$otherFood<-NULL
    d$otherFoodn<-NULL
    #d$predPreySize<-NULL
    d$vulneraIdx<-NULL
    d$suitIdx<-NULL
    d$partM2Idx<-NULL
    d$ovelapIdx<-NULL
    d$overlap<-NULL
    d$stomObsVarIdx<-NULL
  }
  
  
  ### multi species
  if (smsConf>0) {
    d$consum<-data$consum[psAll]; names(d$consum)<-d$allSpNames
    d$meanL<-data$meanL[psAll]; names(d$meanL)<-d$allSpNames
    d$logPropM2<-data$logPropM2[ps]; names(d$logPropM2)<-d$spNames
    d$natMor1<-data$natMor1[ps]; names(d$natMor1)<-d$spNames
    d$otherN<-data$otherN[pso]; names(d$otherN)<-d$othspName
    #d$predPreySize<-data$predPreySize[psAll]; names(d$predPreySize)<-1:npsAll
    d$otherFood<-data$otherFood[psAll]
  
    d$preys<-data$preys[data$preys %in% psAll]
    d$predNames<-rownames(d$info[d$info[,'predator']>0,,drop=FALSE])
    d$preyNames<-rownames(d$info[d$info[,'prey']>0,,drop=FALSE])
    d$othspNames<-rownames(d$info[d$info[,'predator']==2,,drop=FALSE])
    d$nOthSpecies<-length(d$othspNames)
    
    d$usestomObsVar<- !d$info[d$preds, 'stomachVariance'] %in% c(4)
    
    
    stom<-data$stom
   # stom[1,'data'][[1]][[1]]$data[[3]]
    
    stom<-unnest(stom,cols = c(data))
    stom<-filter(stom,pred %in% psAll)
    stom<-unnest(stom,cols = c(data))
    #  filter(stom,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySizeClass )

    oth<-!(stom$prey %in%  d$preys)
    stom[oth,'prey']<-data$nSpecies+1
    stom[oth,'preySizeClass']<-1
    stom[oth,'preySizeW']<-1
    stom[oth,'logPPsize']<-NA
    
    # filter(stom,y==8 & q==1 & pred==1 & predSizeClass==6)
    stom<- stom %>% group_by( area, y, q, yqALK,pred, predSizeClass, predSizeW, stomObsVarIdx, noSampl, phi,prey, preySizeClass, preySizeW ,logPPsize, type ) %>%
       summarize(stomcon=sum(stomcon),.groups = "drop") #%>% ungroup()
    # filter(stom,y==8 & q==1 & pred==1 & predSizeClass==6)
    oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName)) %>% rename(newNo=s)
    stomVarNo<- d$info[d$info[,'predator']>0,'s']
    stom<-left_join(stom,oldNew,join_by(pred==oldNo)) %>% mutate(pred=as.integer(newNo),newNo=NULL,newName=NULL,stomObsVarIdx=match(pred,stomVarNo))
    oldNew<-rbind(filter(oldNew,newName %in% d$preyNames),data.frame(newName=d$otherFoodName,newNo=d$nSpecies+1,oldNo=data$nSpecies+1))
    stom<-left_join(stom,oldNew,join_by(prey==oldNo)) %>% mutate(prey=as.integer(newNo),newNo=NULL,newName=NULL)
    #filter(stom,y==8 & q==1 & pred==1 & predSizeClass==6)    
    
    #parameters
    logStomObsVar<-rep(log(0.5^2),length(d$predNames))
    #names(logStomObsVar)<-d$predNames 
    p$logStomObsVar<-logStomObsVar
    
    
    # vulnerability
    pp<-matrix(0L,nrow=d$nSpecies+1,ncol=d$nSpecies+d$nOthSpecies,dimnames=list(c(d$spNames,d$otherFoodName),c(d$spNames,d$othspNames)))
    pp2<- stom%>% select(pred,prey) %>% unique() %>% arrange(pred,prey)
     ## xtabs(~prey+pred,data=pp2)
    for (i in (1:dim(pp2)[[1]]))  pp[unlist(pp2[i,'prey']),unlist(pp2[i,'pred'])]<-1L
    predPrey<-pp
    pp2$vulneraIdx<-1:dim(pp2)[[1]]
    # print(pp2,n=40)
    
    vulneraIdx<-predPrey
    vulneraIdx[vulneraIdx>0]<-1:sum(vulneraIdx>0)
    d$vulneraIdx<-vulneraIdx # data$vulneraIdx
    p$vulnera<-rep(0,max(vulneraIdx)) #log vulnerability coefficient, parameter 
    
    stom<-left_join(stom,pp2,by = join_by(pred, prey))
    #filter(stom,y==8 & q==1 & pred==1 & predSizeClass==6)
     
    stom<-stom %>% group_by(area,y,q,yqALK, pred,  predSizeClass,prey) %>% 
      mutate(stomconTot=sum(stomcon),minL=min(preySizeClass),nPreyGroups=dplyr::n(),firstL=if_else(minL==preySizeClass,TRUE,FALSE)) %>%
      ungroup() %>% arrange(area,y,q, pred,  predSizeClass,prey,preySizeClass) %>%
      group_by(area,y,q,yqALK, pred,  predSizeClass) %>% mutate(nPreyGroups=sum(firstL),preyIdx=cumsum(firstL), minL=NULL) %>% ungroup()
    #filter(stom,y==8 & q==1 & pred==1 & predSizeClass==6)
    
      d$stom<-stom %>% select(area, y,q,yqALK,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
       nest(data=c(pred,predSizeClass,predSizeW,stomObsVarIdx,noSampl,   phi,nPreyGroups,prey,preySizeClass,preySizeW,logPPsize,type,stomcon, stomconTot ,preyIdx,firstL,vulneraIdx)) %>% 
       rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preySizeClass,preySizeW,logPPsize,type,firstL,preyIdx,stomconTot ,stomcon,vulneraIdx))))
    
     # d$stom[1,'data'][[1]][[1]]$data[[3]]
      
    otherFoodNo<-data$nSpecies+1   
    
    d$overlapIdx<-filter(data$overlapIdx,predNo %in% psAll) %>% filter(preyNo %in% c(d$preys, otherFoodNo))
    if (dim(d$overlapIdx)[[1]] >0) {
      oo<-sort(unique(d$overlapIdx$overlap))
      oo<-data.frame(overlap=oo,newOverlap=1:length(oo))
      d$overlapIdx<-left_join(d$overlapIdx,oo,by = join_by(overlap)) %>% mutate(overlap=newOverlap,newOverlap=NULL)
      
      #oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName))
      oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName)) %>% rename(newNo=s)
      d$overlapIdx<-left_join(d$overlapIdx,oldNew,join_by(predNo==oldNo)) %>% mutate(predNo=as.integer(newNo),newNo=NULL,newName=NULL)
      oldNew<-filter(oldNew,oldNo !=otherFoodNo)
    
      oldNew<-rbind(oldNew,data.frame(newName=data$otherFoodName,oldNo=otherFoodNo,newNo=d$nSpecies+1 ))
      d$overlapIdx<-left_join(d$overlapIdx,oldNew,join_by(preyNo==oldNo)) %>% mutate(preyNo=as.integer(newNo),newNo=NULL,newName=NULL)
      #d$overlapIdx; data$overlapIdx
      p$overlapP<-rep(0,length(unique(d$overlapIdx$overlap))) # log overlap parameter
    } else p$overlapP<- numeric()
    
   
    d$overlap<-data$overlap[,psAll,c(d$spNames,d$otherFoodName)]
    #str(data$overlap);str(d$overlap)
    
    alk<-data$alk
    alk<-unnest(alk,cols = c(data))
    alk<-filter(alk,s %in% psAll)
    alk<-unnest(alk,cols = c(data)) %>% ungroup()
    #filter(alk,y==8 & q==1 & s==1 & sizeClass==6)
    #filter(alk,y==8 & q==1 & s==2 & sizeClass==6)
    #oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName))
    oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName)) %>% rename(newNo=s)
    alk<-left_join(alk,oldNew,join_by(s==oldNo)) %>% mutate(s=as.integer(newNo),newNo=NULL,newName=NULL)
    
    fl<- alk %>% group_by(area,y,q,s,a) %>% summarize(minSize=min(sizeClass), maxSize=max(sizeClass)  , .groups = "drop") #%>% ungroup()
    d$maxSizeCl<-max(fl$maxSize)
    d$alk<- left_join(alk %>% select(area,y,q,s,a,sizeClass,meanLength,alk),fl,by = join_by(area, y, q, s, a)) %>% 
      nest(data=c(a,s,minSize,maxSize,sizeClass,meanLength,alk)) %>% 
      rowwise() %>% mutate(data=list(data %>% nest(data=c(sizeClass,meanLength,alk))))
    
    

    suitIdx<-unnest(d$suitIdx,cols = c(data))
    suitIdx<-filter(suitIdx,predNo %in% psAll) %>% unnest(cols = c(data)) %>% filter(preyNo %in% c(psAll,data$otherFoodn))
    # suitIdx
    #oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName))
    oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName)) %>% rename(newNo=s)
    suitIdx<-left_join( suitIdx,oldNew,join_by(predNo==oldNo)) %>% mutate(predNo=as.integer(newNo),newNo=NULL,newName=NULL)
    oldNew<-rbind(filter(oldNew,newName %in% d$preyNames),data.frame(newName=d$otherFoodName,newNo=d$nSpecies+1,oldNo=data$nSpecies+1))
    suitIdx<-left_join(suitIdx,oldNew,join_by(preyNo==oldNo)) %>% mutate(preyNo=as.integer(newNo),newNo=NULL,newName=NULL,vulneraIdx =NULL)
   
    vul<-stom %>% select(pred,prey,vulneraIdx) %>% unique()
    suitIdx<-left_join(suitIdx,vul,join_by(predNo==pred,preyNo==prey)) 
    
    d$suitIdx<-suitIdx %>% nest(data=c(predNo,predAge,predW,preyNo,preyAge,preyW,logRatio, vulneraIdx)) %>% 
      rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,preyW,logRatio,vulneraIdx))))
    
    
    
    partM2Idx<-unnest(d$partM2Idx,cols = c(data))
    partM2Idx<-filter(partM2Idx,predNo %in% psAll) %>% unnest(cols = c(data)) %>% filter(preyNo %in% c(psAll,data$otherFoodn))
    # partM2Idx
    #oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName))
    oldNew<-left_join(sNew, sOld,by=join_by(newName==oldName)) %>% rename(newNo=s)
    partM2Idx<-left_join( partM2Idx,oldNew,join_by(predNo==oldNo)) %>% mutate(predNo=as.integer(newNo),newNo=NULL,newName=NULL)
    oldNew<-rbind(filter(oldNew,newName %in% d$preyNames),data.frame(newName=d$otherFoodName,newNo=d$nSpecies+1,oldNo=data$nSpecies+1))
    partM2Idx<-left_join(partM2Idx,oldNew,join_by(preyNo==oldNo)) %>% mutate(preyNo=as.integer(newNo),newNo=NULL,newName=NULL,vulneraIdx =NULL)
    
    d$partM2Idx<-partM2Idx %>%
      nest(data=c(predNo,predAge,preyNo,preyAge,partM2,suit,availFood)) %>% 
      rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,partM2,suit,availFood))))
    
    
    d$otherFoodn<-oldNew[oldNew$newName==d$otherFoodName,'newNo']
  }
  
  #parameters
  if (smsConf==0) {  # single species
    p$vulnera<-numeric(0)
    p$overlapP<-numeric(0)
    p$logStomObsVar<-numeric(0)
  }
   

  p$logSdLogObsCatch<-parameters$logSdLogObsCatch[1:max(d$keyCatch[,'keyVarObsCatch'] )]
  
  p$logCatchability<-parameters$logCatchability[1:max(d$keyCatchability)] 
  p$logSdLogObsSurvey<-parameters$logSdLogObsSurvey[1:max(d$keyVarObsSurvey)]
  
  if (max(d$keyLogFstaSd) >0) p$logSdLogFsta<-parameters$logSdLogFsta[1:max(d$keyLogFstaSd)] else p$logSdLogFsta<-numeric(0) 
  if (max(d$keyVarLogN)<0)  p$logSdLogN<-rep(0,0) else p$logSdLogN<-parameters$logSdLogN[1:max(d$keyVarLogN)]  
  
  x<-do.call(c,lapply(d$keylogSeparAgeF,function(x) unique(as.vector(x)))); x<-x[x>0]
  if (length(x)>0) p$logSeparAgeF<-parameters$logSeparAgeF[1:length(x)] else  p$logSeparAgeF<-numeric(0)
  
  
  if (any(d$spNames %in% rownames(parameters$logYearEffectF))) {
   x<-intersect(d$spNames,rownames(parameters$logYearEffectF) )
   p$logYearEffectF<-parameters$logYearEffectF[x,,drop=FALSE]
  }
  else p$logYearEffectF<-numeric(0)
 
  if (sum(d$inclSsbR)==0) p$logSsbRsd<-numeric(0) else p$logSsbRsd<-rep(0.3,sum(d$inclSsbR>0))
  p$rho<-parameters$rho[ps] 
  p$rec_loga<-parameters$rec_loga[ps]         
  p$rec_logb<-parameters$rec_logb[ps] 
  if (sum(foundTC)>0) p$logTechCreep<-parameters$logTechCreep[1:sum(foundTC)] else p$logTechCreep<-numeric(0)
  if (sum(d$nlogN) >0) p$Un<-parameters$Un[1:sum(d$nlogN),,drop=FALSE]   else p$Un<-numeric(0)
  if (sum(d$nlogF) >0) p$Uf<-parameters$Uf[1:sum(d$nlogF),,drop=FALSE]  else p$Uf<-numeric(0)
  
 
  x<-numeric(0)
  for (s in ps) if (data$logNfirstYparamfromTo[s,1] >0) x<-c(x,p$logNfirstYparam[data$logNfirstYparamfromTo[s,1]:data$logNfirstYparamfromTo[s,2]])
  p$logNfirstYparam<-x
  x<-numeric(0)
  for (s in ps) if (data$logNrecruitParamfromTo[s,1] >0) x<-c(x,p$logNrecruitParam[data$logNrecruitParamfromTo[s,1]:data$logNrecruitParamfromTo[s,2]])
  p$logNrecruitParam<-x
  
  
  list(data=d,parameters=p)
}
