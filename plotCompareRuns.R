
my_theme <- function() {
  theme_bw()+
    theme( #panel.grid.major = element_line(linetype = "blank"),
      panel.grid.minor = element_line(linetype = "blank"),
      #axis.text.x = element_text(angle = tAngle, vjust = 0.5),
      legend.title = element_blank(),
      legend.position = "none",
    )  
}

plotSSB<-function(x,sp,mult=1E-3,tit='SSB') {
  x<-filter(x,s==sp) %>% mutate(SSB=SSB*mult)
  a<-ggplot(data=x, aes(x=year, y=SSB,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y='(1000 tonnes)',title=tit)+
    my_theme()
  a
}

SSB_R<-function(SSB,model,a,b) {
  switch( as.character(model),
          "0"=  NA,                                                          ## straight RW
          "1"=  exp(a+log(SSB)-exp(b)*SSB),          ## Ricker
          "2"=  exp(a+log(SSB)-log(1+exp(b)*SSB)),   ## B&H
          "3"=  rep(exp(a),length(SSB)),                                                                ## GM
          "4"=  exp(a-b+log(SSB-(0.5*(SSB-exp(b)+abs(SSB-exp(b)))))),  #Hockey stick
          "6"=  exp(a-b+log(SSB-(0.5*(SSB-exp(b)+abs(SSB-exp(b)))))),  #Hockey stick
          stop(paste0("SR model code ",model," not recognized"))          ## error
  )
}

plotSSBR<-function(x,sp,multN=1E-3,multBIO=1E-3,tit='SSB recruitment',titSp='cod',SR) {
  x<-filter(x,s==sp) %>% mutate(use_rec=as.factor(use_rec+1))
  tit<-paste(tit, titSp[sp])
  sr<-filter(SR,s==sp)
  #a<-ggplot(data=x, aes(x=SSB*multBIO, y=recruit*multN,col=use_rec,shape=use_rec,group=use_rec)) +
  a<-ggplot(data=x, aes(x=SSB*multBIO, y=recruit*multN,colour=use_rec)) +
    
    #geom_line()+
    geom_point()+
    geom_line(mapping = aes(x =SSB*multBIO,y =recruit*multN), data=sr,  colour = 'red',size=1)+
    ylim(0,NA)+
    labs(x='SSB (xxx) ',y='Recruits (XXXXX)',title=tit)+
    my_theme()
  a
}

plotSSBRibbon<-function(x,sp,mult=1E-3,tit='SSB') {
  x<-filter(x,s==sp ) %>% mutate(mid=mid*mult,low=low*mult,high=high*mult)
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y='(1000 tonnes)',title=tit)+
    my_theme()
  a
}


plotLabels<-function(tit='Hej',labels){
  nlabels<-length(labels)
  ll<-max(nchar(labels))
  ggplot(data.frame(x=1,y=1,run=factor(labels,labels)),
         aes(x,y,col=run,shape=run))+
    geom_line(linewidth = 1.5)+
    geom_point(aes(shape=run),size=4)+
    xlim(2,3)+ylim(2,3)+
    #labs(color= tit)+
    #guides( size = FALSE,shape=FALSE)+
    labs(title=tit)+
    theme_void()+
    theme(plot.title = element_text(color="black", size=25, face="bold",margin=margin(t=20, r=0,b =0, l=15, unit="pt")),
          legend.title = element_blank(),
          legend.text=element_text(size=17-ll/5),
          legend.position = "inside",
          legend.position.inside=c(0.15,0.4+0.04*nlabels),
          legend.justification = c("left"),
          legend.box.margin = margin(15, 15, 15, 15),
          legend.key.size = unit(3, "lines"),
          legend.key.height = unit(1.3, "lines"),                
          legend.key.width = unit(2, "lines"), 
          legend.box.background = element_rect(color="black", linewidth=1))
}  

plotRecruits<-function(x,sp,mult=1E-6,tit='Recruits') {
  x<-filter(x,s==sp) %>% mutate(recruit=recruit*mult)
  
  a<-ggplot(data=x, aes(x=year, y=recruit,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y='(billions)',title=tit)+
    my_theme() 
  a
}

plotM2<-function(x,sp) {
  xx<-filter(x,s==sp) 
  showAges<-sort(unique(x$age))
  p<-lapply(showAges,function(a){
    tit<-paste('Age:',a)
    xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
    #print(xtabs(M2~run+year,data=xxx))
    pp<-ggplot(data=xxx, aes(x=year, y=M2,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme()
   })
}


plotFage<-function(x,sp,allInOne=FALSE,titSp='cod',tit='F:') {
  xx<-filter(x,s==sp) 
  if (allInOne) {
    tit<-paste(tit, titSp[sp])
    xx$age<-factor(xx$age)
    pp<-ggplot(data=xx, aes(x=year, y=FisQ,col=age,shape=age,group=age)) +
      geom_line()+
      geom_point()+
      ylim(0,NA)+
      labs(x='Year',y=' ',title=tit)+
      my_theme() +
      theme(legend.title = element_text(),    legend.position = "bottom")
    #theme(    legend.position = "bottom")
    p<-list(pp)
  } else {
    showAges<-sort(unique(x$age))
    p<-lapply(showAges,function(a){
      tit<-paste('Age:',a)
      xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
      pp<-ggplot(data=xxx, aes(x=year, y=FisQ,col=run,shape=run,group=run)) +
        geom_line()+
        geom_point()+
        ylim(0,NA)+
        labs(x='Year',y=' ',title=tit)+
        my_theme() 
    })    
  }
  return(p)
}

plotNage<-function(x,sp) {
  xx<-filter(x,s==sp) 
  showAges<-sort(unique(x$age))
  p<-lapply(showAges,function(a){
    tit<-paste('Age:',a)
    xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
    #print(xtabs(M2~run+year,data=xxx))
    pp<-ggplot(data=xxx, aes(x=year, y=N,col=run,shape=run,group=run)) +
      geom_line()+
      geom_point()+
      ylim(0,NA)+
      labs(x='Year',y=' ',title=tit)+
      #  facet_wrap(~age)+
      my_theme()
  })
}


plotRecruitsRibbon<-function(x,sp,mult=1E-6,tit='Recruits') {
  x<-filter(x,s==sp ) %>% mutate(mid=mid*mult,low=low*mult,high=high*mult)
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y='(billions)',title=tit)+
    my_theme()
  a
}

plotF<-function(x,sp,tit='F') {
  x<-filter(x,s==sp) 
  a<-ggplot(data=x, aes(x=year, y=Fbar,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme() 
  a
}


plotFRibbon<-function(x,sp,tit='F') {
  x<-filter(x,s==sp ) 
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme()
  a
}


# 1. Summary plots (Recruit, SSB and F), option "Summary", 
# 2. Summary pots with uncertainties, option "SummaryConf" 
# 3. M2 by age, option "M2" input rep$res

####
plotCompareRunSummary<-function(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat"),showSpecies=1:100,
                        inpRdata=list('Single','Multi'),
                        labels=c('Single sp','Multi sp'),
                        outFormat=c('screen','pdf','png')[1],
                        showAges=0:4,
                        useF=c("simpleSumqF","fromCandZ")[2],
                        multN=1E-6,multBIO=1E-3,
                        ncols=2,
                        allInOne=FALSE,
                        fileLabel='pl',outDir=data.path,
                        longSpNames=TRUE) {
  
  if (Type %in% c("SummaryConf","Summary"))  allInOne<-FALSE
  
  if (Type %in% c("Summary","SSBrec")) {
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")),verbose=F)
      b<-sms$rep$resSummary %>% mutate(year=as.numeric(year),run=labels[i]) %>% filter(s %in% showSpecies)
      if (Type =="SSBrec") b<- b%>% rowwise%>% mutate(use_rec=sms$data$recruitYears[s,year+sms$data$off.year])
      b
    })) 
    
  } else if (Type %in% c("SummaryConf","SSBrecConf")) {
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-data.frame(value=sms$sdrep$value,sd=sms$sdrep$sd,variable=names(sms$sdrep$value)) 
      nvar<-4
      b<-cbind(expand.grid(s=1:sms$data$nSpecies,year=sms$data$years,dummy=1:nvar),out) %>% filter(s %in% showSpecies) %>%
         mutate(run=labels[i],species=sms$data$spNames[s],mid=exp(value),low=exp(value-2*sd),high=exp(value+2*sd))
      if (Type =="SSBrecConf") b<- b%>% rowwise%>% mutate(use_rec=sms$data$recruitYears[s,year+sms$data$off.year])
      b
    }))
    
  } else if (Type=='M2'){
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-sms$rep$res %>% transmute(year,age,s,species,run=labels[i],M2) %>% 
        filter(s %in% showSpecies & age %in% showAges &M2>0) %>% group_by(s,species,year,age,run) %>%
        summarize(M2=sum(M2),.groups='drop') %>% group_by(s,age) %>% mutate(maxM2=max(M2)) %>% filter(maxM2>1E-5) %>%
        mutate(maxM2=NULL)
    }))
    
  } else if (Type=='F' | Type=="ExpPat"){
  x<-do.call(rbind,lapply(1:length(labels),function(i){
    load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
    if (c('F') %in% colnames(sms$rep$res)) sms$rep$res$FisQ<-sms$rep$res$F
    if (c('FF') %in% colnames(sms$rep$res)) sms$rep$res$FisQ<-sms$rep$res$FF
    out<-sms$rep$res %>% transmute(year,age,s,species,run=labels[i],FisQ) %>% 
      filter(s %in% showSpecies & age %in% showAges & FisQ>0) %>% group_by(s,species,year,age,run) %>%
      summarize(FisQ=sum(FisQ),.groups='drop') %>% group_by(s,age) %>% mutate(maxF=max(FisQ)) %>% filter(maxF>1E-5) %>%
      mutate(maxF=NULL)
  }))
    
  } else if (Type=='N'){
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-sms$rep$res %>% transmute(year,quarter,age,s,species,run=labels[i],N=N*multN) %>% 
        filter(s %in% showSpecies & age %in% showAges & N>0 & (quarter==data$fqa[age+data$off.age])) %>% 
        mutate(quarter=NULL)
    }))
  }
  
  x<- x %>% mutate(run=factor(run,labels))

  if (Type %in% c("SSBrecConf","SSBrec")) {
    SR<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      a<-suppressWarnings(extractParameters(sdrep=sms$sdrep,myMap=sms$map,data=sms$data,parameters=sms$parameters,lu=sms$lu)[[2]]) %>% 
          filter(s %in% showSpecies & name %in% c("rec_loga","rec_logb")) %>%
        transmute(s,sp=Var1,name,estimate,run=labels[i]) %>%
        pivot_wider(names_from = name, values_from = estimate)
      if (!any(grepl('rec_logb',colnames(a)))) a$rec_logb<-NA  
      b<-data.frame(sp=names(sms$data$stockRecruitmentModelCode),model=data$stockRecruitmentModelCode,addInf=data$stockRecruitmentModelCodeAdd)
      b<-left_join(a,b,by = join_by(sp)) %>% mutate(model=as.integer(model))
      b<-left_join(b,x%>% group_by(s) %>% summarize(maxSSB=max(SSB),.groups='drop') %>% mutate(minSSB=1),by = join_by(s))
      
      a<-b %>% mutate(step=round((maxSSB-minSSB)/100),rec_logb=if_else(is.na(rec_logb) & addInf>0,log(addInf),rec_logb)) 
      do.call(rbind,lapply(1:nrow(a),function(i){
        aa<-a[i,]
        SSB<-seq(aa$minSSB,aa$maxSSB,aa$step)
        cat(i,'\n')
        data.frame(s=aa$s,run=aa$run,SSB,recruit=SSB_R(SSB=SSB,model=aa$model,a=aa$rec_loga,b=aa$rec_logb)) 
      }))
      
    }))
  }
 
  
  
  load(file=file.path(data.path,paste0(inpRdata[1],".Rdata")))
  allSpNamesLong<-sms$data$allSpNamesLong
  names(allSpNamesLong)<-sms$data$allSpNames
  showSpecies<-showSpecies[showSpecies %in% sort(unique(x$s))]
   
  if (Type=="ExpPat"){
    fala<-as.data.frame(sms$data$options@avg.F.ages) %>% transmute(s=1:n(),fa=`first-age`+ sms$data$off.age, la=`last-age`+ sms$data$off.age)
    ff<-left_join(fala,x,by = join_by(s)) %>% filter(age>=fa &age <=la)%>% 
      group_by(s,species,year,run) %>% summarize(meanF=mean( FisQ),.groups='drop')
    x<-left_join(x,ff,by = join_by(s, species, year, run)) %>%mutate(FisQ=FisQ/meanF,meanF=NULL) 
  }
  
  
  cleanup()
  for (mysp in showSpecies) {
     tit<-allSpNamesLong[mysp]
     titFile<-tit
     if (Type %in% c("M2")) tit<-paste('M2\n',tit)
     if (Type %in% c("F"))  tit<-paste('F\n',tit)
     if (Type %in% c("N"))  tit<-paste('N\n',tit)
     if (!allInOne) myLab<-plotLabels(tit,labels)
   
     if (Type=="Summary") {
       pSSB<-plotSSB(x,sp=mysp,mult=multBIO)
       if (useF=="fromCandZ") x$Fbar<-x$FbarAnn
       pF<-plotF(x,sp=mysp)
       pRecruit<-plotRecruits(x,sp=mysp,mult=multN)
       
     } else if (Type=="SummaryConf"){
        pSSB<-plotSSBRibbon( x=filter(x,variable=='SSB'),sp=mysp,mult=multBIO)
        pRecruit<-plotRecruitsRibbon(x=filter(x,variable=='recruit'),sp=mysp,mult=multN)
        if (useF=="simpleSumqF")  pF<-plotFRibbon(x=filter(x,variable=='FBAR'),sp=mysp)
        if (useF=="fromCandZ")  pF<-plotFRibbon(x=filter(x,variable=="FBARann"),sp=mysp)

        
     } else if (Type=="SSBrec") {
        pSSBR<-plotSSBR(x,sp=mysp,titSp=allSpNamesLong,multN=multN,multBIO=multBIO,SR=SR)
      
     } else if (Type=="SSBrecConf"){
        pSSB<-plotSSBRibbon( x=filter(x,variable=='SSB'),sp=mysp,mult=1E-3)
        pRecruit<-plotRecruitsRibbon(x=filter(x,variable=='recruit'),sp=mysp,mult=1E-6)
        pF<-plotFRibbon(x=filter(x,variable=='FBAR'),sp=mysp)
          
    }else if (Type=="M2"){
       pAge<-plotM2(x,sp=mysp) 
  
    } else if (Type=="F"){
       pAge<-plotFage(x,sp=mysp,allInOne=allInOne,titSp=allSpNamesLong,tit='Fishing mortality: ') 
       
    } else if (Type=="ExpPat"){
      pAge<-plotFage(x,sp=mysp,allInOne=allInOne,titSp=allSpNamesLong,tit='Exploitation pattern: ') 
      
    } else if (Type=="N"){
      pAge<-plotNage(x,sp=mysp) 
   } 
   if (Type %in% c("Summary","SummaryConf")) pl<-suppressWarnings(print(ggarrange(myLab,pSSB,pRecruit,pF, ncol = 2, nrow = 2)))
   if (Type %in% c("SSBrecConf","SSBrec")) pl<-suppressWarnings(print(ggarrange(pSSBR, ncol = 1, nrow = 1)))
     
   if (Type %in% c("M2","F","N","ExpPat")) {
     if (allInOne)  nc<-1 else nc<-ncols
      nr<-(length(showAges)+1) %/% nc + (length(showAges)+1) %% nc
      if (allInOne) nr<-1
      if (allInOne)  pl<-suppressWarnings(print(ggarrange(plotlist=lapply(pAge,print), ncol = nc, nrow = nr)))
      if (!allInOne) pl<-suppressWarnings(print(ggarrange(myLab,plotlist=lapply(pAge,print), ncol = nc, nrow = nr)))
   }
    if (outFormat=='screen'){
     print(pl)
    } else {
       ggexport(pl, filename = file.path(outDir,paste0(fileLabel,'_',Type,'_',titFile,'.',outFormat)),width = 900,height = 600, pointsize = 16)
       cleanup()
    }
  }
}


plotSeasoData<-function(x,tit='F',ncols=2) {
  a<-ggplot(data=x, aes(x=year, y=value,col=quarter,shape=quarter,group=quarter)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme() +  theme(legend.title=element_text( size = rel(1.25))  ,  legend.position = "right")+
    facet_wrap(vars(age),ncol=ncols, scales="free_y")
    
  a
}



plotSeasonalData<-function(inp,Type=c("N","F","C","M","M1",'M2','Z',"WEST","WECA","propMat","seasFprop","FiProp")[2],
                                CombineSeason=FALSE,
                                showSpecies=1:100,
                                outFormat=c('screen','pdf','png')[1],
                                showAges=0:4,
                                multN=0.001,
                                ncols=2,
                                fileLabel='pl',
                                cummulate=FALSE,  
                                longSpNames=TRUE){
 
  a<-c("N","F",   "C",    "M","M1","M2",'Z',"WEST","WECA","propMat","seasFprop","FiProp")
  b<-c("N","FisQ","canum","M","M1","M2",'Z','west',"weca","propMat","seasFprop","FiProp")
  multiplier<-rep(1,length(a)); multiplier[1]<-multN
  labels<-c("Stock numbers","Fishing mortality","Catch numbers at age","M","M1",'M2','Z, total Mortality','Weight at age in the stock',"Weight at age in the catch","propMat","seasFprop","proportion of F")
  names(b)<-a
  names(labels)<-a
  names(multiplier)<-a
  
  load(file=file.path(data.path,paste0(inp,".Rdata")),verbose=T)
  recSeason<-sms[['data']]$recSeason
  spNames<-sms[['data']]$allSpNames
  spNamesLong<-sms[['data']]$allSpNamesLong
  cat(spNamesLong,'\n')
  x<-sms[['rep']]$res %>% filter(s %in% showSpecies & age %in% showAges) 
  x$value<-unlist(x[,b[Type]])
  x<-x  %>% transmute(species,s,year,quarter=factor(quarter),age=paste('Age',age),value=value*multiplier[Type])
  
  
  if (CombineSeason) {
    if (Type=='N') filter(x,(a>1 & Quarter==1) | (a==1 & quarter==recSeason)) 
    ## more to come ....  
  }
  if (cummulate) {
    x<-x %>% arrange(s,year,age,quarter) %>% group_by(s,year,age) %>% mutate(value=cumsum(value)) 
  }
  cleanup()
  p<-by(x,x$s,function(x){
    pl<-plotSeasoData(x,tit=paste0(spNamesLong[unlist(x[1,'s'])],', ',labels[Type]),ncol=ncols) 
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = paste0(fileLabel,'_',Type,'_',titFile,'.',outFormat),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
 
}  
 


