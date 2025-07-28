

parToShow<- c("logCatchability","logSdLogFsta","logSdLogN","logSdLogObsCatch","logSdLogObsSurvey",  "logSeparAgeF" ,    
              "logFSeasonal","logYearEffectF","logNrecruitParam","overlapP", "rec_loga","rec_logb", "rho","stomObsVar","vulnera") [9]

inpRdata=list('RSMS_SMS_like_raw','RSMS_SMS_like')
labels=c('RSMS_SMS_like_raw','RSMS_SMS_like')


showSpecies<-as.integer(1:12)
addConf=TRUE

outFormat<-c('screen','pdf','png') [1]
longSpNames<-TRUE
nRowPerPlot=3
isBetterRound<-3


x<-do.call(rbind,lapply(1:length(labels),function(i){
  load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
  if (longSpNames) spNames<<-sms$data$allSpNamesLong else spNames<<-sms$data$allSpNames
  
  extractParameters(sdrep=sms$sdrep,myMap=sms$map,data=sms$data,parameters=sms$parameters,lu=sms$lu,showNotUsed=TRUE)[[2]] %>% 
    filter(name %in% parToShow  & s %in% showSpecies ) %>%  mutate(run=labels[i]) 
})) 


cleanup()

for (ToShow in parToShow) {
  cat(ToShow,'\n')
  # test ToShow<-"logYearEffectF"
  switch(ToShow,
    "logSdLogObsSurvey" = {  
        xx<- x %>% mutate(fleet= unlist(lapply(strsplit(Var1,' '),function(x) paste(x[2:length(x)],collapse=' ')))) %>%
               transmute(run=factor(run,labels),s,var1=SpNames[s],Var2=age,Var3=fleet,estimate,estimate.sd)
        xLab<-'Age'; yLab<-'log(survey obs. sd) +-2*sd'
           },
    "logSdLogN" =  {  
      
      xx<- x  %>%  transmute(run,s,var1=sms$data$allSpNames[s],Var2=age,Var3=Var1,estimate,estimate.sd)
      xLab<-'Age'; yLab<-'Process N, log(sd)'
         },
    "logYearEffectF" ={
      xx<-filter(x,name=="logYearEffectF") %>% transmute(run,s, Var1=factor(s,labels=spNames),Var2=min_age ,
                    upper=exp(estimate+2*estimate.sd),lower=exp(estimate-2*estimate.sd),estimate=exp(estimate),pl=floor((s-1)/nRowPerPlot)) 
      xLab<-'Year'; yLab<-'Year effect in F model'
    },
    "logNrecruitParam"={
      xx<-filter(x,name=="logNrecruitParam") %>% transmute(run,s, Var1=factor(s,labels=spNames),Var2=min_age,
                                                         upper=exp(estimate+2*estimate.sd),lower=exp(estimate-2*estimate.sd),estimate=exp(estimate),pl=floor((s-1)/nRowPerPlot)) 
      xLab<-'Year'; yLab<-'Recruitment' 
    },
    stop(paste("Parameter", parToShow,"is not found:"))
  )

  if (ToShow %in% c("logYearEffectF","logNrecruitParam")){
    p<-by(xx,list(xx$pl), function(xx){
    s<-unlist(xx[1,'s'])
    p<- ggplot(xx, aes(x=Var2, y=estimate, group=run, color=run)) + 
    geom_line(lwd=1) +
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5,lwd=0.7,
                  position=position_dodge(0.40))+
    my_theme()+  
    facet_wrap(nrow=4,Var1 ~ ., scales="free_y")+ ylim(0,NA)+
    labs(x=xLab, y=yLab)
    
    if (length(labels) >1) p<-p+theme(legend.title=element_blank())
    
    print(p)
   })
  }  
}






x2<-do.call(rbind,lapply(1:length(labels),function(i){
  load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
  extractParameters(sdrep=sms$sdrep,myMap=sms$map,data=sms$data,parameters=sms$parameters,lu=sms$lu,showNotUsed=TRUE)[[2]] %>% 
  filter(name==parToShow & s %in% showSpecies ) %>%
    mutate(run=labels[i],ages=paste0('Ages:',min_age,'-',max_age), fl=paste(Var1,ages)) 
})) 


a<-xtabs(exp(estimate)~fl+run,data=x2)
a
if (length(labels) >1) { 
  lab1<-labels[1]; lab2<-labels[2]

  a<-cbind(a,diffrence=a[,lab1]-a[,lab2])  
  a<-cbind(a,ratio=a[,lab1 ]/a[,lab2]) 
  b<-round(a,isBetterRound);b
  bb<-rbind(theSame=sum(b[,'diffrence']==0),firstIsBetter=sum(b[,'diffrence']<0),secondIsBetter=sum(b[,'diffrence']>0))
  rownames(bb)<- c('the same',paste(lab1,'is better'),paste(lab2,'is better'));
  bb
  hist(b[,3],main=paste("difference (",lab1,'-',lab2,')'))
        
  bb<-rbind(theSame=sum(b[,'ratio']==1),firstIsBetter=sum(b[,'ratio']<1),secondIsBetter=sum(b[,'ratio']>1))
  rownames(bb)<- c('the same',paste(lab1,'is better'),paste(lab2,'is better'));
  bb
  hist(b[,'ratio'],main=paste("ratio (",lab1,':',lab2,')'))
}

