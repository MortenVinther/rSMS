

sdrep=sms$sdrep
my.map=sms$map
data=sms$data
parToShow<- c("logCatchability","logSdLogFsta","logSdLogN","logSdLogObsCatch","logSdLogObsSurvey", "overlapP", "rec_loga","rec_logb", "rho","stomObsVar","vulnera")[3]

inpRdata<-list('Single','Multi')
labels=c('Single sp','Multi sp')
showSpecies<-1:10


outFormat<-c('screen','pdf','png') [1]
longSpNames<-TRUE
isBetterRound<-3


x<-do.call(rbind,lapply(1:length(labels),function(i){
  load(file=file.path(data.path,paste0(inpRdata[2],".Rdata")))
  extractParameters(sdrep=sms$sdrep,my.map=sms$map,data=sms$data)[[1]] %>% filter(name==parToShow & s %in% showSpecies ) %>%
  mutate(run=labels[i]) 
})) 

switch(parToShow,
  "logSdLogObsSurvey" = {  
      x<- x %>% mutate(fleet= unlist(lapply(strsplit(Var1,' '),function(x) paste(x[2:length(x)],collapse=' ')))) %>%
             transmute(run,s,var1=data$allSpNames[s],Var2=age,Var3=fleet,estimate,estimate.sd)
      xLab<-'Age'; yLab<-'log(survey obs. sd) +-2*sd'
         },

  "logSdLogN" =  {  # Dirichlet distribution (with input of maximum 'concentration' parameter)
    x<- x  %>%  transmute(run,s,var1=data$allSpNames[s],Var2=age,Var3=Var1,estimate,estimate.sd)
    xLab<-'Age'; yLab<-'Process N, log(sd)'
       },
       
  stop(paste("Parameter", parToShow,"is not found:"))
)

x<- x %>% mutate(run=factor(run,labels),Var2=factor(Var2))

by(x,list(x$s), function(xx){
  s<-unlist(xx[1,'s'])
  tit<-data$allSpNamesLong[s]
  p<- ggplot(xx, aes(x=Var2, y=estimate, group=run, color=run)) + 
  geom_line(lwd=1) +
  geom_point()+
  geom_errorbar(aes(ymin=estimate-abs(2*estimate.sd), ymax=estimate+abs(2*estimate.sd)), width=0.5,lwd=1,
                position=position_dodge(0.20))+
  facet_wrap(nrow=2,Var3 ~ ., scales="free_y")+
  labs(x=xLab, y=yLab,title=tit)+
  theme_bw()
})



x2<-do.call(rbind,lapply(1:length(labels),function(i){
  load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
  extractParameters(sdrep=sms$sdrep,my.map=sms$map,data=sms$data)[[2]] %>% filter(name==parToShow & s %in% showSpecies ) %>%
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

