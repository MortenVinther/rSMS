load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE); 

a<-sms[['data']]$stom
sp<-1L:length(sms$data$allSpNames)
names(sp)<-sms$data$allSpNames

spl<-sms$data$allSpNamesLong
splPrey<-spl
splPrey[sms$data$otherFoodn]<-'Other'

b<- a  %>% unnest(cols = c(data))  %>%  unnest(cols = c(data)) %>% mutate(predC=spl[pred],preyC=splPrey[prey])


aa<-extractParameters(sms$sdrep,sms$map,sms$data,sms$parameters,lu=sms$lu,showNotUsed=FALSE)[[2]]
aa<-filter(aa,name=="logStomObsVar") %>% transmute(pred=s,StomObsVar=exp(estimate))

bb<-left_join(aa,b,by = join_by(pred)) %>%  
  transmute(year=y-sms$data$off.year,q,pred=paste(formatC(pred,w=2),predC),size=predSizeClass,yq=paste(year,q,sep='-'),phi,alfa0=noSampl*StomObsVar - 1.0) %>% unique()
bb

cat("phi/alfa0\n")
by(bb,bb$pred,function(x) {
 a1<-xtabs(phi~size+yq,data=x) 
 a2<-xtabs(alfa0~size+yq,data=x) 
 round(a1/a2,2)
})


cat("phi\n")
by(bb,bb$pred,function(x) {
  a1<-xtabs(phi~size+yq,data=x) 
  round(a1,2)
})


cat("alfa0\n")
by(bb,bb$pred,function(x) {
  a2<-xtabs(alfa0~size+yq,data=x) 
  round(a2,2)
})


bbb<-bb %>% group_by(pred,size) %>% summarize(minAlfa0=min(alfa0),maxAlfa0=max(alfa0),.groups='drop')
round(xtabs(minAlfa0~pred+size,data=bbb),1)
round(xtabs(maxAlfa0~pred+size,data=bbb),1)
