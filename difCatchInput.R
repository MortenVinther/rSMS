load(paste0(runName,'.Rdata'),verbose=TRUE)
str(sms,1)

a<-sms$rep$residCatch
head(a)
b<-read.table(file='catch_survey_residuals.out',header=TRUE) %>% filter(data=='catch' & observed> -99) %>%
   mutate(y=Year+sms$data$off.year,a=Age+sms$data$off.age,s=Species.n-15,q=if_else(Quarter==9,1,Quarter)) %>%
  transmute(s,y,q,a,SMSobs=observed,SMSmodel=model,SMSs2=s2)
head(b)

ab<-full_join(a,b,by = join_by(s, y, q, a)) %>% mutate(logSMSobs=log(SMSobs))
head(ab)

a<-ab %>% mutate(difC=abs(exp(logObs)-exp(logSMSobs))) %>% arrange(desc(difC))
head(a)
tail(a)
misC<-filter(a,is.na(difC))
ftable(xtabs(~s+q+a,data=misC))


##### survey


a<-sms$rep$residSurv

a<-left_join(a,a %>% group_by(s) %>% summarize(minf=min(f),.groups='drop')) %>% mutate(f=f-minf+1,minf=NULL,RSMSobs=exp(logObs))  %>% arrange(s,f,y,q,a)
a
b<-read.table(file='catch_survey_residuals.out',header=TRUE) %>% filter(data=='survey' & observed> -99) %>%
  mutate(y=Year+sms$data$off.year,a=Age+sms$data$off.age,s=Species.n-15,q=if_else(Quarter==9,1,Quarter)) %>%
  transmute(s,f=fleet,y,q,a,SMSobs=observed,SMSs2=s2) %>% as_tibble() %>% arrange(s,f,y,q,a)
b


filter(a,s==2)
filter(b,s==2) %>% mutate(logObs=log(SMSobs))

ab<-full_join(a,b,by = join_by(s, f,y, q, a)) %>% mutate(logSMSobs=log(SMSobs))
head(ab)

a<-ab %>% mutate(difC=abs(exp(logObs)-exp(logSMSobs))) %>% arrange(desc(difC))

round(ftable(xtabs(difC~s+f+a,data=a)),0)
filter(a,s==11 & f==4 & a==8)
head(a)
tail(a)
misC<-filter(a,is.na(difC))
ftable(xtabs(~s+q+a,data=misC))
ftable(xtabs(~f+q+a,data=misC))
