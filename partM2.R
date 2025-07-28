dir<-data.path
load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE); 

a<-sms[['rep']]$partM2Idx 
sp<-1L:length(sms$data$allSpNames)

names(sp)<-sms$data$allSpNames

spl<-sms$data$allSpNamesLong
splPrey<-spl
splPrey[sms$data$otherFoodn]<-'Other'

b<- a  %>% unnest(cols = c(data))  %>%  unnest(cols = c(data)) %>% mutate(pred=spl[predNo],prey=splPrey[preyNo])
b

filter(b,preyNo==7 &y==1 & q==1 & preyAge==2)
bb<-filter(b,y %in% c(1,10,49))
round(ftable(xtabs(partM2~prey+preyAge+y,data=bb)),3)

bb<-filter(b,y %in% c(10))
round(ftable(xtabs(partM2~prey+preyAge+q,data=bb)),3)


old<-read_csv(file=file.path("../",compareSMS,"who_eats_whom_level1.csv"))
head(old)
old<-old %>% transmute(y=Year+sms$data$off.year,Year,q=Quarter,Predator,predNo=sp[Predator],predAge=Predator.age+sms$data$off.age,
                  Prey,preyNo=sp[Prey],preyAge=Prey.age+sms$data$off.age,M2SMS=Part.M2)
ab<-full_join(b,old,by = join_by(y, q, predNo, predAge,preyNo, preyAge)) %>% arrange(y,q,predNo,predAge,preyNo,preyAge)
ab
ab %>% filter(!is.na(M2) & !is.na(M2SMS))

# no old SMS but new
ab %>% filter(is.na(M2) & !is.na(M2SMS))

# old but no new
ab %>% filter(!is.na(M2) & is.na(M2SMS))

