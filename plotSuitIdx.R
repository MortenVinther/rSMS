dir=data.path
load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE); 
selPred=c(2L)


# read rsms data
suitIdx<-unnest(sms$data$suitIdx,cols = c(data))
suitIdx
suitIdx$data[[1]]
s<-filter(suitIdx,predNo %in% selPred) %>% unnest(cols = c(data)) %>% filter(preyNo !=sms$data$otherFoodn)  %>% 
  mutate(prey=factor(sms$data$spNames[preyNo]))
s 
summary(s)
by(s,s$predNo,function(x)ftable(xtabs(~prey+preyAge+q,data=x)))


#### compare data for calc of M2 from old SMS and rsms
spl<-sms$data$allSpNamesLong
sp<-1L:length(spl)
names(sp)<-spl

sps<-sms$data$allSpNames
sps<-1L:length(sps)
names(sps)<-sms$data$allSpNames

old<-read_csv(file=file.path("../",compareSMS,"who_eats_whom_level1.csv"))
head(old)
old<-old %>% transmute(y=Year+sms$data$off.year,Year,q=Quarter,Predator,predNo=sps[Predator],predAge=Predator.age+sms$data$off.age,
                       Prey,preyNo=sps[Prey],preyAge=Prey.age+sms$data$off.age,M2SMS=Part.M2) %>% filter(predNo %in% selPred)

dim(old); dim(s)
ab<-full_join(s,old,by = join_by(y, q, predNo, predAge,preyNo, preyAge)) %>% mutate(preyW=NULL,predW=NULL,vulneraIdx=NULL) %>% 
         arrange(y,q,predNo,predAge,preyNo,preyAge)
ab
dim(ab)

ab %>% filter(!is.na(prey) & !is.na(Prey))  # both
ab %>% filter(is.na(prey) & !is.na(Prey))  #old SMS only
ab %>% filter(!is.na(prey) & is.na(Prey))  #rsms  only
###

oldPreds<-c(tail(sms$data$predNames,-5),head(sms$data$predNames,5))
minpp<-scan(file=file.path(data.path,"min_max_size_pref.out"),comment.char="#",n=sms$data$nSpecies*length(sms$data$preds))
minpp<-matrix(minpp,nrow=length(sms$data$preds),ncol=sms$data$nSpecies,byrow = TRUE,
       dimnames=list(oldPreds,sms$data$spNames))
minpp<-array2DF(minpp) %>% transmute(pred=Var1,prey=Var2,minpp=Value)

maxpp<-scan(file=file.path(data.path,"min_max_size_pref.out"),comment.char="#",skip=24)
maxpp<-maxpp[1:(sms$data$nSpecies*length(sms$data$preds))]
maxpp<-matrix(maxpp,nrow=length(sms$data$preds),ncol=sms$data$nSpecies,byrow = TRUE,
              dimnames=list(oldPreds,sms$data$spNames))
maxpp<-array2DF(maxpp) %>% transmute(pred=Var1,prey=Var2,maxpp=Value)
pp<-left_join(maxpp,minpp,by = join_by(pred, prey)) %>% mutate(predNo=sps[pred],preyNo=sps[prey],pred=NULL,prey=NULL)
pp<-filter(pp,maxpp>0)

ab2<-left_join(ab,pp,by = join_by(predNo, preyNo)) %>% mutate(minpp=log(minpp),maxpp=log(maxpp),within=logRatio>=minpp & logRatio<=maxpp,
                                                        tooBig=logRatio<minpp,tooSmall=logRatio>maxpp) 
ab2

filter(ab2,is.na(within)) # no rsms range (min max data)
filter(ab2,within)        # within rsms range
filter(ab2,!within)       # outside rsms range
filter(ab2,tooBig)
filter(ab2,tooSmall)
ab2
by(ab2,ab2$predNo,function(x)ftable(xtabs(tooBig~preyNo+preyAge+predAge,data=filter(x,tooBig & q==3))))
# no cod obs by(ab2,ab2$predNo,function(x)ftable(xtabs(tooSmall~preyNo+preyAge+predAge,data=filter(x,tooSmall & q==3))))

####################################
###  Compare input to SMS and RSMS

s1<-read_table("summary_stom.out") %>% filter(Prey.no>0) %>%
  transmute(Year,y=Year+sms$data$off.year,q=Quarter.no,predator=oldPreds[Predator.no],predNo=sps[predator],predSizeClass=Predator.length.class,Predator.weight,
            preyNo=Prey.no-15,preySizeClass=Prey.length.class,Prey.weight,stomconSMS=stom.input)  %>% arrange(y,q,predNo,predSizeClass,preyNo,preySizeClass) 

stom<-unnest(sms$data$stom,cols = c(data)) %>% unnest(cols = c(data)) %>% filter(prey !=sms$data$otherFoodn)  %>% 
  transmute(area,y,q,preyNo=prey,prey=sms$data$spNames[preyNo],predNo=pred,pred=sms$data$predNames[predNo],predSizeClass,predW=predSizeW,preySizeClass,preyW=preySizeW,logRatio=logPPsize,stomcon) %>%
  arrange(y,q,predNo,predSizeClass,preyNo,preySizeClass) 

stom 
s1
sort(unique(paste(formatC(s1$predNo,w=2),s1$predator,sep='_')))
sort(unique(paste(formatC(stom$predNo,w=2),stom$pred,sep='_')))  

xtabs(~y+q,data=filter(s1,predator=='MAC'))
xtabs(~y+q,data=filter(stom,pred=='MAC'))

xtabs(~y+q,data=filter(s1,predator=='GSE'))
xtabs(~y+q,data=filter(stom,pred=='GSE'))

summary(stom)

b<-full_join(stom,s1,join_by(y, q, preyNo, predNo, predSizeClass, preySizeClass)) %>% 
  select(Year,y,q ,pred,  predator, predNo, predSizeClass, predW,Predator.weight,preyNo, prey,  preySizeClass,   preyW ,Prey.weight,logRatio , stomcon , stomconSMS) %>% 
  arrange(y,q,predNo,predSizeClass,preyNo,preySizeClass) 
b

summary(b)  # some RSMS stomcon are missing
print(filter(b,is.na(stomcon)),n=200)

stomDiff<-b %>% mutate(relDiff= abs((stomcon-stomconSMS)/stomcon),absDiff=abs(stomcon-stomconSMS)) %>% arrange(desc(absDiff))
print(stomDiff,n=300)

###  end Compare input to SMS and RSMS
stom<-unnest(data$stom,cols = c(data)) 
stom
stom<-filter(stom,pred %in% selPred) %>% unnest(cols = c(data)) %>% filter(prey !=data$otherFoodn)  %>% 
  transmute(area,y,q,preyNo=prey,prey=data$spNames[preyNo],predNo=pred,pred=data$predNames[predNo],predW=predSizeW,logRatio=logPPsize)
stom 
summary(stom)
sort(names(stom));sort(names(s)); 

ss<-rbind(stom %>% select(area,y,q,preyNo,predNo,predW,logRatio) %>% mutate(type='stom'),
            s %>% select(area,y,q,preyNo,predNo,predW,logRatio) %>% mutate(type='M2'))  %>%
     mutate(prey=data$spNames[preyNo],pred=data$predNames[predNo])
ss

# UNDGÅ DETTTE LIB library(plyr)
sss<-filter(ss,predNo==1)
by(sss,sss$prey,function(x) by(x,x$type,function(xx) range(xx$logRatio)))
by(sss,sss$prey,function(x) by(x,x$type,function(xx) exp(range(xx$logRatio))))

mu <- ddply(sss, "type", summarise, grp.mean=mean(logRatio))
mu
ggplot(sss, aes(x=logRatio, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)+
  facet_wrap(~prey,  ncol=2, strip.position = 'right')
ggplot(sss, aes(x=logRatio, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),
             linetype="dashed")+
  facet_wrap(~paste(prey,type),  ncol=2, strip.position = 'right',scales='free_y')



out<-by(s,list(s$predNo),function(x) {
  p=unlist(x[1,'predNo'])
  st<-filter(stom,predNo==p)
  tit<- data$predNames[(x$predNo)[1]] 
  a<-ggplot(x,aes(x=predW, y = logRatio, col=prey,shape=prey)) +
    facet_wrap(~paste0('Q:',q),  ncol=2, strip.position = 'right')+
    geom_point( )+
    geom_point(data=st,aes(x=log(predW), y = logRatio), col='black') +
    labs(x='log (Predator size)',y='log(predator/prey size ratio)',title=tit)
    #theme_minimal() +
    #theme( panel.grid.major = element_line(linetype = "blank"),
    #       panel.grid.minor = element_line(linetype = "blank")
    #)
    print(a)
})

cleanup()
out<-by(s,list(s$predNo,s$preyNo),function(x) {
  st<-filter(stom,predNo==unlist(x[1,'predNo']) & preyNo==unlist(x[1,'preyNo']))
  tit<- paste(data$predNames[(x$predNo)[1]], 'eating',data$spNames[(x$preyNo)[1]]) 
  a<-ggplot(x,aes(x=predW, y = logRatio, col=prey,shape=prey)) +
    facet_wrap(~paste0('Q:',q),  ncol=2, strip.position = 'right')+
    geom_point( )+
    geom_point(data=st,aes(x=log(predW), y = logRatio), col='black') +
    labs(x='log (Predator size)',y='log(predator/prey size ratio)',title=tit)
  #theme_minimal() +
  #theme( panel.grid.major = element_line(linetype = "blank"),
  #       panel.grid.minor = element_line(linetype = "blank")
  #)
  print(a)
})

cleanup()
