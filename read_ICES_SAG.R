library(icesSAG)


cod<-getSAG(stock ='cod.27.47d20', year = 2022)
cod<-subset(cod,select=c(Year,recruitment,SSB,catches,F,fishstock))
str(cod)


b<-getSAG(stock ='pok.27.3a46', year = 2023)
str(b)

stocks<-data.frame(stocks=
  c('whg.27.47d',
  'had.27.46a20',
  'pok.27.3a46',
  'mac.27.nea',
  'her.27.3a47d',
  'nop.27.3a4',
  'spr.27.3a4',
  'ple.27.420',
  'sol.27.4'),
  Species.n=c(17,18,19,20,21,24,25,26,27))


b<-lapply(stocks$stocks,function(x) {
  print(x)
  b<-getSAG(stock =x, year = 2023)
  b<-subset(b,select=c(Year,recruitment,SSB,catches,F,fishstock))
  return(b) }
)

stocks<- rbind(data.frame(stocks='cod.27.47d20',Species.n=16), stocks)
a<-do.call(rbind,b)
a<-rbind(a,cod)

b<-merge(x=a,y=stocks,by.x=c('fishstock'),by.y=c('stocks'))
head(b)
bb<-data.frame(Species.n=b$Species.n,Year=b$Year,Rec=b$recruitment, SSB=b$SSB,TSB=0,SOP=b$catches,SOP.hat=0,Yield=0,yield.hat=0,mean.F=b$F,Eaten=0)
head(bb)
xtabs(~Species.n,data=bb)
write.table(bb,file=file.path(data.path,"summary_table_raw_ICES.out"),col.names=TRUE,row.names=FALSE)
            
            
