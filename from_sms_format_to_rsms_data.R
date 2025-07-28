From_SMS_format_to_rsms<-function(otherPredExist=TRUE,catchMultiplier=1,dir=data.path,sms.dat="rsms.dat") {
 
# test   otherPredExist=multi;catchMultiplier=1;dir=data.path;sms.dat="rsms.dat";  )
  

sms<-read.RSMS.control(dir,file=sms.dat)  
  
la<-sms@max.age.all
fa<-sms@first.age
rec.season <-sms@rec.season
years<-c(1,1)
years[1]<-sms@first.year.data
years[2]<-sms@last.year.data
ny<-years[2]-years[1]+1
npr<-sum(sms@species.info[,'predator']>=1)
nsp<-sms@no.species
nq<-sms@last.season
noAreas<-1L #sms@no.areas

first.VPA<-which(sms@species.info[,'predator'] %in% c(0,1))[1]

#############  catch data
scanData<-function(name='canum.in',announce=FALSE){
  if (announce) cat('Reading ',name,'\n')
  x<-scan(file.path(dir,name),comment.char='#',quiet=TRUE)
  stopifnot(tail(x,1)== -999)
  x<-head(x,-1)
}
CATCHN<-scanData('canum.in')
WCATCH<-scanData('weca.in')
Prop.landed<-scanData('proportion_landed.in'); Prop.landed<-Prop.landed[1:length(WCATCH)]
Prop.season.F<-scanData('proportion_of_annual_f.in'); Prop.season.F<-Prop.season.F[1:length(WCATCH)]

b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
b<-data.frame(b,CATCHN=CATCHN,WCATCH=WCATCH,PROP_CAT=Prop.landed,PROP_SEASON_F=Prop.season.F)
b<-subset(b,select=c(year,species.n,quarter,sub_area,age,WCATCH,CATCHN,PROP_CAT,PROP_SEASON_F))
b$CATCHN<-catchMultiplier*b$CATCHN


out<-list(catch=b)
############## bio data
WSEA<-scanData('west.in');WSEA<-WSEA[((first.VPA-1)*noAreas*ny*(la-fa+1)*nq+1):length(WSEA)]
PROPMAT<-scanData('propmat.in')
M<-scanData('natmor.in')
#M<-scanData('natmorm1m2.out')
M1<-scanData('natmor1.in')
PROP_M2<-scanData('n_proportion_m2.in')

b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:(years[2]),quarter=1:nq,age=fa:la)
b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
length(WSEA);length(PROPMAT);length(M);length(M1);length(PROP_M2);dim(b)[[1]]
b<-data.frame(b,WSEA=WSEA, PROPMAT=PROPMAT,M=M,M1=M1,PROP_M2=PROP_M2)
b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))

out<-c(out,list(bio=b))
################  other_sp

if (otherPredExist) {
  WSEA<-scanData('west.in')
  WSEA<-WSEA[1:((first.VPA-1)*noAreas*ny*(la-fa+1)*nq)]
  length(WSEA)
  N<-scanData('other_pred_n.in')
  b<-expand.grid(sub_area=1:noAreas,species.n=1:(first.VPA-1),year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  length(N)
  if ( length(WSEA) !=length(N) ) stop('different number of observations in west.in for "other predators") and other_pred_n.in')

  b<-data.frame(b,WSEA=WSEA, N=N)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,species.n,age,quarter,sub_area,WSEA,N))
  b$N<-b$N*catchMultiplier

  out<-c(out,list(other=b))

  ###################### mean l

  l<-scanData('lsea.in')
  b<-expand.grid(sub_area=1:noAreas,species.n=1:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
   b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,mean_l=l)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,age,quarter,sub_area,mean_l))
  out<-c(out,list(mean_l=b))

  ###########################  consum

  CONSUM<-scanData('consum.in')
  b<-expand.grid(sub_area=1:noAreas,species.n=1:npr,year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,CONSUM=CONSUM)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,sub_area,CONSUM))
  out<-c(out,list(consum=b))
}
 return(out)
}


