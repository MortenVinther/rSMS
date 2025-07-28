make_rsms_data<-function(dir=data.path,rsmsFile='rsms.dat',fleetIn="fleet_info.dat",multi=FALSE,stom_data=file.path(dir,'stom.dat'),alk_data=file.path(dir,'alk.dat')) {
 
 # TEST single species and one species
#  dir=data.path; rsmsFile='rsms.dat';multi=FALSE
  
  
 #TEST   dir=data.path; rsmsFile='rsms.dat';multi=FALSE; fleetIn="fleet_info.dat"; multi=TRUE; stom_data=file.path(data.path,'stom_input','STOM_model_new.dat') ;alk_data=file.path(data.path,'stom_input','ALK_stom_new.dat')

if (TRUE) { 
sms<-read.RSMS.control(dir=dir,file=rsmsFile,test=FALSE)
# write.RSMS.control(sms,file=file.path(data.path,outfile=rsmsFile))

info<-sms@species.info
first.VPA<-which(info[,'predator'] %in% c(0,1))[1]

oldAllSpeciesNames<-rownames(info)
oldPredNames<-rownames(info[info[,'predator']>=1,,drop=FALSE])
oldPredNo<-match(oldPredNames,rownames(info))


incl_other<- any(info[,'predator']==2)

if (incl_other) {
  off.species<-as.integer(-first.VPA+1) 
  off.oths<-as.integer(first.VPA+1)  #other species offset # in rsms, other species come after VPA (analytical)  species
  info<-rbind(info[first.VPA:sms@no.species,,drop=FALSE],info[1:(first.VPA-1),,drop=FALSE])
  nOthSpecies<-as.integer(first.VPA-1)
} else {
  off.species<-0L
  off.oths<- -9L
  info<-info[first.VPA:sms@no.species,,drop=FALSE]
  nOthSpecies<-0L
}

info<-cbind(s=1L:dim(info)[[1]],info)

off.season<-0L # not used 

nsp<-sms@no.species
nSpecies<-sms@no.species-first.VPA+1L  # species with analytically assessment
otherFoodn<-nSpecies+1L
nAllSp<-nSpecies+nOthSpecies

nPredator<-sum(info[,'predator']>0)

predators<-info[info[,'predator']>0,'s']

preys<-info[info[,'prey']>0,'s']
isPrey<-info[1:nSpecies,'prey']>0
if (!multi) isPrey[]<-FALSE

oldNewSp_n<-c((1:nOthSpecies)+nSpecies,1:nSpecies)

ages<-sms@first.age:sms@max.age.all
nAges<-length(ages)
recSeason<-as.integer(sms@rec.season)
nSeasons<-as.integer(sms@last.season)
extendModel<- sms@extend.model

years.model<-sms@first.year.model:(sms@last.year.model+extendModel); nYears.model<-length(years.model)
years.data<-sms@first.year.data:sms@last.year.data; nYears.data<-length(years.data)

yq_idx<-matrix(1:(nYears.model*nSeasons),nrow=nYears.model,ncol=nSeasons,byrow=TRUE)
off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)

fq<-1L
lq<-nSeasons
fa<-sms@first.age
fqa<-rep(1L,max(info[,'last-age'])+off.age); fqa[1]<-recSeason
faq<-rep(1L,nSeasons); if (recSeason>1) faq[1:(recSeason-1)]<-2

info<-cbind(info,fModel=c(sms@fModel,rep(-9,nOthSpecies)))
fSepar<-if_else(sms@fModel %in% c(2,3),sms@firstAgeYearEffect+off.age,99L)
info<-cbind(info,fSepar=c(fSepar,rep(-9,nOthSpecies)))

inputCbySeason<-info[1:nSpecies,"C by timestep"]

spNames<-dimnames(info)[[1]][1:nSpecies]
if (nOthSpecies>0) {
  othspNames<-dimnames(info)[[1]][(nSpecies+1):nAllSp] 
  allSpNames<-c(spNames,othspNames)
} else {
  allSpNames<-spNames
  othspNames<-as.character(NULL)
}
allSpNamesLong<-sms@species.names.long[ match(allSpNames,oldAllSpeciesNames) ]

predNames<- dimnames(info[info[,'predator']>0,,drop=FALSE])[[1]]
preyNames<- dimnames(info[info[,'prey']>0,,drop=FALSE])[[1]]
otherFoodName<-'OTH'

# zero-catch years
zy<-data.frame(s=1:nSpecies,y=0L)
bz<-expand.grid(s=(first.VPA:nsp) +off.species,y=1:nYears.data,q=1:nSeasons)
bz<-bz[order(bz$s,bz$y,bz$q),]

if (sms@zero.catch.year.season==1) {
  fname<-'zero_catch_year_season.in'
  z<-head(scan(file.path(dir,fname),comment.char='#',quiet = TRUE),-1)
  bz$z<-z
} else bz$z<-1
zeroCatchYearSeason<-bz
zeroCatchYear<-bz %>%group_by(s,y) %>% summarize(z=sum(z),.groups = "drop") %>%  filter(z==0) %>% mutate(z=NULL)
zy<-rbind(zy,zeroCatchYear)



zeroCatchYear<-by(zy,list(zy$s),function(x) { x<-x$y[x$y>0]})
zeroCatchYearExistsSp<-unlist(lapply(zeroCatchYear,function(x) length(x)))
zeroCatchYearExistsSp[zeroCatchYearExistsSp>0] <-1L

zy<-zy %>% mutate(zero=TRUE) %>% filter(y>0)
zeroCatchYear.df<-zy
if (dim(zy)[[1]]>0) zeroCatchYearExists<-1L else zeroCatchYearExists<-0L

zeroCatchYearSeason<-zeroCatchYearSeason %>% arrange(s,y,q,z)
#xtabs(z~y+q,data=filter(zeroCatchYearSeason,s==8))

#zero-catch season-age
bz<-expand.grid(s=(first.VPA:nsp) +off.species,a=1:nAges,q=1:nSeasons)
bz<-bz[order(bz$s,bz$q,bz$a),]

if (sms@zero.catch.year.season==1) {
  fname<-'zero_catch_season_ages.in'
  z<-head(scan(file.path(dir,fname),comment.char='#',quiet = TRUE),-1)
  bz$z<-as.integer(z) 
} else  bz$z<-1  
bz<-left_join(bz,data.frame(s=1:nSpecies,FFA=info[1:nSpecies,"first-age F>0"]+off.age),by = join_by(s)) %>% mutate(z=if_else(a>=FFA,z,0),FFA=NULL) %>%
        mutate(z=if_else(q<fqa[a] & a==fa+off.age,0,z))


zeroCatchSeasonAge.df<-bz
zeroCatchSeasonAge<-tapply(bz$z,list(bz$s,bz$q,bz$a),sum)
# ftable(zeroCatchSeasonAge[1,,])


stockRecruitmentModelCode<-sms@SSB.R
names(stockRecruitmentModelCode)<-spNames
stockRecruitmentModelCodeAdd<-sms@SSB.R.add
names(stockRecruitmentModelCodeAdd)<-spNames



use_rho<-as.integer(sms@use.rho)
#use_rho<-factor(use_rho, levels = c(0L,1L,2L),labels = c("no_correlation", "compound_symmetry", "AR(1)"))

fbarRange<-matrix(as.integer(sms@avg.F.ages),ncol=2,dimnames=dimnames(sms@avg.F.ages))+off.age;
colnames(fbarRange)<-c('first-a','last-a')

seasonalCatches<-as.integer(sms@combined.catches==3)    # are seasonal catches used in the catch observation likelihood
FfromSeparableModel<-as.integer(sms@combined.catches %in% c(3,4)) # are seasonal F derived from seasonal catch observations

## configuration of parameter keys

rho<-rep(0.7,nSpecies)
# 
# ##### states at age for F random walk,
# # use sms@keyLogFsta
# #  first ages in each group by species 
# 1 2 3 5	#  COD 
# 
### First age in separable model (in case of fModel==2|3), (option firstAgeYearEffect))
#    COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
#3    99    99    99    99    99    99    99    99    99    99    99 

upAge<-fSepar-off.age;
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  i<-1L
  la<-info[s,"last-age"]
  aa<-sort(unique(c( sms@keyLogFsta[[s]][sms@keyLogFsta[[s]]<=upAge[s]] ,la+1)));
  for (j in (1:(length(aa)-1))) {
    for (age in aa[j]:(aa[j+1]-1)) {
      x[s,age+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=age,a=age+off.age,keyLogFsta=i))
    }
    i<-i+1L
  }
}

keyLogFsta<-x
keyLogFsta.df<-xx

keyLogFsta.list<-by(keyLogFsta.df,keyLogFsta.df$s,function(x) as.vector(x$keyLogFsta)-x[1,"keyLogFsta"]+1 )

# not used  keyLogFstaSeason<-lapply(1:nSeasons,function(q) keyLogFsta)


useFrandomWalk<- ! info[1:nSpecies,'fModel']==3

nlogF<-apply(keyLogFsta,1,max)
#nlogF<-unlist(lapply(sms@keyLogFsta,length))
nlogF[!useFrandomWalk]<-0
nlogFto<-cumsum(nlogF)
nlogFto[!useFrandomWalk]<-NA
nlogFfrom<- c(1,head(cumsum(nlogF),-1)+1 )
nlogFfromTo<-matrix(c(as.integer(nlogFfrom),as.integer(nlogFto)),ncol=2)
nlogFfromTo[is.na(nlogFfromTo[,2]),]<- 0

keyLogFstaSd<-keyLogFsta
keyLogFstaSd[keyLogFstaSd<0]<- -9999L
keyLogFstaSd<-keyLogFstaSd+nlogFfrom-1
keyLogFstaSd[keyLogFstaSd<0]<- -1L

keyLogFstaSd[!useFrandomWalk,] <--1L

if (max(keyLogFstaSd)>0) logSdLogFsta<-rep(log(0.2),max(keyLogFstaSd)) else logSdLogFsta<-numeric(0)  # parameter


x<-as.integer(!useFrandomWalk)
idxLogYearEffectF<-cumsum(x)
idxLogYearEffectF[useFrandomWalk]<-0L

logYearEffectF<-  matrix(log(1),ncol=nYears.model,nrow=max(idxLogYearEffectF))
rownames(logYearEffectF)<-spNames[!useFrandomWalk]
colnames(logYearEffectF)<-as.character(years.model)


# useFrandomWalk


##### year and Ages in separable F-Model
# use sms@catch.sep.year and sms@catch.sep.season

x<-map2(sms@keyLogFsta,sms@catch.sep.year,function(x,y) {
  xx<-expand.grid(age=x,year=y)
})

for (s in 1:nSpecies) x[[s]]<-data.frame(x[[s]],s=s,fModel=info[s,'fModel'],fSepar=fSepar[s])
ay<-do.call(rbind,x) %>% filter(fModel %in% c(2,3)) %>%arrange (s,age,year) 
#filter(ay,s==1)
if (dim(ay)[[1]]>0) {
  ay<-  ay %>%  mutate(y=year-sms@first.year.model+1,a=age-sms@first.age+1) %>% filter(a>=fSepar) 
  ay$index<-1:dim(ay)[[1]]
}

#filter(ay,s==1)
x<-lapply(1:nSpecies,function(x) {
  matrix(-1L,ncol=info[x,'last-age']-sms@first.age+1L,nrow=nYears.model,dimnames=list(years.model,paste('age',sms@first.age:info[x,'last-age'])))
})

if (dim(ay)[[1]]>0) for (i in (1:dim(ay)[[1]])) {
  x[[ay[i,'s']]][ay[i,'y']:nYears.model,ay[i,'a']]<-ay[i,'index'] 
}

x<-lapply(x,function(xx){
  for (a in 2:ncol(xx)) if (xx[1,a]== -1) xx[,a]<- xx[,a-1]
  xx
})
#x[[1]]

keylogSeparAgeF<- x
if (dim(ay)[[1]]>0) logSeparAgeF<-rep(log(0.1),max(ay$index)) else logSeparAgeF<-numeric(0)


####################

x<-map2(sms@catch.sep.season,sms@catch.sep.year,function(x,y) {
  xx<-expand.grid(age=sort(unique(c(fa,x))),year=y,q=fq:lq)
})

#  filter(x[[6]],age==0)
for (s in 1:nSpecies) {
  x[[s]]<-data.frame(x[[s]],s=s,seasonalC=seasonalCatches[s],FfromSeparable=FfromSeparableModel[s]) %>%
    filter(FfromSeparable==1) %>%   mutate(y=year+off.year,seasonalC=NULL,FfromSeparable=NULL,a=age+off.age) 
}
#head(zeroCatchSeasonAge.df)
xx<-do.call(rbind,x)

xx<-left_join(xx,zeroCatchSeasonAge.df,by = join_by(q, s, a)) %>% mutate(qGrp=if_else(z==1,1L,0L)) %>%
    arrange (s,y,q,age) %>%   mutate(qGrp=cumsum(qGrp)) %>%mutate(qGrp=if_else(z==1,qGrp,0L)) %>%rename(fy=y,fa=age)

seasonalFdataOver<- xx %>% filter(qGrp>0) %>%transmute(s,fy,q,a,qGrp) %>% arrange(qGrp)


# ftable(xtabs(qGrp~year+q+a,data=filter(xx,s==7)))

if (dim(xx)[[1]] >0) p<-max(xx$qGrp) else p<-0
if (p>0)  logFSeasonal<-rep(log(1/nSeasons),p)  else  logFSeasonal<-numeric(0)  # parameter


allY<-do.call(rbind,lapply(1:length(sms@catch.sep.year),function(s){
  y<-sms@catch.sep.year[[s]]
  y<-sort(unique(c(y,nYears.model-off.year+1L)))+off.year
  do.call(rbind,lapply(1:(length(y)-1L),function(i) expand.grid(s=s,fy=y[i],y=y[i]:(y[i+1]-1L))))
}))
#sms@catch.sep.year[[9]]
#filter(allY,s==9)

allA<-do.call(rbind,lapply(1:length(sms@catch.sep.season),function(s){
  a<-sort(unique(c(fa,sms@catch.sep.season[[s]])))
  a<-sort(unique(c(a,info[s,'last-age']+1L)))
  do.call(rbind,lapply(1:(length(a)-1L),function(i) expand.grid(s=s,fa=a[i],age=a[i]:(a[i+1]-1L))))
}))

x2<-full_join(xx,allY,by = join_by(s, fy),relationship = "many-to-many")
x2<-full_join(x2,allA,by = join_by(s, fa),relationship = "many-to-many") %>%
        filter(!is.na(q)) %>% mutate(fa=NULL,fy=NULL,a=age+off.age,age=NULL) %>%
        arrange(s,y,q,a)
#head(x2)

keylogSeasonF<-lapply(1:nSpecies,function(x) {
  y<-filter(x2,s==x)
  tapply(y$qGrp,list(y$y,y$q,y$a),sum)
})

if (FALSE) {
  keylogSeasonF[[6]][1,,]
  keylogSeasonF[[7]][1,,]
  keylogSeasonF[[8]][1,,]
  keylogSeasonF[[9]][1,,]
  keylogSeasonF[[10]][1,,]
}


##### states at age for N ,

nlogN<-as.integer(info[1:nSpecies,"last-age"]-sms@first.age+off.age)

#### Process Noise N
# use sms@keyVarLogN
inclProcess<-sms@incl.process.noise; names(inclProcess)<-spNames
# inclProcess[c(2:3,8)]<-2 ; inclProcess[c(5,7)]<-3; inclProcess 
#    1= Both recruits and older ages have process noise estimated
#    2= Recruits only (N at age older than recruits in the first year become parameters)
#    3= Do not include process noise at all (N in the first year and recruits become parameters)


chk<-stockRecruitmentModelCode==0 & inclProcess==3
if (any(chk)) {
  cat("You cannot have recruitment as a fixed parameter (incl.process.noise=3) and a random walk recruits (recruitment model=0)\n",
      "incl.process.noise:\t",inclProcess,"\n",
      "recruitmet model:  \t",stockRecruitmentModelCode,'\n',
      "problem:",chk,'\n')
  stop()
}
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age"]
  aa<-sort(unique(c(sms@keyVarLogN[[s]],la+1)))
  if (inclProcess[s]==2) aa<-aa[1] else if (inclProcess[s]==3)  aa<-numeric(0)
  if (inclProcess[s]==1) {
    for (j in (1:(length(aa)-1))) {
      for (age in aa[j]:(aa[j+1]-1)) {
        x[s,age+off.age ]<-i;
      }
      i<-i+1L
    }} else if (inclProcess[s]==2) {
      x[s,off.age ]<-i;
      i<-i+1L  
   }
}
#inclProcess;x
keyVarLogN<-x


aa<-sort(unique(as.vector(keyVarLogN))) ;aa<-aa[aa>0]
logSdLogN<-rep(log(0.2),length(aa))  


nlogN[inclProcess==2]<-1L # recruits only
nlogN[inclProcess==3]<-0L

nlogNto<-cumsum(nlogN)  
nlogNfrom<- c(1,head(cumsum(nlogN),-1)+1L )
nlogNfromTo<-matrix(c(as.integer(nlogNfrom),as.integer(nlogNto)),ncol=2)
nlogNfromTo[nlogNfromTo[,1]>nlogNfromTo[,2]]<-c(0L,0L);
rownames(nlogNfromTo)<-spNames


inclSsbR<-(sms@incl.Ssb.R>0) & inclProcess==3
logSsbRsd<-rep(0.3,sum(inclSsbR>0))
tmp<-cumsum(inclSsbR)
tmp[!inclSsbR]<-0L
inclSsbR<-tmp
SsbRweight<-sms@incl.Ssb.R

#### N in first year or recruits parameter, if needed
nlogN2<-as.integer(info[1:nSpecies,"last-age"]-sms@first.age+off.age)

logNfirstYparam<- numeric(0)
logNrecruitParam<-numeric(0)
logNfirstYparamfromTo<-matrix(0L,nrow=nSpecies,ncol=2); rownames(logNfirstYparamfromTo)<-spNames
logNrecruitParamfromTo<-logNfirstYparamfromTo
iniN<-log(1E7)

for (s in 1:nSpecies) {
  if (inclProcess[s] %in% c(2,3))  {
     for (a in (2L:nlogN2[s])) logNfirstYparam<-c(logNfirstYparam,iniN-a)
     logNfirstYparamfromTo[s,1]<-max(logNfirstYparamfromTo)+1L
     logNfirstYparamfromTo[s,2]<-max(logNfirstYparamfromTo)+nlogN2[s]-2L
  }
  if (inclProcess[s] %in% c(3))  {
    logNrecruitParam<-c(logNrecruitParam,rep(iniN,nYears.model))
    logNrecruitParamfromTo[s,1]<-max(logNrecruitParamfromTo)+1
    logNrecruitParamfromTo[s,2]<-max(logNrecruitParamfromTo)+nYears.model-1

  }  
}



#### variance of catch and surveys     sms@seasonal.catch.s2

# FfromSeparableModel
# seasonalCatches

# first catches
a<-do.call(rbind,lapply(1:nSpecies,function(s) data.frame(s=s,vari=TRUE,a=sort(unique(c(sms@catch.s2.group[[s]]+off.age))))))
b<-do.call(rbind,lapply(1:nSpecies,function(s) data.frame(s=s,array2DF(zeroCatchSeasonAge[s,,])) %>% transmute(s,q=as.integer(Var1),a=as.integer(Var2),incl=Value)))

ab<- left_join(b,a,by = join_by(s, a))

a<-data.frame(s=1:nSpecies,seasonal.catch.s2=sms@seasonal.catch.s2)
ab<- left_join(ab,a,by = join_by(s))

ft<-data.frame(s=1:nSpecies,faf=info[1:nSpecies,'first-age F>0']+off.age,laf=info[1:nSpecies,'last-age']+off.age,seasonalCatches=seasonalCatches)

ab<-left_join(ab,ft,by = join_by(s)) %>% filter(a>=faf & a<=laf & incl) %>% filter( !(a==1 & seasonalCatches & q <fqa[1])) %>%
  group_by(s,a) %>% mutate(firstQ=min(q)) %>%
  filter((q==firstQ & seasonalCatches==0 ) | (q==firstQ & seasonal.catch.s2==0) |( seasonalCatches==1 & seasonal.catch.s2==1)) %>% mutate(q=if_else( seasonalCatches==1 & seasonal.catch.s2==1,q,1)) %>% 
  transmute(s,q,a,vari=if_else(is.na(vari),FALSE,TRUE)) %>% ungroup() %>%
  arrange(s,q,a) %>% mutate(varGroup=cumsum(vari))

keyVarObsCatch.df<-ab
logSdLogObsCatch<-rep(log(0.5),max(keyVarObsCatch.df$varGroup)) # parameter

#####



###########  survey

a<-readFleetInfo(dir=dir,of=fleetIn,off.age=off.age,off.year=off.year,nSpecies=nSpecies,verbose=F) 
fleetControl<-a


minCVsurvey<-a$cv
keySurvey<-a[['k']]
notUse<-keySurvey[keySurvey[,'useFleet']==0,'f']
catchability<-a[['catchability']]
catchability.var<-a[['catchability.var']]

if (length(notUse>0)) {
  catchability[notUse]<-NULL
  catchability.var[notUse]<-NULL
  keySurvey<-keySurvey[-notUse,]
  keySurvey[,'f']<-1L:dim(keySurvey)[[1]]
}

nFleets<-nrow(keySurvey)
fleetNames<-rownames(keySurvey)

keySurvey.df<-as.data.frame(keySurvey) %>% tibble::rownames_to_column("fName")

## survey variance

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(catchability.var[[f]],keySurvey[f,'maxage']+1)+off.age))
  s<-keySurvey[f,'s']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,a=j,age=j-off.age,keyVarObsSurvey=i))
    }
    i<-i+1L
  }
}

keyVarObsSurvey<-x
keyVarObsSurvey.df<-xx

logSdLogObsSurvey<-rep(log(0.4^2),max(keyVarObsSurvey))

## survey catchability
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(catchability[[f]],keySurvey[f,'maxage']+1)+off.age))
  s<-keySurvey[f,'s']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,a=j,age=j-off.age,keyCatchability=i))
    }
    i<-i+1L
  }
}
keyCatchability<-x
keyCatchability.df<-xx

logCatchability<-rep(0.0,max(keyCatchability)) 

## Catchability power age
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<- 0L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keySurvey[f,'s']
  for (j in (keySurvey[f,'mina']:keySurvey[f,'maxa'])) {
    if (j==keySurvey[f,'PowerAge']) {i<-i+1L; x[f,j ]<-i} 
    xx<-rbind(xx,data.frame(f=f,s=s,species.n=s,age=j-off.age,a=j,keyPowerQ=i))
  }
}

keyQpow<-x
keyQpow.df<-xx
if (max(keyQpow) >0) logQpow<-rep(0.0,max(keyQpow)) else  logQpow<-rep(0.0,0)


### survey observations
cpue<-readFleetCatch(fleetControl,dir,of=fleetIn,so='fleet_catch.in', condense=TRUE,nSpecies=nSpecies,verbose=FALSE)
if (any(sms@zero.survey.catch==1)) {
  cpue<-left_join(cpue,data.frame(s=1L:nSpecies,subZero= sms@zero.survey.catch),by = join_by(s)) %>% mutate(isZero=obs==0)
  filter(cpue,isZero)
  minC<- filter(cpue,obs>0) %>% group_by(s,f,age) %>% summarize(minC=min(obs),.groups='drop')
  cpue<-left_join(cpue,minC,join_by(s, f,age )) %>% mutate(obs=if_else(isZero & subZero==1,minC,obs),subZero=NULL,isZero=NULL,minC=NULL)
}

eff_fl<-keySurvey[keySurvey[,'type'] %in% c(4,5),,drop=FALSE]

if (dim(eff_fl)[[1]] > 0) { #inflate data with effort for each age group for type=effort
  ef<-  eff_fl[,'f']
  cpue_ok<-filter(cpue,!(f %in% ef)) 
  cpue_eff<-filter(cpue,f %in% ef) 
  a<-lapply(ef,function(ff){
    a<-filter(cpue_eff,f==ff)
    a<-a %>% mutate(obs=obs/mean(obs)) 
    ages<-data.frame(age=-999,newage=keySurvey[ff,"minage"]:keySurvey[ff,"maxage"])
    left_join(a,ages,by = join_by(age),relationship = "many-to-many")  %>% mutate(age=newage,newage=NULL)
  })
  cpue_eff<-do.call(rbind,a)
  cpue<-rbind(cpue_ok,cpue_eff)
}

if (FALSE) {
  out<-by(cpue,cpue$f,function(x){
    (xtabs(obs~age+year,data=droplevels(x)))
  })
  out[[1]]
  out[[2]]
}
cpue<-filter(cpue,!(is.na(obs)))
cpue<-cpue %>% filter(year<=sms@last.year.model+sms@extend.model)   # SKAL RETTES


cpue<-cpue %>% mutate(y=year+off.year,a=age+off.age,species.n=s,season=q,fleet.no=f) %>%
  dplyr::select(year,y,season,q,age,a,species.n,s,fleet.no,f,obs) %>%
  filter(obs>0) %>% arrange(s,fleet.no,y,a) 
 cpue$obs.no<-1L:dim(cpue)[[1]]
#head(cpue)
 if (FALSE) {
   out<-by(cpue,cpue$f,function(x){
     (xtabs(obs~age+year,data=droplevels(x)))
   })
   out[[1]]
   out[[2]]
   out[[18]]
   out[[19]]
   out[[20]]
 }

fobs<-cpue %>% dplyr::select(s,f,y,a) %>% group_by(s,f) %>%summarize(n=dplyr::n(),.groups = "drop") %>% mutate(last=cumsum(n),first=lag(last)+1L)
fobs[1,"first"]<-1L


keySurvey.overview<-cbind(keySurvey,n=as.vector(fobs$n),first=as.vector(fobs$first),last=as.vector(fobs$last))

sumCreep<-sum(keySurvey.overview[,'techCreep'])
logTechCreep<-rep(0.0,sumCreep)

k<-full_join(keyVarObsSurvey.df,keyCatchability.df,by = join_by(f, s, a, age)) %>% dplyr::select(f, s, a, keyVarObsSurvey, keyCatchability) %>% arrange(f,s,a)
k<-full_join(k,dplyr::select(keyQpow.df,s,f,a,keyPowerQ),by = join_by(f, s, a))

cpue<-left_join(cpue,k,by = join_by(a, s, f)) %>% arrange(f,s,a)
cpue$obs.no<-1:dim(cpue)[[1]]


keySurvey<-  cpue %>% dplyr::select(obs.no,f,s,y,a,q ,keyVarObsSurvey, keyCatchability,keyPowerQ)
#summary(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))
# print(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))

stopifnot(dim(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))[[1]]==0)
keySurvey<-keySurvey %>% mutate_if(is.numeric,as.integer)

keySurvey<-as.matrix(keySurvey)

logSurveyObs<-log(cpue$obs)


x<-scan(file=file.path(dir,'recruitment_years.in'),comment.char='#',quiet=TRUE)
stopifnot(tail(x,1)== -999)
rec_years<-matrix(head(x,-1),nrow=nSpecies,byrow=TRUE)==1

if (extendModel>0) rec_years<-cbind(rec_years,rep(FALSE,nSpecies))

source(file.path(rsms.root.prog,"from_sms_format_to_rsms_data.R"))


d<-From_SMS_format_to_rsms(dir,otherPredExist=multi,catchMultiplier=1) 


#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  d$other$WSEAoth<-d$other$WSEA; d$other$WSEA<-NULL; d$other$species.n.1<-NULL
  b2<-full_join(b2,d$other,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

a<-data.frame(species.n=first.VPA:nsp,annualC=sms@combined.catches,FfromSeparableModel=FfromSeparableModel,inputCbySeason)
b<-left_join(b,a,by = join_by(species.n))
#head(b)

if (extendModel==1) {
  years.model<-sms@first.year.model:(sms@last.year.model+extendModel); nYears.model<-length(years.model)
  years.data<-sms@first.year.data:sms@last.year.data; nYears.data<-length(years.data)
 if (sms@last.year.model+extendModel > sms@last.year.data ) {
    ext<-filter(b,year==sms@last.year.model) 
    for (y in sms@last.year.data:(sms@last.year.model+extendModel)) {
      ext$year<-y+1
     b<-rbind(b,ext)
    }
 }
}

#x<-b$age==fa & b$quarter < recSeason & b$annualC==0
 x<-b$age==fa & b$quarter < recSeason & b$inputCbySeason==1
b[x,'CATCHN']<-0
b[x,'WCATCH']<-0

#    1=annual catches in likelihood, with seasonal proportions of annual F at age from number of seasons (1/nseasons)
b1<-filter(b,annualC==1) %>%  dplyr::select(year, species.n, quarter, sub_area,age, annualC) %>% mutate(seasFprop=1/nSeasons) %>%
     mutate(seasFprop=if_else(age==fa & quarter >= recSeason,0.5, seasFprop)) %>%
     mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop))

#    2=annual catches in likelihood with seasonal proportions of annual F at age from file proportion_of_annual_f.dat 
b2<-filter(b,annualC==2) %>%  dplyr::select(year, species.n, quarter, sub_area,age, annualC,PROP_SEASON_F) %>% mutate(seasFprop=PROP_SEASON_F,PROP_SEASON_F=NULL) %>%
  filter((age==fa & quarter>=recSeason) | age>fa)

#    3=seasonal catch in likelihood with estimation of seasonal F from separable F-model 
b3<-filter(b,annualC==3) %>% group_by(year,species.n, sub_area, quarter,age) %>% mutate(seasFprop=CATCHN/sum(CATCHN,na.rm=T)) %>% ungroup() %>%
  mutate(seasFprop=if_else(is.na(seasFprop),0,seasFprop)) %>% 
  mutate(seasFprop=if_else(seasFprop>0,1,0)) %>%
  mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop)) %>%
  dplyr::select(year, species.n, quarter, sub_area,   age, annualC, seasFprop) 

#   4=annual catches in likelihood with estimation of seasonal F from separable F-model 
b4<-filter(b,annualC==4) %>% group_by(year,species.n, sub_area, quarter,age) %>% mutate(seasFprop=CATCHN/sum(CATCHN,na.rm=T)) %>% ungroup() %>%
  mutate(seasFprop=if_else(is.na(seasFprop),0,seasFprop)) %>% 
  mutate(seasFprop=if_else(seasFprop>0,1,0)) %>%
  mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop)) %>%
  dplyr::select(year, species.n, quarter, sub_area,   age, annualC, seasFprop) 



#filter(b,species.n==22)
bb<-rbind(b1,b2,b3,b4) 
head(bb)

b<-b %>% mutate(PROP_SEASON_F=NULL)

b<-left_join(b,bb,by = join_by(year, species.n, quarter, sub_area, age,annualC))


b<-b %>% mutate(y=year+off.year,q=quarter,a=age+off.age)  

#filter(b,species.n==21)

if (nSpecies>1 & incl_other) b<-b %>% mutate(s=oldNewSp_n[species.n]) else b<-b %>% mutate(s=species.n)  # change species numbering
b<-left_join(b,data.frame(s=info[,'s'],la=info[,"last-age"]+off.age),by = join_by(s)) %>% filter(a<=la)

if (multi) {b[b$s>nSpecies,'WSEA']<-b[b$s>nSpecies,'WSEAoth']; b$WSEAoth<-NULL}

propMat         <-by(b,b$s,function(x) {y<-tapply(x$PROPMAT,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
stockMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WSEA,list(x$y,x$q,x$a),sum) ; y[is.na(y)]<-0; y})
catchMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WCATCH,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
natMor          <-by(b,b$s,function(x) {y<-tapply(x$M,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
catchNumber     <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})

if (multi) {
  consum          <-by(b,b$s,function(x) {y<-tapply(x$CONSUM,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  meanL           <-by(b,b$s,function(x) {y<-tapply(x$mean_l,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  propM2          <-by(b,b$s,function(x) {y<-tapply(x$PROP_M2,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  natMor1         <-by(b,b$s,function(x) {y<-tapply(x$M1,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  otherN          <-by(b,b$s,function(x) {y<-tapply(x$N,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
} else {consum<-NA;meanL<-NA;propM2=NA; natMor1=NA;  otherN<-NULL}


trimIt<-function(x,used=c(1,2),spnames=spNames){
  x<-x[used]
  names(x)<-spnames
  return(x)
}

propMat         <-trimIt(x=propMat,used=1:nSpecies)
stockMeanWeight <-trimIt(x=stockMeanWeight,used=1:nSpecies)
natMor          <-trimIt(x=natMor,used=1:nSpecies)
catchMeanWeight <-trimIt(x=catchMeanWeight,used=1:nSpecies)
catchNumber     <-trimIt(x=catchNumber,used=1:nSpecies)

#  ftable(catchNumber[[7]][3,,])

if (multi) {
  consum          <-trimIt(x=consum,used=1:nAllSp,spnames=allSpNames)
  meanL           <-trimIt(x=meanL,used=1:nAllSp,spnames=allSpNames)
  propM2          <-trimIt(x=propM2,used=1:nSpecies)
  natMor1         <-trimIt(x=natMor1,used=1:nSpecies)
  otherN          <-trimIt(x=otherN,used=1:nAllSp,spnames=allSpNames)
  
  logPropM2<-propM2
  for (s in (1:nSpecies)) logPropM2[[s]][propM2[[s]]>0] <-log(logPropM2[[s]][propM2[[s]]>0])
  
  NN<-filter(b,!is.na(N)) %>% transmute(y,q,predNo=s,predAge=a,N) %>%filter(N>0)%>%mutate(N=NULL)
  
}
zero<-propMat; for (s in 1:nSpecies) {zero[[s]][,,]<-0}


#### catch observations


catch<-left_join(b %>% mutate(PROP_CAT=NULL,     WSEA=NULL, PROPMAT=NULL,      M=NULL,   M1=NULL, PROP_M2=NULL, species.n=NULL),
                 data.frame(s=1:nAllSp,fa=sms@first.age,combCat=c(sms@combined.catches,rep(-9,nOthSpecies)),recSeason=sms@rec.season,faf=info[,"first-age F>0"]+off.age),
                 by = join_by(s)) %>%
       filter(combCat %in% c(1,2) | (combCat %in% c(3,4) & (age>fa | (age==fa & q>=recSeason )))) 


if (FALSE) {
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==1 & y==1)))) #annual catch, faf=1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==2 & y==1))))  #annual catch, faf=0 and age 0 catch in Q1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==7 & y==1))))  #seasonal catch, no catch in Q1 and Q4
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==9 & y==1))))  #seasonal catch, no age 0 catch in Q1 and Q2
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==10 & y==2)))) #seasonal catch,no 0-group catch
}


bz<-zeroCatchSeasonAge.df
# filter(bz,s==10)
#head(catch);head(bz)
catch<-left_join(catch,bz,by = join_by(s,a, q)) %>% mutate(CATCHN   =if_else(z==0 & inputCbySeason==1,0,CATCHN),
                                                           seasFprop=if_else(z==0 & inputCbySeason==1,0,seasFprop),z=NULL)
if (FALSE) {
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==1 & y==1))),2) 
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==2 & y==1))),2)  
}  

zy<-zeroCatchYear.df
#filter(zy,s==6)
#head(catch);head(zy)
catch<-left_join(catch,zy,by = join_by(y, s)) %>% filter(is.na(zero)) %>% mutate(zero=NULL)
#filter(catch,s==2 & y==5)

if (FALSE) {
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==1 & y==1)))) #annual catch, faf=1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==2 & y==1))))  #annual catch, faf=0 and age 0 catch in Q1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==7 & y==1))))  #seasonal catch, no catch in Q1 and Q4
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==8 & y==49))))  #seasonal catch, no catch in Q1 and Q4
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==9 & y==1))))  #seasonal catch, no age 0 catch in Q1 and Q2
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==10 & y==2)))) #seasonal catch,no 0-group catch
  
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==1 & y==1))),2) #annual catch, faf=1
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==2 & y==1))),2)  #annual catch, faf=0 and age 0 catch in Q1
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==7 & y==1))),2)  #seasonal catch, no catch in Q1 and Q4
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==9 & y==1))),2)  #seasonal catch, no age 0 catch in Q1 and Q2
  round(ftable(xtabs(seasFprop~year+q+age,data=filter(catch,s==10 & y==2))),2) #seasonal catch,no 0-group catch
}

catch<-left_join(catch,zeroCatchYearSeason,by = join_by(y, q,s)) %>% filter(z==1) %>% mutate(z=NULL)

#  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch2,s==2))))
#round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==8 & y==49))))  #seasonal catch, no catch in Q1 and Q4

catch2<-catch

catch<-catch %>% mutate(q=if_else(combCat %in% c(3),q,1))
#summary(filter(catch,s==1))
#summary(filter(catch,s==6))



catch<-aggregate(CATCHN~year+ sub_area+ age+y+q+ s+ a,FUN=sum,data=filter(catch,CATCHN>0 & a>=faf)) %>% arrange(s,y,q,a)

if (FALSE) {
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==1 & y==1)))) #annual catch, faf=1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==2 & y==1))))  #annual catch, faf=0 and age 0 catch in Q1
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==7 & y==1))))  #seasonal catch, no catch in Q1 and Q4
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==8 & y==49))))  #seasonal catch, no catch in Q1 Q3 and Q4
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==9 & y==1))))  #seasonal catch, no age 0 catch in Q1 and Q2
  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch,s==10 & y==2)))) #seasonal catch,no 0-group catch
}

catch$obs.no<-1:dim(catch)[[1]]
logCatchObs<-log(catch$CATCHN)

keyCatch<-catch %>% mutate(CATCHN=NULL,species.n=NULL) %>% mutate_if(is.numeric,as.integer)

## seasonal catches with only catches from one season only. No seasonal effect and logFSeasonal should be fixed 
#head(filter(keyCatch,s==1),10)
#head(filter(keyCatch,s==2),10)

a<-left_join(keyCatch %>% transmute(q,s,a) %>% unique(),data.frame(s=1:nSpecies,FfromSeparableModel),by = join_by(s)) %>% filter(FfromSeparableModel==1) %>% mutate(FfromSeparableModel=NULL)
a<-a %>% group_by(s,a) %>% summarize(n=dplyr::n(),.groups='drop') %>% filter(n==1) %>% mutate(n=NULL)
a$fix<-T

bfix<-left_join(seasonalFdataOver,a,by = join_by(s, a)) %>% mutate(fix=if_else(is.na(fix),FALSE,fix)) %>% arrange(qGrp)

if (any(bfix$fix)) {
  fixMapseasonalF<-b$qGrp[bfix$fix]
  logFSeasonal[fixMapseasonalF]<-0.0
} else fixMapseasonalF<- -9


# add keyLogFsta
keyCatch<-left_join(keyCatch,keyLogFsta.df,by = join_by(age, s, a))

# add variance group
# catches by season, but common variance
tmp<-left_join(data.frame(s=1:nSpecies,seasonal.catch.s2=sms@seasonal.catch.s2,seasonalCatches),expand.grid(s=1:nSpecies,q=fq:lq),by = join_by(s),relationship = "many-to-many") %>%
  filter(seasonal.catch.s2==0 & seasonalCatches==1) %>% transmute(s,q)
qq<-left_join(keyVarObsCatch.df,tmp,by = join_by(s),relationship = "many-to-many") %>% mutate(q=if_else(is.na(q.y),q.x,q.y),q.x=NULL,q.y=NULL, vari=NULL)
                                                                                   
#print(filter(qq,s==9 | s==6),n=70)

keyCatch<-left_join(keyCatch,qq,by = join_by(s, q,a)) %>% mutate(species.n=NULL,age=NULL) %>% rename(keyVarObsCatch=varGroup)  %>% arrange(obs.no)
#filter(keyCatch,s==9 & year==2022)

#print('Check it Morten !')
stopifnot(dim(filter(keyCatch,is.na(keyLogFsta)))[[1]]==0)

keyCatch <-keyCatch %>% mutate(year=NULL,keyLogFsta=NULL)

keyCatch<-as.matrix(keyCatch)

#print( tail(keyCatch[keyCatch[,'s'] %in% 8,],10))


#presence of catches by y,q and a, and if so the proportion of F, does not include years with no catches


#round(ftable(xtabs(seasFprop~y+q+a,data=filter(catch2,s==1 & y==1))),2) 
#round(ftable(xtabs(seasFprop~y+q+a,data=filter(catch2,s==2 & y==1))),2) 
y<-tapply(catch2$seasFprop,list(catch2$s,catch2$y,catch2$q,catch2$a),sum) # does not include years with no catches
y[is.na(y)]<-0
# y[1,1,,]; y[2,1,,]

y<-array2DF(y) %>% mutate_if(is.character,as.integer) %>% transmute(s=Var1,y=Var2,q=Var3,a=Var4,seasFprop=Value)


if (extendModel>=1) {
  ext<-filter(y,y==nYears.data) 
  for (yy in ((nYears.data+1):nYears.model)) {
    ext$y<-yy
    y<-rbind(y,ext)
  }
}

#round(ftable(xtabs(seasFprop~y+q+a,data=filter(y,s==1 & y==1))),2) 
#round(ftable(xtabs(seasFprop~y+q+a,data=filter(y,s==2 & y==1))),2)

seasFprop <-by(y,y$s,function(x) {y<-tapply(x$seasFprop,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y}) 

# seasFprop[[s]][1,,]; seasFprop[[s]][2,,]

seasFprop2 <-trimIt(x=seasFprop,used=1:nSpecies)

# seasFprop2[[s]][1,,]; seasFprop2[[s]][2,,]

logSeasFprop<-seasFprop
for (s in (1:nSpecies)) logSeasFprop[[s]][seasFprop[[s]]>0] <-log(logSeasFprop[[s]][seasFprop[[s]]>0])


catchNoSeason   <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
catchNoSeasonUsed   <-by(catch2,catch2$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})



info<-cbind(info,
    la=info[,"last-age"]+off.age,
    faf=info[,"first-age F>0"]+off.age)
 
##### recruitment

recEst<-unlist(lapply(catchNumber,function(x) sum(unlist(x))/nYears.model*2))
names(recEst)<-NULL

#      1=Ricker, 2=Beverton & Holt, 3=Geom mean,
#      4= Hockey stick, 5=hockey stick with smoother,
#      6= hockey stick with known breakpoint (given as input in SSB.R.add )

rec_loga <- log(recEst)           # geom mean
rec_loga[stockRecruitmentModelCode==1] <- log(200)  # Ricker
rec_loga[stockRecruitmentModelCode==2] <- 4         # B&W
rec_loga[stockRecruitmentModelCode %in% c(4,5,6)] <- log(1E6) #hockey slope

rec_logb <- rep(1,nSpecies) 
rec_logb[stockRecruitmentModelCode==1] <- -12  # Ricker
rec_logb[stockRecruitmentModelCode==2] <- -12        # B&W
rec_logb[stockRecruitmentModelCode %in% c(4,5) ] <- log(50000)     
rec_logb[stockRecruitmentModelCode ==6 ] <- log(stockRecruitmentModelCodeAdd[stockRecruitmentModelCode ==6] ) #hockey with known inflection point
}


###############   stomach data

if (multi) {
  
  oldNewPred <-function(opt,label,info1=info) {
    x<-  slot(sms,opt)
    names(x)<-oldPredNames
    x<-x[predNames]
    info1<-cbind(info1,xx= -9)
    info1[rownames(info1) %in% predNames,'xx']<-x
    colnames(info1)<-c(head(colnames(info1),-1),label)
    return(info1)
  }

  info<-oldNewPred(opt='size.selection',label='sizeSlct')
  info<-oldNewPred(opt='sum.stom.like',label='sumStomLike')
  info<-oldNewPred(opt='stom.max.sumP',label='stomMaxSumP')
  info<-oldNewPred(opt='max.stom.sampl',label='stomMaxSamp')
  info<-oldNewPred(opt='size.other.food.suit',label='sizeOtherFoodSuit')
  info<-oldNewPred(opt='stomach.variance',label='stomachVariance')

    
  s<-read_delim(stom_data,show_col_types = FALSE)   %>% arrange(year,quarter,pred.no,pred.size.class,prey.no,prey.size.class)

 
 if (FALSE) { #use stomach data from SMS output   
  l<-length(oldAllSpeciesNames)+1L
  s1<-read_table("summary_stom.out") %>% mutate(Prey.no=if_else(Prey.no==0,l,Prey.no),stom.type=if_else(stom.type==1 | is.na(stom.type),'obs','mid'),)  %>% 
    transmute(SMS_area=SMS.area,  year=Year,quarter=Quarter.no, pred=oldAllSpeciesNames[Predator.no], pred.no=Predator.no, pred.size=as.character(Predator.length), pred.size.class=Predator.length.class,
        pred.mean.length=Predator.length.mean, haul.no=N.haul, stom.no=N.samples, phi=Diri.sum.p,
        prey=c(oldAllSpeciesNames,otherFoodName)[Prey.no], prey.no=if_else(Prey.no==l,otherFoodn,Prey.no),  prey.size=as.character(Prey.length.mean),           
        prey.size.class=Prey.length.class,      prey.mean.length=Prey.length.mean,     prey.mean.length.ALK=-9, 
        stomcon,type=stom.type, mean.weight=Prey.weight,haul.prey.no =-9,       
        calc.prey.number=-9,used.prey.number=-9)  %>% arrange(year,quarter,pred.no,pred.size.class,prey.no,prey.size.class)
    
  #glimpse(s)
  #glimpse(s1)
  # s11<-full_join(s,s1,by = join_by(year,quarter,pred.no,pred,pred.size.class,prey.no,prey,prey.size.class))
  # glimpse(s11)
  
  s<-s1
}
  
  
# sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

#a<-filter(s,pred %in% c('RAZ','GLT'))
#round(ftable(tapply(a$stomcon,list(a$pred,a$prey,a$prey.size),sum)),1)

# North Sea fix
s[s$pred=='W_H','pred']<-'WHM'
s[s$pred=='N_H','pred']<-'NHM'

preds<-info[info[,'predator']>0,'s']
# sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey),prey.mean.length.ALK=NULL,pred.no=as.integer(pred.no),prey.no=as.integer(prey.no))
#xtabs(~tpred+tprey,data=s)  #OK but old species order

#s<-s %>% mutate(pred.no=oldNewSp_n[pred.no])
s$pred.no<-oldNewSp_n[s$pred.no]

# sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

excl<-!(s$prey %in% c(preyNames,otherFoodName))
#summary(excl)
soth<-as.data.frame(head(filter(s,prey==otherFoodName),1))
if (any(excl)) {
  s[excl,'prey']<-otherFoodName
  s[excl,'prey.no']<-soth$prey.no
  s[excl,'prey.size']<-soth$prey.size
  s[excl,'prey.size.class']<-soth$prey.size.class
  s[excl,'prey.mean.length']<-soth$prey.mean.length
  s[excl,'mean.weight']<-soth$mean.weight
  s[excl,'calc.prey.number']<-soth$calc.prey.number
  s[excl,'used.prey.number']<-soth$used.prey.number
  s[excl,'haul.prey.no']<-soth$haul.prey.no
  #names(s)
  #summary(excl)
}
#sort(unique(s$type))
# sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

a<-left_join(data.frame(pred=oldPredNames,typ=sms@stom.type.include), 
  rbind(
  data.frame(typ=1,type=c('obs')),
  data.frame(typ=2,type=c('obs','mid')),
  data.frame(typ=3,type=c('obs','mid','tail'))
 ),relationship = "many-to-many",by = join_by(typ))

#sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

s<-inner_join(a,s,by = join_by(pred, type))


#sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))
if (!('phi' %in% names(s))) s$phi<--9

s<-s %>% mutate(noSampl=if_else(haul.no>info[pred,'stomMaxSamp'],info[pred,'stomMaxSamp'],haul.no),
                    phi=if_else(haul.no>info[pred,'stomMaxSumP'],info[pred,'stomMaxSumP'],phi))

#filter(s,pred=='COD')
s<-s %>% group_by(SMS_area,year,quarter,pred,pred.no,pred.size,pred.size.class,pred.mean.length,noSampl,phi,
               prey,prey.no,prey.size,prey.size.class,prey.mean.length,type,mean.weight,haul.prey.no,calc.prey.number,used.prey.number) %>%
               summarize(stomcon=sum(stomcon),.groups = "drop")
                         
#sort(unique(paste(formatC(s$pred.no,w=2),s$pred)))

x2<-filter(s,prey=='OTH')
s<-filter(s,prey!='OTH') %>% mutate(prey.no=oldNewSp_n[prey.no])
s<-rbind(s,x2)
s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey))
# xtabs(~tpred+tprey,data=s) #OK


# vulnerability
pp<-matrix(0L,nrow=nSpecies+1,ncol=nSpecies+nOthSpecies,dimnames=list(c(spNames,'OTH'),c(spNames,othspNames)))
pp2<- s%>% dplyr::select(pred,pred.no,prey,prey.no) %>% mutate(prey.no=as.integer(ifelse(prey.no==0,nSpecies+1,prey.no))) %>%unique() 
#xtabs(~prey.no+pred.no,data=pp2)
for (i in (1:dim(pp2)[[1]])) {
  #cat(i,'\n'); 
  pp[unlist(pp2[i,'prey.no']),unlist(pp2[i,'pred.no'])]<-1L
}
predPrey<-pp
# predPrey

pp2$no<-1:dim(pp2)[[1]]

vulneraIdx<-predPrey
vulneraIdx[vulneraIdx>0]<-1:sum(vulneraIdx>0) 
# vulneraIdx
vulnera<-rep(0,max(vulneraIdx)) #log vulnerability coefficient, parameter 

vulneraIdx.df<-array2DF(vulneraIdx) %>% rename(prey=Var1,pred=Var2) %>% filter(Value>0)

# pred prey overlap
ov<-scan(file.path(dir,"season_overlap.in") ,comment.char = "#",quiet=TRUE)
ov<-head(ov,-1)
ov<-matrix(ov,ncol=nSpecies+1,nrow=nAllSp*4,dimnames=list(paste0(rep(c(spNames,othspNames),each=4),' Q',1:4), c(spNames,otherFoodName)),byrow=TRUE )
ov<-array(ov,dim=c(4,nAllSp,nSpecies+1),dimnames=list(paste0('Q',1:4),c(spNames,othspNames), c(spNames,otherFoodName)))                                                        
ovd<-array2DF(ov)
names1<-c( spNames,othspNames)
names2<-c(spNames,otherFoodName)
ovdt<-ovd%>%mutate(q=as.integer(substr(Var1,2,2)),predNo=match(Var2,names1),preyNo=match(Var3,names2)) %>% 
   transmute(q,predNo,preyNo,overlap=Value)
overlap<-ov
# ftable(overlap)
overlap[overlap>0] =log(overlap[overlap>0])
overlapIdx<-filter(ovdt,overlap<0) %>% unique() %>% mutate(overlap= -overlap)
overlapP<-rep(0,length(unique(overlapIdx$overlap))) # log overlap parameter

# put weight on from length
lw<-scan(file.path(dir,"length_weight_relations.in") ,comment.char = "#",quiet=TRUE)
lw<-head(lw,-1)
lw<-matrix(lw,ncol=2,nrow=nSpecies+nOthSpecies,dimnames=list(c(spNames,othspNames),c('a','b')),byrow=TRUE )
lwdf<-data.frame(lw) %>% mutate(sp=1L:dim(lw)[[1]])

stomRange<-data.frame(pred.no=oldPredNo,srL=sms@size.range.lower/100,srU=sms@size.range.upper/100,srMethod=sms@size.range.method)
stomRange$pred.no<-oldNewSp_n[stomRange$pred.no]
s<-left_join(s,stomRange,by = join_by(pred.no))

minStom<-data.frame(pred=predNames,minStom=sms@min.stom.cont)
s<-left_join(s,minStom,by = join_by(pred))

#observed predator prey size ratio
ss<-s%>% left_join(x=s,y=lwdf,by = join_by(pred.no == sp)) %>% 
  mutate(pred.size.w=a*pred.mean.length**b,a=NULL,b=NULL,prey.size.w=mean.weight,pred.prey.size.w=pred.size.w/prey.size.w)  
# filter(ss,pred=='COD' & year==1981 & quarter==1) %>% select(pred.no, pred.size ,pred.size.class ,pred.mean.length ,pred.size.w,prey.size.w,pred.prey.size.w)

a<-filter(ss,type=='obs' & prey !='OTH'  & stomcon>minStom) %>%  
  mutate(log.pred.prey.size.w=log(pred.size.w/prey.size.w),log.pred.size.w=log(pred.size.w)) %>%
  select(SMS_area,pred.no, pred,tpred,prey.no,prey,tprey,log.pred.size.w,log.pred.prey.size.w,srL, srU,srMethod) 

#summary(exp(a$log.pred.size.w))

k<-a %>% group_by(pred.no,pred,prey.no,prey,srMethod) %>% 
  summarize(min_s=quantile(log.pred.prey.size.w,probs=srL[1]),max_s=quantile(log.pred.prey.size.w,probs=srU[1]),n=dplyr::n(),.groups = "drop")  %>%
  filter(srMethod==0) 

# round(xtabs(exp(min_s)~pred.no+prey.no,data=k),1)
# round(xtabs(exp(max_s)~pred.no+prey.no,data=k),3)


#if (TRUE) { #quantile regression

# sms@size.range.method
  # constraint uniform size selection
  # first the predator size independent prey range 
  
  aa<-a %>% group_by(pred.no,pred,tpred,prey.no,prey,tprey) %>% summarize(min_s=min(log.pred.prey.size.w),max_s=max(log.pred.prey.size.w),.groups = "drop") 
  
  
  pp<-array(0,dim=c(4,nSpecies,nSpecies+nOthSpecies),dimnames=list(c('lower_intc','lower_slope','upper_intc','upper_slope'),spNames,c(spNames,othspNames)) )                                                         
  for (i in (1:dim(aa)[[1]])) {
    pred<-unlist(aa[i,'pred'])
    prey<-unlist(aa[i,'prey'])
    pp[1,prey,pred]<-unlist(aa[i,'min_s'])
    pp[3,prey,pred]<-unlist(aa[i,'max_s'])
  }
  #  round(ftable(pp),2)

  
  # add parameters from quantile regressions
  
  library(quantreg)
  #q1<-0.025 #lower quantile
  #q2<-0.975 #higher quantile
  
  
  if(FALSE) by(a,list(a$pred.no),function(x) {
    
    pp<-ggplot(data=x, aes(x=log.pred.size.w, y=log.pred.prey.size.w)) +
      geom_point()+
      labs(title=x[1,'pred'],y='log(predator size/prey size)',x='log(predator size)')+
      geom_smooth(method = "lm", se =F)+
      facet_wrap(vars(prey),ncol=3)
  })
    
  
  if (FALSE) {
    a1<-filter(a,pred.no %in% c(1:9) )
    aa<-by(a1,list(a1$prey.no,a1$pred.no),function(x,minNobs=5) {
      titl<-paste(x[1,'pred'],'eating',x[1,'prey'])
      print(titl)
      
      if (dim(x)[[1]]<minNobs) {
        print('too few observations')
      } else {
        ru<-rq(log.pred.prey.size.w~log.pred.size.w ,data=x, tau=q1)
        ru<-ru[["coefficients"]]
        rl<-rq(log.pred.prey.size.w~log.pred.size.w ,data=x, tau=q2)
        rl<-rl[["coefficients"]]
        
        
        pp<-ggplot(data=x, aes(x=log.pred.size.w, y=log.pred.prey.size.w)) +
          geom_point()+
          labs(title=titl,y='log(predator size/prey size)',x='log(predator size)')+
          geom_smooth(method = "lm", se =F)+
          geom_abline(slope=ru[2],intercept=ru[1])+
          geom_abline(slope=rl[2],intercept=rl[1])
        print(pp)
      } 
      if (dim(x)[[1]]<minNobs) {
        return(list(pred=unlist(x[1,'pred']),prey=unlist(x[1,'prey'])))
      } else return(list(pred=unlist(x[1,'pred']),prey=unlist(x[1,'prey']),ru=ru,rl=rl,plt=pp))
      
    })
  }
  
  # sort(unique(a$tpred))
  a1<-filter(a,srMethod==1 )
  doQuant<-TRUE
  if (dim(a1)[[1]]>1) {
    n<-a1 %>% group_by(pred.no,prey.no) %>% summarize(n=dplyr::n(),.groups = "drop")
    an<-left_join(a1,n,by = join_by(pred.no, prey.no))
    
    aa<-filter(an,n>55) %>% 
      nest_by(pred.no,pred,prey.no,prey,n,srL,srU) %>%  
      mutate(rl=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=srL,data=data)),
             ru=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=srU,data=data)))
    
    
    for (i in (1:dim(aa)[[1]])) {
      pred<-unlist(aa[i,'pred'])
      prey<-unlist(aa[i,'prey'])
      cof<-aa[i,'rl'][[1]][[1]][['coefficients']]
      pp[1,prey,pred]<-cof[1]
      pp[2,prey,pred]<-cof[2]
      cof<-aa[i,'ru'][[1]][[1]][['coefficients']]
      pp[3,prey,pred]<-cof[1]
      pp[4,prey,pred]<-cof[2]
    }  
    #round(ftable(pp),2)  
    
    predPreySizeParm<-pp
    
    names1<-c( spNames,othspNames)
    names2<-c(spNames,otherFoodName)
    
    k.quantReg<-array2DF(pp) %>% transmute(preyNo=match(Var2,names2),predNo=match(Var3,names1),Var1,Value) %>% 
      pivot_wider(names_from = Var1, values_from = Value) %>% filter(upper_intc >0) 
  } else {
    k.quantReg<-NULL
    doQuant<-FALSE
  }
  
  if (FALSE){
    u<-scan(file.path(stom_data,"pred_prey_size_range_param.csv"),comment.char = "#",sep=',') 
    u<-u[!is.na(u)]
    u<-head(u,-1)
    u<-matrix(u,nrow=(nSpecies)*4,ncol=nSpecies++nOthSpecies,byrow=TRUE)
    round(ftable(u),1)
    
    uu<-array(0,dim=c(4,nSpecies,nSpecies+nOthSpecies),dimnames=list(c('lower_intc','lower_slope','upper_intc','upper_slope'),spNames,c(spNames,othspNames)) )                                                         
    round(ftable(uu),1)
    for (i in (1:4)) {
      ul<-nSpecies*(i-1)+1
      up<-nSpecies*(i) 
      uu[i,,]<-u[ul:up,]
    }
    
    pp<-uu
    round(ftable(pp),1)
  }
#}  ########## end  quantile regression

#k.quantReg %>% rename(pred.no=predNo,prey.no=preyNo)
if ( doQuant) {
  kk<-full_join(k.quantReg %>% rename(pred.no=predNo,prey.no=preyNo),k,by = join_by(prey.no, pred.no)) %>%
      mutate(srMethod=if_else(is.na(srMethod),1L,0L)) %>%
      arrange(pred.no,prey.no)
  #print(kk,n=500)
} else kk<-k

kk<-  kk %>% mutate(pred=if_else(is.na(pred),allSpNames[pred.no],pred), prey=if_else(is.na(prey),allSpNames[prey.no],prey))

#sort(unique(pp2$prey))
#sort(unique(pp2$prey.no))
# the full set of possible combinations of pred, pred age, prey and prey age
pp3<-by(pp2,list(pp2$no),function(x){
  predAge<-1L:info[x$pred.no,'la']
  preyAge<-1L:info[x$prey.no,'la']
  expand.grid(y=1:nYears.model,q=1L:nSeasons,predNo=x$pred.no,predAge=predAge,preyNo=x$prey.no,preyAge=preyAge)
})
pp3<-do.call(rbind,pp3)
dim(pp3)

#filter(pp3,preyAge==1 & preyNo==1 & q==3 & y==1)
#xtabs(~preyNo+preyAge+q,data=pp3)

pp3<-filter(pp3,(predAge > 1) | (q>=recSeason & predAge==1)) %>% 
     filter((    preyAge>1)   | (q>=recSeason &  preyAge==1)) %>% 
  filter(preyNo<otherFoodn | (preyNo==otherFoodn & preyAge==2))  %>%
  arrange(y,q,preyNo,preyAge,predNo,predAge)


# all combinations
b7<-left_join(b,lwdf%>%transmute(s=sp,l_a=a,l_b=b),by = join_by(s))

aa2<-left_join(pp3%>% transmute(y,q,predNo,predAge,preyNo,preyAge),b7%>%  transmute(y,q,predNo=s,predAge=a,predW=WSEA,predL=mean_l,predWl=l_a*predL**l_b,N),by = join_by(y,q, predNo, predAge))
aa2<-left_join(aa2,b7%>%  transmute(y,q,preyNo=s,preyAge=a,preyW=WSEA, preyL=mean_l,preyWl=l_a*preyL**l_b),by = join_by(y,q, preyNo, preyAge))

#filter(aa2,y==1 & q==1 & predNo==1 & predAge==4)


aa2<-aa2 %>% filter(preyNo !=otherFoodn)
# sort(unique(aa2$preyNo))
#xtabs(~preyNo+preyAge+q,data=aa2)

aa2<- aa2 %>% filter(predNo<=nSpecies | predNo>nSpecies & N>0)  %>% mutate(N=NULL)


#aa2<- aa2%>%mutate(preyW=if_else(preyNo==otherFoodn,1,preyW))

#filter(aa2,preyNo==1)
#aa2<-filter(aa2,predW>0 & preyW>0 & preyNo <=nSpecies+1)
#filter(aa2,preyNo==13)

# sort(unique(paste(formatC(k$pred.no,w=2),k$pred)))
# actual comb
aa3<-left_join(aa2,kk,by = join_by(predNo==pred.no, preyNo==prey.no)) %>% as_tibble()
#  sort(unique(paste(formatC(aa3$predNo,w=2),aa3$pred)))
#filter(aa3,preyNo==13)
#aa3<-aa3 %>% mutate(predW=predWl,preyW=preyWl)
aaa<-aa3 %>%  mutate(logRatio=log(predW/preyW),predW=log(predW),preyW=log(preyW),predNo=as.integer(predNo),preyNo=as.integer(preyNo)) 

a4<-filter(aaa,srMethod==0 ) %>% 
       mutate(incl=if_else(logRatio >=min_s & logRatio <=max_s,TRUE,FALSE )) %>%
       mutate(incl=if_else(preyNo==otherFoodn,TRUE,incl)) %>% mutate( min_s=NULL,max_s=NULL,n=NULL,lower_intc=NULL, lower_slope=NULL, upper_intc=NULL, upper_slope=NULL )
if (doQuant) {
  a5<-filter(aaa,srMethod==1) %>% 
    mutate(incl=if_else(logRatio >=lower_intc+ lower_slope*predW  & logRatio <=upper_intc+ upper_slope*predW  ,TRUE,FALSE )) %>%
    mutate(incl=if_else(preyNo==otherFoodn,TRUE,incl)) %>% mutate( min_s=NULL,max_s=NULL,n=NULL,lower_intc=NULL, lower_slope=NULL, upper_intc=NULL, upper_slope=NULL)
    aaa<-rbind(a4,a5) %>% filter(incl) %>%  mutate(incl=NULL,pred=NULL,prey=NULL,srMethod=NULL) %>% arrange(y,q,preyNo,preyAge,predNo,predAge)
} else aaa<-a4 %>% filter(incl) %>%  mutate(incl=NULL,pred=NULL,prey=NULL,srMethod=NULL) %>% arrange(y,q,preyNo,preyAge,predNo,predAge)
dim(aaa)
# filter(aaa,predNo==1 & y==1 &q==1 & predAge==4)
#  
# ftable(xtabs(~q+preyNo+preyAge,data=aaa))
# ftable(xtabs(~predAge+preyNo+preyAge,data=filter(aaa,predNo==1 & q==3)))
# ftable(xtabs(~preyNo+preyAge,data=filter(aaa,predNo==1 & q==1 & y==1 & predAge==4)))


aaa<-aaa %>% select(y,q,predNo,predAge,predW,preyNo,preyAge,preyW,logRatio,preyNo,preyAge,preyW,logRatio)
aaa<-rbind(aaa, aaa %>% transmute(y,q,predNo,predAge,predW) %>% unique() %>% mutate(preyNo=otherFoodn,preyAge=2L,preyW=0,logRatio=0))

#add vulnerability param index
aaa<-left_join(aaa, pp2 %>% transmute(area=1L,predNo=pred.no,preyNo=prey.no,vulneraIdx=no),by = join_by(predNo, preyNo))

#print(filter(aaa,predNo==1 & y==1 &q==1 & predAge==4),n=200)
# sort(unique(paste(formatC(aaa$predNo,w=2),aaa$pred)))

#filter(aaa,preyNo==otherFoodn)
#filter(aaa,predNo==1 & predAge==2 &y==1 & q==3)

a4<-aaa %>% select(y,q,predNo,predAge,predW,preyNo,preyAge,preyW,logRatio, vulneraIdx,preyNo,preyAge,preyW,logRatio) %>% mutate(area=1L)


#xtabs(~predNo+preyNo,data=a4)

if (FALSE) {  #test for data structure
  yqa=1; yqapa=3
  
  tst<-filter(a4,predNo==1,predAge==4  & q==1 & y==1);tst;dim(tst)
  # remember the number of rows
  
  system.time(suitIdx<-a4 %>% nest(data=c(predNo,predAge,predW,preyNo,preyAge,preyW,logRatio, vulneraIdx)) %>% 
    rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,preyW,logRatio,vulneraIdx)))))
  object.size(suitIdx)/1E6
  str(suitIdx,1)
  tst<-suitIdx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge);tst;dim(tst)
  # too many rows. error ?!
  
  system.time(suitIdx<-a4 %>% group_by(y,q,area)%>% nest()  %>% 
    rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,preyW,logRatio,vulneraIdx)))))
  object.size(suitIdx)/1E6
  str(suitIdx,1)
  tst<-suitIdx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge);tst;dim(tst)
  # too many rows. error ?!
  
  
  system.time(suitIdx<-a4 %>% nest_by(y,q,area)  %>% 
               mutate(data=list(data %>% nest(data=c(preyNo,preyAge,preyW,logRatio,vulneraIdx)))))
  object.size(suitIdx)/1E6
  str(suitIdx,1)
  tst<-suitIdx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge);tst;dim(tst)
  # too many rows. error ?!
  
  
  system.time(suitIdx<-a4 %>% nest_by(y,q,area)  %>% 
                mutate(data=list(data %>% nest_by(predNo,predAge,predW))))
  object.size(suitIdx)/1E6
  str(suitIdx,1)
  tst<-suitIdx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge);tst;dim(tst)
  # no rows as there should be
  # takes longer time to produce, but is correct!
  
  x1<- suitIdx[yqa,]
  y<-x1$y
  q<-x1$q
  
  x2<-x1$data[[1]][yqapa,]
  predNo<-x2$predNo
  predAge<-x2$predAge
  predW<-x2$predW
  
  v1<-x2$data[[1]] %>% arrange(preyNo,preyAge)
  v2<-suitIdx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge)
  v3<-filter(a4,predNo==1,predAge==4  & q==1 & y==1) %>% arrange(preyNo,preyAge)
  dim(v1);dim(v2);dim(v3)
}  

suitIdx<-a4 %>% nest_by(y,q,area) %>% mutate(data=list(data %>% nest_by(predNo,predAge,predW)))


if (FALSE) {  # compare with the old SMS data
  selPred<-c(1:100)
  s<-unnest(suitIdx,cols = c(data))
  s
  s$data[[1]]
  s<-filter(s,predNo %in% selPred) %>% unnest(cols = c(data)) %>% filter(preyNo !=otherFoodn)  %>% 
    mutate(prey=factor(spNames[preyNo]))
  s 

  summary(s)
  by(s,s$predNo,function(x)ftable(xtabs(~prey+preyAge+q,data=x)))
  
  
  #### compare data for calc of M2 from old SMS and rsms
  spl<-allSpNamesLong
  sp<-1L:length(spl)
  names(sp)<-spl
  
  sps<-allSpNames
  sps<-1L:length(sps)
  names(sps)<-allSpNames
  
  old<-read_csv(file="who_eats_whom_level1.csv")
  head(old)
  old<-old %>% transmute(y=Year+off.year,Year,q=Quarter,Predator,predNo=sp[Predator],predAge=Predator.age+off.age,
                         Prey,preyNo=sp[Prey],preyAge=Prey.age+off.age,M2SMS=Part.M2) %>% filter(predNo %in% selPred)
  
  for (p in sort(unique(s$predNo))) {
    cat(p,allSpNames[p],"\nold:",dim(filter(old,predNo==p))[[1]],"\nNew:",dim(filter(s,predNo==p))[[1]],'\n\n') 
  }

  
  
  ab<-full_join(s,old,by = join_by(y, q, predNo, predAge,preyNo, preyAge)) %>% mutate(preyW=NULL,predW=NULL,vulneraIdx=NULL) %>% 
    arrange(y,q,predNo,predAge,preyNo,preyAge)
  ab
  dim(ab)
  
  ab %>% filter(!is.na(prey) & !is.na(Prey))  # both
  ab %>% filter(is.na(prey) & !is.na(Prey))  #old SMS only
  ab %>% filter(!is.na(prey) & is.na(Prey))  #rsms  only
  ###
  
  oldPreds<-c(tail(predNames,-5),head(predNames,5))
  minpp<-scan(file=file.path(data.path,"min_max_size_pref.out"),comment.char="#",n=nSpecies*length(preds))
  minpp<-matrix(minpp,nrow=length(preds),ncol=nSpecies,byrow = TRUE,
                dimnames=list(oldPreds,spNames))
  minpp<-array2DF(minpp) %>% transmute(pred=Var1,prey=Var2,minpp=Value)
  
  maxpp<-scan(file=file.path(data.path,"min_max_size_pref.out"),comment.char="#",skip=24)
  maxpp<-maxpp[1:(nSpecies*length(preds))]
  maxpp<-matrix(maxpp,nrow=length(preds),ncol=nSpecies,byrow = TRUE,
                dimnames=list(oldPreds,spNames))
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
}
if (FALSE) {
    partM2Idx<-filter(a4,preyNo !=otherFoodn) %>% mutate(partM2=0,suit=0, availFood=0,predW=NULL,preyW=NULL,logRatio=NULL, vulneraIdx=NULL) %>%
    nest(data=c(predNo,predAge,preyNo,preyAge,partM2,suit,availFood)) %>% 
    rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,partM2,suit,availFood))))
  object.size(partM2Idx)/1E6
}


partM2Idx<-filter(a4,preyNo !=otherFoodn) %>% mutate(partM2=0,suit=0, availFood=0,predW=NULL,preyW=NULL,logRatio=NULL, vulneraIdx=NULL) %>%
  nest_by(y,q,area) %>% mutate(data=list(data %>% nest_by(predNo,predAge)))
object.size(partM2Idx)/1E6

yqa=1; yqapa=3
partM2Idx[yqa,]$data[[1]][yqapa,]$data[[1]] %>% arrange(preyNo,preyAge)
#partM2Idx
#partM2Idx$data[[1]]
#partM2Idx$data[[1]]$data[[1]]


# suitIdx
# suitIdx[,'y']
# suitIdx[1,'data'][[1]][[1]]
# #suitIdx[1,'data'][[1]][[1]][['data']]
# suitIdx[1,'data'][[1]][[1]][['data']][[2]]



## stomach contents

oth<-ss$prey==otherFoodName

ss[oth,'prey.no']<-otherFoodn
ss[oth,'prey.size.class']<-1
ss[oth,'prey.size.w']<-1
ss[oth,'pred.prey.size.w']<-1 # reset later on

# glimpse(filter(ss,oth))
aBasis<-ss%>% transmute(area=as.integer(SMS_area),
                        year=as.integer(year),
                        y=as.integer(year+off.year),
                        q=as.integer(quarter+off.season),
                        predC=pred,
                        pred=pred.no,
                        predSize=pred.size,
                        predSizeClass=as.integer(pred.size.class),
                        predMeanLength=as.integer(pred.mean.length),  
                        predSizeW=pred.size.w,
                        noSampl=as.integer(noSampl),
                        phi=phi,
                        preyC=prey,
                        prey=as.integer(prey.no),
                        preySize=prey.size,
                        preySizeClass=as.integer(prey.size.class),
                        preyMeanLength=as.integer(prey.mean.length),
                        preySizeW=prey.size.w,
                        logPPsize=log(pred.prey.size.w),
                        type=type,
                        stomcon=stomcon )
aBasis[aBasis$preyC=='OTH','logPPsize']<-NA

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize) %>% select(-preySizeClass,-preyMeanLength,-predMeanLength)
# sum up to 1 (into proportions)
aBasis<-aBasis %>% group_by(area,year,y,q,predC,pred,predSize,predSizeClass,predMeanLength,predSizeW,noSampl,
                       phi,preyC,prey,preySize,preySizeClass,preyMeanLength,preySizeW,logPPsize,type) %>% summarize(stomcon=sum(stomcon),.groups = "drop") %>%  
  group_by(area,year,y,q,predC,pred,predSize,predSizeClass,predMeanLength,predSizeW,noSampl,phi) %>% mutate(stomcon=stomcon/sum(stomcon)) %>% ungroup()
   
# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize) %>% select(-preySizeClass,-preyMeanLength,-predMeanLength)

nstom<-aBasis %>% group_by(year,q,pred, predC, predSizeClass, predSize) %>% summarize(nstom=min(noSampl),.groups='drop') %>% arrange(nstom,pred,predSizeClass) 
write.csv(nstom,file=file.path(data.path,'no_stom.csv'),row.names = FALSE)

delStom<-read_csv((file=file.path(data.path,'stom_input','stomach_delete.csv')),show_col_types = FALSE) %>% filter(nstom<0) %>% mutate(nstom=NULL,delete=TRUE)
aBasis<-left_join(aBasis,delStom,by = join_by(year, q, predC, pred, predSize, predSizeClass)) %>% 
         filter(!delete | is.na(delete)) %>% mutate(delete=NULL)


minNoSampl<-aBasis %>% group_by(pred) %>% summarize(minNoSampl=min(noSampl),.groups='drop')


#vulneraIdx.df
aBasis<-left_join(aBasis,select(vulneraIdx.df,pred,prey,vulneraIdx=Value),by=join_by(predC==pred,preyC==prey))

aBasis<-aBasis %>% group_by(area,year,y,q, pred,  predSizeClass,prey) %>% 
  mutate(stomconTot=sum(stomcon),minL=min(preySizeClass),nPreyGroups=dplyr::n(),firstL=if_else(minL==preySizeClass,TRUE,FALSE)) %>%
  ungroup() %>% arrange(area,year,y,q, pred,  predSizeClass,prey,preySizeClass) %>%
  group_by(area,year,y,q, pred,  predSizeClass) %>% mutate(nPreyGroups=sum(firstL),preyIdx=cumsum(firstL), minL=NULL) %>% ungroup()

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize)%>% select(-preySizeClass,-preyMeanLength,-predMeanLength)


#tmp<-filter(aBasis2,y==8&q==1&pred==1 & predSizeClass >5) %>% select(predSizeClass,prey,preySizeClass,stomcon,stomconTot,nPreyGroups,firstL,preyIdx) %>% arrange (predSizeClass , prey ,preySizeClass )
#print(tmp,n=20)


logStomObsVar<-1L:length(predNames)
names(logStomObsVar)<-predNames 
aBasis$stomObsVarIdx<-logStomObsVar[aBasis$predC]


logStomObsVar[]<-log(0.5^2)

#alfa0<-wDiri*stomVar - 1.0  # Dirichlet alfa0 from number of samples
minNoSamplV<-minNoSampl$minNoSampl
minNoSamplV<-log(1/minNoSamplV+0.1)
idx<-info[, 'stomachVariance'][predNames] %in% c(3)
logStomObsVar[idx]<-minNoSamplV[idx]

names(logStomObsVar)<-NULL

usestomObsVar<- !info[, 'stomachVariance'][predNames] %in% c(4)
#glimpse(aBasis)

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize)%>% select(-preySizeClass,-preyMeanLength,-predMeanLength)
if (FALSE) {
  newS <-aBasis %>% mutate(preyIdx=NULL,stomObsVarIdx=NULL,vulneraIdx=NULL) %>%   
          filter(preyC !='OTH') %>% arrange(y,q,pred,predSizeClass,prey,preySizeClass) 
  names(newS)
  glimpse(newS)
  s1<-read_table("summary_stom.out") %>% filter(Prey.no>0) %>% 
        mutate(SMS.area=NULL,Size.model=NULL,Predator.available.food=NULL,Diri.like=NULL,multinom.like=NULL,Obj.N.contrib=NULL,Prey.number.prob=NULL,Prey.avail=NULL,Prey.avail.part=NULL,avil.Prey.number.sum=NULL,Residual=NULL) %>%
       mutate(N.haul=NULL,N.haulstomcon.hat=NULL,avail.prey.number.sum=NULL,Stom.var=NULL,obj.contrib=NULL,Stom.used.avail=NULL,stom.used.like=NULL,L.N.bar=NULL,)
  glimpse(filter(s1,Predator.no==16))
   glimpse(newS)

  
  oldS<-s1 %>%transmute(year=Year,y=as.integer(Year+off.year),q=as.integer(Quarter.no),predC=oldPredNames[Predator.no],pred=oldNewSp_n[Predator.no],predSizeOld=Predator.length,
                        predSizeClass=Predator.length.class,predMeanLengthOld=Predator.length.mean,predSizeWOld=Predator.weight,
                        noSamplOld=N.samples,phiOld=Diri.sum.p,
                        preyC=oldAllSpeciesNames[Prey.no],
              prey=as.integer(Prey.no-15),preySizeClass=as.integer(Prey.length.class),Prey.weight,stomconSMS=stom_data)  %>% 
             arrange(y,q,pred,predSizeClass,prey,preySizeClass) 
  glimpse(oldS)
  glimpse(newS)
  setdiff(names(newS),names(oldS))

  aa<-full_join(oldS,newS,by = join_by(year, y, q, predC, pred, predSizeClass, preyC, prey, preySizeClass))  %>% 
       mutate(area=NULL,firstL=NULL,nPreyGroups =NULL,stomconTot=NULL)
  glimpse(aa)
  names(aa)
  on<-aa %>% transmute(year,y,q,pred,predC,predSizeClass ,predSize,predMeanLength,predMeanLengthOld,predSizeW,predSizeWOld,noSampl,noSamplOld,phi,phiOld,
                   preyC,prey,preySizeClass,preySize,preySizeW,Prey.weight,stomcon,stomconSMS) %>% arrange(y,q,pred,predSizeClass,prey,preySizeClass) 
  
  oldNew<-function(x,old,new,lim=0.01,ignore=1E-4) {
   x$dif<-abs((old-new)/old)
   x$old<-old;x$new=new
   x<-x %>% transmute(year,y,q,pred,predC,predSizeClass,predSize,prey,preyC,preySizeClass,preySize,old,new,dif) %>% filter(dif>lim & old>ignore) %>% arrange(desc(dif)) %>%
     mutate_if(is.character,as.factor)
   x
  }

 a<- oldNew(x=on,old=on$stomconSMS,new=on$stomcon,lim=0.01,ignore=1E-5)
 a
 summary(a)
 xtabs(~predC+preyC,data=a)
  
 a<-oldNew(x=on,old=on$noSamplOld,new=on$noSampl,lim=0.01)
 a
  summary(a)
 xtabs(~predC+year,data=a)
 

  a<-oldNew(x=on,old=on$phiOld,new=on$phi,lim=0.01)
  a
  summary(a)
  xtabs(~predC+year,data=a)
  
  
  oldNew(x=on,old=on$predMeanLengthOld,new=on$predMeanLength,lim=0.01)
  oldNew(x=on,old=on$Prey.weight,new=on$preySizeW,lim=0.001)
}  
  
if (FALSE) {
  chk<- aBasis %>% select(area, y,q,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
    nest(data=c(pred,predSizeClass,predSizeW,stomObsVarIdx,noSampl,   phi,nPreyGroups,prey,preySizeClass,preySizeW,logPPsize,type,stomcon, stomconTot ,preyIdx,firstL,vulneraIdx))
  chk[1,'data'][[1]]
  
  chk2<-  chk%>%rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preySizeClass,preySizeW,logPPsize,type,firstL,preyIdx,stomconTot ,stomcon,vulneraIdx))))
  chk2[1,'data'][[1]]
  chk2[1,'data'][[1]][[1]]$data
}


###### stom ALK
s<-read_delim(alk_data,show_col_types = FALSE)
# sort(unique(s$prey))
d2010<-filter(s,year==2013) %>%mutate(year=2010)  #snyd
s<- rbind(filter(s,year !=2010),d2010)
s<-s%>%mutate(prey.no=as.integer(prey.no),prey.age=as.integer(prey.age),quarter=as.integer(quarter),year=as.integer(year),prey.size.class=as.integer(prey.size.class))


# new number order
s$prey.no<-oldNewSp_n[s$prey.no]
s<-filter(s,prey %in% union(preyNames,setdiff(predNames,othspNames)))  #just preys and analytical predators


alkBasis<-s%>% transmute(area=as.integer(SMS_area),
                        year,
                        y=as.integer(year+off.year),
                        q=as.integer(quarter+off.season),
                        sp=prey,
                        s=prey.no,
                        age=as.integer(prey.age),
                        a=as.integer(prey.age+off.age),
                        size=prey.size,
                        sizeClass=as.integer(prey.size.class),
                        meanLength=ALK.MeanL,  
                        alk=ALK/100
)

#a        <-by(alkBasis,list(alkBasis$y,alkBasis$q),function(x) {y<-tapply(x$alk,list(x$s,x$a,x$sizeClass),sum); y[is.na(y)]<-0; y})
#str(a)
#round(ftable(a[[1]]),3)

alkBasis<-left_join(alkBasis,data.frame(s=info[,'s'],la=info[,'la']),by = join_by(s))
plusG<-alkBasis$a>alkBasis$la
# summary(plusG)
alkBasis$wf<-1
alkBasis[plusG,'wf']<-exp(-(alkBasis[plusG,'a']-alkBasis[plusG,'la'])*0.5) #quick and dirty plusgroup
alkBasis[plusG,'a']<-alkBasis[plusG,'la']
tst1<-filter(alkBasis,s==2&a>6)
#tst1


alkBasis<-alkBasis %>% group_by(area,year,y ,q, sp,s, a ,size,sizeClass) %>% 
  summarize(meanLength=weighted.mean(x=meanLength, w=wf),alk=weighted.mean(x=alk, w=wf),.groups = "drop") %>%
  group_by(area,year,y ,q, sp,s, a ) %>%
  mutate(alk=alk/sum(alk))  %>% ungroup()
#tst1
#filter(alkBasis,s==2&a>6)
  
fl<- alkBasis %>% group_by(area,y,q,s,a) %>% summarize(minSize=min(sizeClass), maxSize=max(sizeClass),.groups = "drop")

#tapply(fl$maxSize,list(fl$s,fl$y),max)
#tapply(fl$minSize,list(fl$s,fl$y),min)
maxSizeCl<-max(fl$maxSize)
alk<- left_join(alkBasis %>% select(area,y,q,s,a,sizeClass,meanLength,alk),fl,by = join_by(area, y, q, s, a)) %>% mutate(s=as.integer(s))
if (FALSE) {
  chk<-filter(alk,y==40,q==3,s==2)
  chk<-xtabs(alk~a+sizeClass,data=chk)
  round(cbind(chk,rowSums(chk)),3)
}
#alk$alkObsNo<-1L:dim(alk)[[1]]
#aBasis$stomObsNo=1L:dim(aBsis)[[1]]

#check that ALK and STOM match
alk <- alk %>% mutate(yqALK=paste(y,q,sep='-'))

alkyq<-alk %>% select(area,y,q,yqALK) %>% unique() %>% arrange(area,y,q) 
aBasis<-left_join(aBasis,alkyq,by = join_by(area, y, q))
#glimpse(aBasis)


a1<-sort(unique(alk$yqALK))
a2<-sort(unique(aBasis$yqALK))
stopifnot(length(setdiff(a2,a1))==0) # missing ALKs
a1a2<-setdiff(a1,a2) #There are ALK for y and q combinations that are not used
#a1a2
if (length(a1a2)>0) {
  alk<-alk %>% filter(yqALK %in% a2)
}


# check there is at least on ALK observation for a predator size class
stomChkPred<-aBasis %>% filter(pred<=nSpecies) %>%select("area","year","y","q","predC","pred","predSizeClass") %>% unique() %>% mutate(fromStom=TRUE)
alkChkPred<-alk %>%  select("area","y","q","s","sizeClass") %>% unique() %>% mutate(fromALK=TRUE)
chkPred<-full_join(x=stomChkPred,y=alkChkPred,by=join_by(area==area,y==y,q==q,pred==s, predSizeClass==sizeClass))
# summary(stomChkPred);summary(alkChkPred); summary(chkPred)
#print(filter(chkPred,is.na(fromStom)),n=20) #  ALK but no STOM, not a problem
problem<-filter(chkPred,is.na(fromALK))
print(problem,n=20) #  STOM but no ALK, that is a problem

aBasis<-anti_join(aBasis,problem%>%mutate(fromStom,fromALK)) 
#stopifnot(dim(problem)[[1]]==0)

# check there is at least on ALK observation for a prey size class foud in the stomachs
stomChkPrey<-aBasis %>% select("area","year","y","q","preyC","prey","preySizeClass") %>% filter(preyC !=otherFoodName) %>%unique()%>% mutate(fromStom=TRUE)
alkChkPrey<-alk %>%  select("area","y","q","s","sizeClass") %>% unique() %>% mutate(fromALK=TRUE)
chkPrey<-full_join(x=stomChkPrey,y=alkChkPrey,by=join_by(area==area,y==y,q==q,prey==s, preySizeClass==sizeClass))
problem<-filter(chkPrey,is.na(fromALK))
#summary(stomChkPrey);summary(alkChkPrey);summary(chkPrey)
print(problem,n=200) #  STOM but no ALK, that is a problem
#dim(problem)
stopifnot(dim(problem)[[1]]==0)

alkyq<-alk %>% select(area,y,q,yqALK) %>% unique() %>% arrange(area,y,q) %>% mutate(yqALK=1L:dplyr::n())
aBasis$yqALK<-NULL
aBasis<-left_join(aBasis,alkyq,by = join_by(area, y, q))
#summary(aBasis)
filter(aBasis,is.na(yqALK) & pred <=nSpecies ) %>% select(area,y,year,q,pred,predC) %>% unique()

if (FALSE) {
  stom<- aBasis %>% select(area, y,q,yqALK,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
    nest(data=c(pred,predSizeClass,predSizeW,stomObsVarIdx,noSampl,   phi,nPreyGroups,prey,preySizeClass,preySizeW,logPPsize,type,stomcon, stomconTot ,preyIdx,firstL,vulneraIdx)) %>% 
    rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preySizeClass,preySizeW,logPPsize,type,firstL,preyIdx,stomconTot ,stomcon,vulneraIdx))))
  
  object.size(stom)/1E6
}

stom<- aBasis %>% select(area, y,q,yqALK,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
  nest_by(area,y,q,yqALK) %>% mutate(data=list(data %>% nest_by(pred, predSizeClass, predSizeW, stomObsVarIdx, noSampl, phi, nPreyGroups)))  



if (FALSE){
  str(stom$data[[1]],2)
  str(stom$data[[1]][['data']],2)
  stom[,'y']
  stom[1,'data'][[1]][[1]]
  stom[1,'data'][[1]][[1]][['data']]
  stom[1,'data'][[1]][[1]][['data']][[1]]
  stom[1,'data'][[1]][[1]][['data']][[2]]
  stom[1,'data'][[1]][[1]]$data[[3]]
  

  stst<-unnest(stom,cols = c(data))
  stst<-unnest(stst,cols = c(data))
  names(stst)
  tst<-select(stst,area,y,q,pred,predSizeClass,prey,preySizeClass)
  dim(stst);dim(tst); dim(unique(tst))
  
}

if (FALSE){
  alk2<-alk%>%  nest(data=c(a,s,minSize,maxSize,sizeClass,meanLength,alk)) %>% 
    rowwise() %>% mutate(data=list(data %>% nest(data=c(sizeClass,meanLength,alk))))
  object.size(alk2)/1E6
}


alk<-alk %>% nest_by(area,y,q,yqALK)  %>% mutate(data=list(data %>% nest_by( a,s, minSize, maxSize )))
#object.size(alk3)/1E6

#str(alk$data[[1]],2)
#alk[1,'data'][[1]][[1]]
#alk[1,'data'][[1]][[1]][['data']]




#other food
of<-scan(file.path(dir,"other_food.in"),comment.char = "#",quiet=TRUE)
of<-head(of,-1)
names(of)<-c(spNames,othspNames)
of
} #end if multi

#Un<-matrix(0.0, nrow=sum(info[,"last-age"]-sms@first.age+1L),ncol=nYears )
Un<-matrix(0.0, nrow=sum(nlogN),ncol=nYears.model )

Uf<-matrix(0.0, nrow=length(unlist(sms@keyLogFsta)),ncol=nYears.model)


#return values
rl<-list(
   data=list(    
     silent=TRUE,
     options=sms,
     sms.mode=sms@VPA.mode,                            #
     useNbar=sms@use.Nbar,
     M2iterations=sms@M2.iterations,
     info=info,                                        # Various information for each species with analytically assessment
     nSpecies=nSpecies,                                # Number of species with analytically assessment
     preds=predators,                                  # predator species
     preys=preys,                                      # prey species
     isPrey=isPrey,
     nOthSpecies=nOthSpecies,                          # Number of other predators (used for multispecies mode) 
     nSpeciesAll=nAllSp,
     otherFoodn=otherFoodn,                            # Species number for other food
     nYears.model=nYears.model,                        # Number of years used in the model
     extendModel=extendModel,
     nAges=nAges,
     nSeasons=nSeasons,                                 # number of seasons
     fq=fq,
     lq=lq, 
     fa=fa,
     fqa=fqa,                                          # fist quarter for each age
     faq=faq,                                          # first age for each quarter
     yqIdx=yq_idx,                                     # index for year-season combinations
     minAge=as.integer(sms@first.age),                 # A vector of integers defining the the lowest age class in the assessment for each species.
     #maxAge=as.integer(sms@max.age.all),              # Maximum age for all species. The actual age range by species is given in "info"
     recAge=as.integer(sms@first.age),                 # recruitment age
     years.model=sms@first.year.model:(sms@last.year.model+extendModel),   # A vector of the years used in the model
     spNames=spNames,                                  # Species names of species with analytically assessment
     othspNames=othspNames,                            # Species names of "other" species (without analytically assessment)
     allSpNames=allSpNames, 
     allSpNamesLong=allSpNamesLong,
     predNames=predNames,
     preyNames=preyNames,
     otherFoodName=otherFoodName,
     fleetNames=fleetNames,                            # names of survey fleets
     spawnSeason=1L,
     recSeason=recSeason,                              # recruitment season
     fbarRange=fbarRange,                              # Minimum and maximum age by species used to calculate Fbar 
     off.age=off.age,                                  # offset between age and the age index used in RSMS (where fist age index  have to be 1 )
     off.year=off.year,                                # offset between year and the year index used in RSMS (where fist year index  have to be 1 )
     #off.season=off.season,
     #off.species=off.species,
     #off.oths=off.oths,                          
     seasonalCatches=seasonalCatches,                  # are seasonal catches used in the catch observation likelihood
     catchSepYear=sms@catch.sep.year,
     FfromSeparableModel=FfromSeparableModel,         # are seasonal F derived from seasonal catch observations
     fixMapseasonalF=fixMapseasonalF,
     bfix=bfix,
     useFrandomWalk=useFrandomWalk,
     idxLogYearEffectF=idxLogYearEffectF,
     doProcessNoise=inclProcess,
     doProcessN_any=inclProcess %in% c(1,2),
     doProcessN_old=inclProcess %in% c(1),
     doProcessN_none=inclProcess==3,
     useRho=use_rho,                                   # use correlation in F at age
     stockRecruitmentModelCode=stockRecruitmentModelCode,
     stockRecruitmentModelCodeAdd=stockRecruitmentModelCodeAdd,
     zeroCatchYearExists=zeroCatchYearExists,
     zeroCatchYearExistsSp=zeroCatchYearExistsSp,
     zeroCatchYear=zeroCatchYear,
     zeroCatchSeasonAge=zeroCatchSeasonAge,
     keyLogFsta=keyLogFsta,                         # A matrix of integers. The number of rows is equal to the number of species fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of the fishing mortality states. '-1' is used for entries where no fishing mortality applies. For the valid entries consecutive integers starting 1 must be used, because they are used as indices in the corresponding state vector. If the same number is used for two fleet-age combinations, then the fishing mortality for those are assumed equal (linked to the same state).
     nlogF=nlogF,                                   # A vector with the sum of number of F state age groups by species. Each species number corresponds to the  number of unique (and not -1 value) in each (species) row of keyLogFsta
     keyLogFstaSd=keyLogFstaSd,                     # A matrix which links each F state age group by species to the logSdLogFsta=
     nlogFfromTo=nlogFfromTo,
     nlogN=nlogN,                                   # A vector with the sum of number of N state age groups by species
     nlogNfromTo=nlogNfromTo,
     logNfirstYparamfromTo=logNfirstYparamfromTo,
     logNrecruitParamfromTo=logNrecruitParamfromTo,
     inclSsbR=inclSsbR,
     SsbRweight=SsbRweight,
     keyVarLogN=keyVarLogN,
     keyVarObsCatch.df=keyVarObsCatch.df,
     minVarObsCatch=sms@min.catch.CV,
     minVarObsSurvey=minCVsurvey,
     keyVarObsSurvey=keyVarObsSurvey,
     keyCatchability=keyCatchability,
     propMat=  propMat,      
     stockMeanWeight=stockMeanWeight,
     catchMeanWeight=catchMeanWeight,
     catchNumber=catchNumber,
     seasFprop=seasFprop,
     logSeasFprop=logSeasFprop,
    # ageQCatchMis=ageQCatchMis,
     natMor= natMor,
     propF=zero,
     propM=zero,
     keyCatch=keyCatch,
     keylogSeasonF=keylogSeasonF,    
#     catchNoSeason=catchNoSeason,
 #    catchNoSeasonUsed=catchNoSeasonUsed,
     logCatchObs=logCatchObs,
     recruitYears=rec_years,
     #nCatchObs=length(logCatchObs),
     keylogSeparAgeF=keylogSeparAgeF,
     keySurvey=keySurvey,
     keySurvey.overview=keySurvey.overview,  
     logSurveyObs=logSurveyObs
     #nSurveyObs=length(logSurveyObs)
   ),
   
  parameters=list(
      logSdLogObsCatch=logSdLogObsCatch,
      logFSeasonal=logFSeasonal,
      logCatchability=logCatchability,
      logSdLogObsSurvey=logSdLogObsSurvey,
      logQpow=logQpow,
      logSeparAgeF=logSeparAgeF,
      logYearEffectF=logYearEffectF,
      logTechCreep=logTechCreep,
      logSdLogFsta=logSdLogFsta,
      logSdLogN=logSdLogN,
      logNrecruitParam=logNrecruitParam,
      logNfirstYparam=logNfirstYparam,
      logSsbRsd=logSsbRsd,
      rho=rho,
      rec_loga=rec_loga, 
      rec_logb=rec_logb,
       Un=Un,
      Uf=Uf
    ))

if (multi){
  rl[['data']]<-c(rl[['data']],
     list(
      #multi species
      consum=consum,
      meanL=meanL,
      logPropM2=logPropM2,          
      natMor1=natMor1,         
      otherN=otherN, 
      stom=stom,
      alk=alk,
      maxSizeCl=maxSizeCl,
      otherFood=of,
     # predPreySizeParm=predPreySizeParm,
      vulneraIdx=vulneraIdx,
      suitIdx=suitIdx,
      partM2Idx=partM2Idx,
      overlap=overlap, 
      overlapIdx=overlapIdx,
      usestomObsVar=usestomObsVar)
  )
      
  rl[['parameters']]<-c(rl[['parameters']],
   list(  
    vulnera=vulnera,
    overlapP=overlapP,
    logStomObsVar=logStomObsVar
   ))
}
return(rl)
}


