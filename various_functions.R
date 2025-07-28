cleanrun<-function(silent=TRUE) {
  dels<-c('obj','opt','sdrep','myRep','myMap')
  for (i in dels) {
    varName<-i
    if (!silent) cat(varName,exists(varName),'\n')
    if (exists(varName)) rm(list=as.character(substitute(varName)),pos=1) 
  }
}

cleanup<-function(){for(i in dev.list()) if (names(dev.off())=='null device') break()}

likelihoodCompare<-function(inpRdata,labels,doPrint=TRUE) {
  ll<-lapply( 1:length(inpRdata),function(i) {
    load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
    a<-sms$rep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a)) ;return(round(a,1))
  })
  if (doPrint) lapply(ll,print)
  invisible(ll)
}
  

getParameters<-function(runName,printIt=TRUE,details=2){
 load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=FALSE)
  cat('Hesssian:',sms$sdrep$pdHess,'\n')
  a<-extractParameters(sdrep=sms$sdrep, myMap=sms$map,data=sms$data,parameters=sms$parameters,lu=sms$lu)[[details]]
  if (printIt) print(a,n=500)
  invisible(a)
}

writeSeasonalF<-function(inp,outfile="proportion_of_annual_f.in",dir=data.path) {
  load(file=file.path(data.path,paste0(inp,".Rdata")),verbose=T)
  x<- sms[['rep']]$res %>% select(s,year,quarter,age,FisQ) %>% 
    group_by(s,year,age) %>% 
    mutate(FiProp=FisQ/sum(FisQ)) %>% mutate(FiProp=if_else(is.na(FiProp),0,FiProp),FisQ=NULL)
  summary(x)
  xx<-tapply(x$FiProp,list(x$s,x$year,x$quarter,x$age),sum)  
  xx[is.na(xx)]<-0
  xx[1,1,,]
  d<-dim(xx)
  dn<-dimnames(xx)
  ofile<-file.path(dir,outfile)
  spNames<-sms[['data']]$spNames
  cat('### proportion of annual F\n',file=ofile)
  for (s in seq_len(d[1])) {
    cat("##  ",spNames[s], "## \n",file=ofile,append=TRUE)
    for (y in  dn[[2]]) {
      cat("#",spNames[s]," year:",y,'\n',file=ofile,append=TRUE)
      p<-xx[s,y,,]
      if (round(p[3,2],2) !=0.25) p<-round(p,4)
      write.table(p,col.names=FALSE,row.names = FALSE,append=TRUE,quote=FALSE,file=ofile)
    }
  }
  cat('-999  # checksum\n',file=ofile,append=TRUE)
}  

#writeSeasonalF(inp=runName)


logLik.RSMS<-function(sms, ...){
  ret<- -sms$opt$objective
  attr(ret,"df")<-length(sms$opt$par)
  attr(ret,"nall")<- length(sms$data$logCatchObs)+ length(sms$data$logSurveyObs)  #stomach obs are still missing
  attr(ret,"nobs")<- length(sms$data$logCatchObs)+ length(sms$data$logSurveyObs)
  class(ret)<-"logLik"
  ret
}


AIC<-function(ll){
  K<- attributes(ll)$df  # no of parameters
 # n<- attributes(ll)$nobs
  -2*ll[1] + 2*K
}


AICc<-function(ll){
  K<- attributes(ll)$df
  n<- attributes(ll)$nobs
  -2* ll[1] + 2*K*(n /(n-K-1))
}

AICCompare<-function(inpRdata,labels, pval=TRUE) {
  ll<-lapply( 1:length(inpRdata),function(i) {
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      ll<-logLik.RSMS(sms)
      AIC(ll)
      AICc(ll)
      data.frame(Model=labels[i],logL=-ll[1],npar=attributes(ll)$df, AIC=AIC(ll),AICc=AICc(ll))
  })
  out<-do.call(rbind,ll)

  o <- dim(out)[[1]]
  if(dim(out)[[1]]>=2 & pval){
    out <-out %>% transmute(Model,logL,npar,AIC,AICc,D=NA,df=NA,Pval =NA)
    out<-arrange(out,desc(npar))
    for (x1 in (1:(o-1))) {
      x2<-x1+1
      if(out[x1,'npar']>out[x2,'npar']) {
        df <- out[x1,'npar'] - out[x2,'npar']
        D <- 2*(out[x1,'logL']-out[x2,'logL'])
        P <- 1- pchisq(D,df)
        out[x2,'Pval']<- P
        out[x2,'df']<-df
        out[x2,'D']<-D
      }
    }
  }
  out
}


write_proportion_of_annual_f<-function(data) {
  outFile<-'proportion_of_annual_f.dat'
  cat("# proportion af annaual F by season\n",file=outFile)
  for (s in seq_len(data$nSpecies)) {
    cat("################   ", data$spNames[s], "   ################\n",file=outFile,append=TRUE)
    x<-data$seasFprop[[s]]
    cat(paste0("#  ",data$speciesNames[s],"\n"),file=outFile,append=TRUE)
    for (y in seq_len(data$nYear)){
      cat('# Year',y-data$off.year,'\n',file=outFile,append=TRUE)
      write.table(x[y,,],file=outFile,append=TRUE,row.names=FALSE,col.names=FALSE)
    }
  }
  cat(-999, "  # check value",file=outFile,append=TRUE)
}
#write_proportion_of_annual_f(inp_all[['data']])

transExternalData<-function(inp='summary.out',outSet='SMS_old_detailed',
                               spNames='COD',
                               exSpeciesNames=c("FUL","GLT","HEG","KTW","GBG","GNT","PUF","RAZ","RAJ","GUR","WHM","NHM","GSE","HBP","HKE","COD","WHG","HAD","POK","MAC","HER","NSA","SSA","NOP","SPR","PLE","SOL")) {
  a<-read.table(file=inp,header=TRUE) %>%as_tibble() %>%
    transmute(year=Year,quarter=Quarter,age=Age,species=exSpeciesNames[Species.n],s=match(species,spNames),M2,N,FF=`F`) %>%
    filter(!is.na(s) & M2>=0) 
  
  sms<-list(rep=list(res=a),data=list(allSpNames=spNames,allSpNamesLong=spNames))
  save(sms,file=file.path(data.path,paste0(outSet,".Rdata")))
}

# transExternalData(inp='summary.out',outSet='SMS_old_detailed')


transExternalSummary<-function(inp='summary_table_raw.out',outSet='ICES',
                               spNames='COD',
                              exSpeciesNames=c("FUL","GLT","HEG","KTW","GBG","GNT","PUF","RAZ","RAJ","GUR","WHM","NHM","GSE","HBP","HKE","COD","WHG","HAD","POK","MAC","HER","NSA","SSA","NOP","SPR","PLE","SOL")) {
  a<-read.table(file=inp,header=TRUE) %>%as_tibble() %>%
    transmute(year=Year,species=exSpeciesNames[Species.n],SSB=SSB,recruit=Rec,Fbar=mean.F, FbarAnn=mean.F,yield=Yield,s=match(species,spNames)) %>%
    filter(!is.na(s))

 sms<-list(data=list(allSpNames=spNames,allSpNamesLong=spNames),rep=list(resSummary=a))
 save(sms,file=file.path(data.path,paste0(outSet,".Rdata")))
 #load(file=file.path(data.path,paste0(outSet,".Rdata")))
}

# transExternalSummary(inp='summary_table_raw.out',outSet='ICES_single_sp')
                              

calcCV<-function(sdrep,data) {
  cv_par<-function(par,key,type){
    species<-rownames(key)
    rownames(key)<-paste0('var_',1:dim(key)[[1]])
    s<-data.frame(species,Var1=rownames(key),sp.no=1:length(species))
    a<-array2DF(key)  %>%filter(Value>0) %>% 
      mutate(Var=as.numeric(substr(Var1,5,6)),age=as.numeric(substr(Var2,5,6)),param=par) %>%
      arrange(Var,age)  %>%group_by(Var) %>% mutate(dup=duplicated(Value)) %>%filter(!dup) %>%
      mutate(nexta=lead(age),group=if_else(is.na(nexta),paste0(age,'-'),paste0(age,'-',nexta-1)),type=type) %>% ungroup() 
    a<-   left_join(a,s,by = join_by(Var1)) %>%
      select(sp.no,param,Var,age,group,Value,type)
   if (type=='fleet')   a<-a %>% mutate(sp.no= data$keySurvey.overview[sp.no,'s'])
   a
  }
  
  k<-rbind(
    cv_par(par='logSdLogObsCatch',key=data$keyVarObsCatch,type='species'),
    cv_par(par='logSdLogFsta',key=data$keyLogFstaSd,type='species'),
    cv_par(par="logCatchability",key=data$keyCatchability,type='fleet'),
    cv_par(par="logSdLogObsSurvey",key=data$keyVarObsSurvey,type='fleet'),
    cv_par(par="logSdLogN",key=data$keyVarLogN,type='species'),
    data.frame(param="rho", sp.no=1:data$nSpecies, Var=1:data$nSpecies,age=-9,group='-9',Value=as.integer(data$useRho),type='species') %>% filter(Value>0) %>%
      mutate(Value=cumsum(Value))%>%as_tibble()
  )
  
  pf<-sdrep[['par.fixed']]
  pcv<-data.frame(param=names(pf),parValue=pf,sd=sqrt(diag(sdrep[['cov.fixed']])),gradient=sdrep$gradient.fixed) %>% 
    mutate(cv=round(abs(sd/parValue),3),mark=if_else(cv>0.9,'check','')) %>%
    group_by(param) %>% mutate(n=dplyr::row_number()) %>% ungroup() %>% data.frame()
  pcv
 
  cv<-left_join(pcv,k,by=join_by(param==param,n==Value)) %>% 
    mutate(vari=if_else(type=='species',data$spNames[Var],rownames(data$keySurvey.overview)[Var])) %>%
    select(sp.no,param,parValue,sd, vari,gradient, ages=group) %>%arrange(desc(sd)) %>% as_tibble()
  cv
}


# print(calcCV(sdrep),n=40)



readFleetCatch<-function(fleetControl,dir=data.path,of='fleet_info.dat', so='fleet_catch.in',condense=TRUE,nSpecies=1,verbose=FALSE) {
  if (missing(fleetControl)) {
    a<-readFleetInfo(dir=dir,of=of,off.age=1,off.year=1,nSpecies,verbose)
  } else a<-fleetControl
  
  k<-a[['k']]
  notUse<-k[k[,'useFleet']==0,'f']
  
  nFleets<-nrow(k)
  nSpecies<-length(unique(k[,'s']))
  fleetNames<-rownames(k)
  b<-scan(file=so,quiet=TRUE,comment.char = "#")
  i<-1
  out<-lapply(1:nFleets,function(x) {
    ny<-k[x,"maxy"]-k[x,'miny']+1L
    na<-k[x,"maxage"]-k[x,'minage']+1L
    
    if (k[x,'type'] %in% c(2L,3L)) na<-1L
    if (k[x,'type'] %in% c(4L,5L)) na<-0L
    ii<-ny*(na+1)
    #cat(x,i,ii,ny,na,k[x,'type'],'\n')
    bb<-b[i:(i+ii-1)]
    bb<-matrix(bb,nrow=ny,ncol=na+1,byrow=T)
    
    if (na>0) colnames(bb)<-c('effort',paste('age',k[x,'minage']:k[x,'maxage'],sep='_')) else colnames(bb)<-'effort'
    rownames(bb)<-paste('y',k[x,'minyear']:k[x,'maxyear'],sep='_')
    
    if (condense) {
      if (na>0) bb<- bb[,2:(na+1),drop=FALSE]  /bb[,1]  #cpue
      bb<-array2DF(bb) %>% 
        transmute(year=as.integer(substr(Var1,3,6)),q=k[x,'q'],age=substr(Var2,5,6),s=k[x,'s'],
                  f=k[x,'f'],type=k[x,'type'],obs=Value) %>% as_tibble()
    } 
    i<<-i+ii 
    return(bb)
  })
  if (length(notUse)>0) {
    out<-out[-notUse]
    for (f in (1:length(out))) {
      out[[f]]$f<-f
    }
  }
  if (condense) {
    out<-do.call(rbind,out) 
    out[out$age=='rt','age']<- "-999"
    out$age<-as.integer(out$age)
  }
  out
}

# readFleetCatch(condense=TRUE,verbose=TRUE)
#readFleetCatch(condense=FALSE)



index_old_new<-function(dir,outDir=dir,sms.dat='sms.dat',of="fleet_info.dat",verbose=FALSE) {
  if (FALSE) {
    dir=data.path; sms.dat='sms.dat'; seasFfrom<-'catch'; of="fleet_info_template.dat"
  }
 
    Init.function(dir=dir,sms.dat=sms.dat) # initialize SMS environment
  #cat(SMS.control@species.names,'\n') # just checking
  
  sms<-SMS.control  # just shorter name
  info<-sms@species.info
  multi.environment<- any(info[,'predator']==2)
  
  nSpecies<-sms@no.species-first.VPA+1L  # species with analytical assessment
  nOthSpecies<-0L
  ages<-sms@first.age:sms@max.age.all
  nAges<-length(ages)
  recSeason<-as.integer(sms@rec.season)
  nSeasons<-as.integer(sms@last.season)
  years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)
  info<-sms@species.info[first.VPA:sms@no.species,c("last-age", "first-age F>0", "last-age-selec", "last-age-likelihood", "+group", "SSB/R","RecAdd2"),drop=FALSE]
  spNames<-dimnames(info)[[1]]
  othspNames<-dimnames(info)[[1]][1:(first.VPA-1)]
  
  
  # read survey data and options into FLR objects
  indices<-SMS2FLIndices(control=sms,path=file.path(dir))
  #indices[[1]]@range.SMS
  
  nFleets<-length(indices)
  spNo<-unlist(lapply(indices,function(x) x@range.SMS["species"]))
  PowerAge<-unlist(lapply(indices,function(x) x@range.SMS["power.age"]))
  season<-unlist(lapply(indices,function(x) x@range.SMS["season"]))
  PowerAge[PowerAge<0]<- -9
  PowerAge[PowerAge>=0]<- PowerAge[PowerAge>=0] +1L
  fleetNames<-unlist(lapply(indices,function(x) x@name))
  q.age<-unlist(lapply(indices,function(x) x@range.SMS['q.age']))
  var.age<-lapply(indices,function(x) x@range.SMS['var.age.group'])
  a<-do.call(data.frame,lapply(indices,function(x) x@range))
  a<-rbind(a,s=spNo,q.age=q.age,f=1:nFleets,PowerAge=PowerAge,q=season)
  a["plusgroup",]<-0L
  ra<-rownames(a)
  ra[1]<-'minage'
  ra[2]<-'maxage'
  rownames(a)<-ra
  
  a<- a%>% mutate_if(is.numeric,as.integer)
  a<-t(a)[,c('f','s','minyear',"maxyear","minage","maxage","q.age","plusgroup","PowerAge","q","startf","endf")]
  a<-cbind(a,type=1L,techCreep=0L)
  
  if (is.vector(a)) {lena<-length(a); nama=names(a); a<-matrix(a,nrow=1,ncol=lena); colnames(a)<-nama}
  rownames(a)<-paste(1:nFleets,fleetNames)
  
  keySurvey<-a
  keySurvey.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName") %>% mutate(species=spNames[s])
  
  ## fleets per species
  a<-keySurvey.df %>% group_by(s,species) %>% summarize(n=dplyr::n()) %>% data.frame()
  
  cat("# File for configuration of survey observations\n",file=of,append=F)
  cat(0.2, " # minimum CV of CPUE observations\n",file=of,append=T)
  cat("### number of fleets by species\n",file=of,append=T) 
  for (i in (1:dim(a)[[1]])) cat(a[i,'n']," # ",a[i,'species'],'\n',file=of,append=T)
  
  cat("####### survey options ###### \n",
      "#      Fleet name\n",
      "# 1:   Use the fleet asssessment, 1=yes, 0=no",
      "# 2:   Index type and input:\n",
      "#        1 = Effort and index catch in numbers by age;\n",
      "#        2 = Effort and index for total exploitable biomass, (including all ages with F>0);\n",
      "#        3 = Effort and index for toal spawning stock biomass;\n",
      "#        4 = Commercial effort used with the assumption that Effort is a proxy for F at age;\n",
      "#        5 = Commercial effort used with the assumption that Effort is a proxy for Fbar ;\n",
      "# 3-4:  First and last year\n",
      "# 5:    Season for survey, use 1 for annual data\n",
      "# 6.7:  Alpha and beta - the start and end of the survey period for the fleet given as fractions of the season\n",
      "#                  (or of the year if annual data are used)\n",
      "# 8-9:  First and last age, use -1 if not relevant\n",
      "# 10:    Last age is a plus group? (0=no, 1=yes)\n",
      "# 11:   Last age for stock size dependent catchability (power model), -1 shows no ages uses power model\n",
      "# 12:   Technical creep in effort (only used for type 4 index)\n",
      "# 13:   Number  of cathability groups\n",
      "#         First ages in group of ages  for estimates of catcahability at age\n",
      "#       Number of variance of cathability groups\n",
      "#         First ages in group of ages  for estimates of variance of catcahability at age\n",
      "#       Check value (must be -999)\n",
      "########################",file=of,append=T)
  
  txt<-FALSE 
  ff<-0
  for (sp in (1:dim(a)[[1]])) {
    cat("#########  ",a[sp,'species'],'\n',file=of,append=T)
    fls<-filter(keySurvey.df,sp==s) 
    for (i in (1:dim(fls)[[1]])) {
      ff<-ff+1
      cat("##  ",fls[i,'species'],"# Fleet:",i,fleetNames[fls[i,'f']], "  fleet",ff,'\n',file=of,append=T)
      cat(fls[i,'useFleet'],"\t\t# Use the fleet\n",file=of,append=T)
      cat(fls[i,'type'],"\t\t# Index type\n",file=of,append=T)
      cat(fls[i,'minyear'],fls[i,'maxyear'],"\t# First and last year\n",file=of,append=T)
      cat(fls[i,'q'],"\t\t# Season for survey\n",file=of,append=T)
      cat(fls[i,'startf'],fls[i,'endf'],"\t\t# Alpha and beta\n",file=of,append=T)
      cat(fls[i,'minage'],fls[i,'maxage'],"\t\t# First and last age\n",file=of,append=T)
      cat(0,"\t\t# Is last age a plus group (0=no, 1=yes)\n",file=of,append=T)
      cat(fls[i,'PowerAge'],"\t\t# Last age power model\n",file=of,append=T)
      cat("0 \t\t# Technical creep\n",file=of,append=T)
      nag<-fls[i,'minage']:fls[i,'q.age']
      if (fls[i,'q.age']<fls[i,'maxage']) nag<-c(nag,tail(nag,1)+1)
      cat(length(nag),"\t\t# Number of cathability groups\n",file=of,append=T)
      cat("  ",nag,if_else(length(nag)<3,"\t\t","\t"),"# First ages, catchability groups\n",file=of,append=T)
      va<-var.age[[ff]][[1]]
      cat(length(va),"\t\t# Number of variance of catchability groups\n",file=of,append=T)
      cat("  ",va,if_else(length(va)<3,"\t\t","\t"),"# First ages, variance of catchahability at age\n",file=of,append=T)
      cat("-999 \t\t# check value (must be -999)\n",file=of,append=T)
    }
  }
 cat("file: ",of,"has been written\n")
}

# index_old_new(dir,outDir=dir,sms.dat='sms.dat',of=file.path(root,dir,"fleet_info_new.dat")) 



readFleetInfo<-function(dir=data.path,of='fleet_info.dat',off.age,off.year,nSpecies=1,verbose=FALSE){
  
  #   a<-readFleetInfo(dir=dir,of=of,off.age=1,off.year=1)

  # dir=data.path; of='fleet_info.dat'; off.age=1;off.year=-1973; verbose=TRUE

  a<-scan(file=file.path(dir,of),comment.char = "#",quiet=TRUE)
  cv<-a[1]
  a<-a[-1]
  nf<-a[1:nSpecies]
  a<-a[-(1:nSpecies)]
  
  
  nFleets<-sum(nf)
  # Fleetnames
  s<-readLines(file.path(dir,'fleet_names.in'), n=1000)
  s<-gsub('_',' ',s)
  s<-sub('[[:space:]]+$', '', s)
  s<-s[substr(s,1,1)!='#']
  flNames<-s[1:nFleets]
  
  
  k<-matrix(-999L,nrow=nFleets,ncol=18)
  colnames(k)<-c("f","s","minyear","maxyear","miny","maxy","minage","maxage","mina","maxa","plusgroup","PowerAge", "q","startf","endf","type","techCreep","useFleet")
  rownames(k)<-flNames
  
  k[,'f']<-1:nFleets
  fl<-0
  catchability<-list()
  catchability.var<-list()
  
  for (s in (1:nSpecies)) {
    if (verbose) cat('reading species:',s,'\n')
    for (nfl in(1:nf[s]) ) {
      if (verbose) cat('\treadingfleet:',nfl)
      i<-1
      fl<-fl+1
      k[fl,'s']<-s
      k[fl,'useFleet']<-a[i]; i<-i+1
      k[fl,'type']<-a[i]; i<-i+1
      k[fl,'minyear']<-a[i]; i<-i+1
      k[fl,'maxyear']<-a[i]; i<-i+1
      k[fl,'q']<-a[i]; i<-i+1
      k[fl,'startf']<-a[i]; i<-i+1
      k[fl,'endf']<-a[i]; i<-i+1
      k[fl,'minage']<-a[i]; i<-i+1
      k[fl,'maxage']<-a[i]; i<-i+1
      k[fl,'plusgroup']<-a[i]; i<-i+1
      k[fl,'PowerAge']<-a[i]; i<-i+1
      k[fl,'techCreep']<-a[i]; i<-i+1
      n.catcha<-a[i];i<-i+1
      catchability<-c(catchability,list(a[i:(i+n.catcha-1)])); i<-i+n.catcha
      n.catcha<-a[i];i<-i+1
      catchability.var<-c(catchability.var,list(a[i:(i+n.catcha-1)])); i<-i+n.catcha
      chk<-a[i]; 
      if (verbose) cat(' done, check=:',chk,'\n')
      stopifnot(chk== -999)
      a<-a[-(1:i)]
    }
  } 
  
  k[,'mina']<-k[,'minage']+off.age
  k[,'maxa']<-k[,'maxage']+off.age
  k[,'miny']<-k[,'minyear']+off.year
  k[,'maxy']<-k[,'maxyear']+off.year
  names(catchability)<-rownames(k)
  names(catchability.var)<-rownames(k)
  #mode(k)<-'integer'
  list(k=k,cv=cv,catchability=catchability,catchability.var=catchability.var)
}



writeFleetInfo<-function(dir=dir,of='fleet_info.dat',key,speciesNames){
  
  a<-key[['k']]
  catcb<-key[['catchability']]
  catcbVar<-key[['catchability.var']]
  key.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName") %>% mutate(species=speciesNames[s])
  
  ## fleets per species
  b<-key.df %>% group_by(s,species) %>% summarize(n=dplyr::n(),.groups='drop') %>% ungroup() %>% data.frame()
  
  cat("# File for configuration of survey observations\n",file=of,append=F)
  cat(key$cv, " # minimum CV of CPUE observations\n",file=of,append=T)
  cat("### number of fleets by species\n",file=of,append=T) 
  for (i in (1:dim(b)[[1]])) cat(b[i,'n']," # ",b[i,'species'],'\n',file=of,append=T)
  
  cat("####### survey options ###### \n",
      "#      Fleet name\n",
      "# 1:   Use the fleet asssessment, 1=yes, 0=no\n",
      "# 2:   Index type and input:\n",
      "#        1 = Effort and index catch in numbers by age;\n",
      "#        2 = Effort and index for total exploitable biomass, (including all ages with F>0);\n",
      "#        3 = Effort and index for toal spawning stock biomass;\n",
      "#        4 = Commercial effort used with the assumption that Effort is a proxy for F at age;\n",
      "#        5 = Commercial effort used with the assumption that Effort is a proxy for Fbar ;\n",
      "# 3-4:  First and last year\n",
      "# 5:    Season for survey, use 1 for annual data\n",
      "# 6.7:  Alpha and beta - the start and end of the survey period for the fleet given as fractions of the season\n",
      "#                  (or of the year if annual data are used)\n",
      "# 8-9:  First and last age, use -1 if not relevant\n",
      "# 10:    Last age is a plus group? (0=no, 1=yes)\n",
      "# 11:   Last age for stock size dependent catchability (power model), -1 shows no ages uses power model\n",
      "# 12:   Technical creep in effort (only used for type 4 index)\n",
      "# 13:   Number  of cathability groups\n",
      "#         First ages in group of ages  for estimates of catcahability at age\n",
      "#       Number of variance of cathability groups\n",
      "#         First ages in group of ages  for estimates of variance of catcahability at age\n",
      "#       Check value (must be -999)\n",
      "########################",file=of,append=T)
  
  txt<-FALSE 
  ff<-0
  for (sp in (1:dim(b)[[1]])) {
    cat("#########  ",b[sp,'species'],'\n',file=of,append=T)
    fls<-filter(key.df,s==sp) 
    for (i in (1:dim(fls)[[1]])) {
      ff<-ff+1
      cat("##  ",fls[i,'species'],"# Fleet:",i,fls[i,'fName'], "  fleet",ff,'\n',file=of,append=T)
      cat(fls[i,'useFleet'],"\t\t# Use the fleet\n",file=of,append=T)
      cat(fls[i,'type'],"\t\t# Index type\n",file=of,append=T)
      cat(fls[i,'minyear'],fls[i,'maxyear'],"\t# First and last year\n",file=of,append=T)
      cat(fls[i,'q'],"\t\t# Season for survey\n",file=of,append=T)
      cat(fls[i,'startf'],fls[i,'endf'],"\t\t# Alpha and beta\n",file=of,append=T)
      cat(fls[i,'minage'],fls[i,'maxage'],"\t\t# First and last age\n",file=of,append=T)
      cat(0,"\t\t# Is last age a plus group (0=no, 1=yes)\n",file=of,append=T)
      cat(fls[i,'PowerAge'],"\t\t# Last age power model\n",file=of,append=T)
      cat("0 \t\t# Technical creep\n",file=of,append=T)
      nag<-length(catcb[[ff]])
      cat(nag,"\t\t# Number of cathability groups\n",file=of,append=T)
      cat("  ",catcb[[ff]],if_else(length(nag)<3,"\t\t","\t"),"# First ages, catchability groups\n",file=of,append=T)
      va<-catcbVar[[ff]]
      cat(length(va),"\t\t# Number of variance of catchability groups\n",file=of,append=T)
      cat("  ",va,if_else(length(va)<3,"\t\t","\t"),"# First ages, variance of catchahability at age\n",file=of,append=T)
      cat("-999 \t\t# check value (must be -999)\n",file=of,append=T)
    }
  }
  cat("file: ",of,"has been written\n")
}


#a<-readFleetInfo(dir=data.path,of='fleet_info.dat',off.age=data$off.age,off.year=data$off.year,nSpecies=data$nSpecies)
# writeFleetInfo(dir=data.path,of='fleet_info3.dat',key=a,speciesNames=data$spNames)



announce<-function(x) {
  cat("\nobjective:",x$objective,"  convergence:",x$convergence, "  ", x$message, "  iterations:",x$iterations, "  evaluations:",x$evaluations,'\n')
}

inputToDF<-function(data) {
  a<-rbind(
  #  list_rbind(lapply(data$catchMeanWeight,array2DF),names_to='species') %>% mutate(var="weca"),
    list_rbind(lapply(data$catchNumber,array2DF),names_to='species') %>% mutate(var="canum"),
    list_rbind(lapply(data$stockMeanWeight,array2DF),names_to='species') %>% mutate(var="west"),
    list_rbind(lapply(data$propMat,array2DF),names_to='species') %>% mutate(var="propMat"),
    list_rbind(lapply(data$propM,array2DF),names_to='species') %>% mutate(var="propM"),
    list_rbind(lapply(data$propF,array2DF),names_to='species') %>% mutate(var="propF"),
    list_rbind(lapply(data$natMor,array2DF),names_to='species') %>% mutate(var="M"),
    list_rbind(lapply(data$seasFprop,array2DF),names_to='species') %>% mutate(var="seasFprop"),
    list_rbind(lapply(data$natMor1,array2DF),names_to='species') %>% mutate(var="M1")


  )
  pivot_wider(a,names_from=var,values_from=Value) %>%
    mutate(year=as.integer(Var1)-data$off.year, quarter=as.integer(Var2),age=as.integer(Var3)-data$off.age, 
           Var1=NULL,Var2=NULL,Var3=NULL,s=match( species,data$spNames))
}


#a<-inputToDF(data)

saveResults<-function(runName='nr1',data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep) {
  a<-inputToDF(data)
  x<- rep$res %>% select(s,year,quarter,age,FisQ) %>% group_by(s,year,age)%>% mutate(FiProp=FisQ/sum(FisQ)) %>% mutate(FiProp=if_else(is.na(FiProp),0,FiProp),FisQ=NULL)
  a<-left_join(a,x,by = join_by(s, year, quarter, age))
  rep$res<-left_join(rep$res,a,by = join_by(s,year, quarter, age, species))
  
  rep$residSurv  <-rep$residSurv %>% mutate(fleetName=data$fleetNames[f],year=y-data$off.year,age=a-data$off.age,species=data$spNames[s],resid=logObs-logPred)
  rep$residCatch <-rep$residCatch %>% mutate(year=y-data$off.year,age=a-data$off.age,species=data$spNames[s],resid=logObs-logPred)
  
  sms<-list(data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=map,random=random,rep=rep,sdrep=sdrep)
  save(sms,file=file.path(data.path,paste0(runName,".Rdata")))
  print(object.size(sms),units="Mb"); #lapply(sms,function(x) print(object.size(x),units="Mb"))
  invisible(sms)
}


outputAnnuToDF<-function(obj) {
  rep<-obj$report()
  a<-rbind(
    list_rbind(lapply(rep$Chat,array2DF),names_to='s') %>% mutate(var="Chat")
  )  
  pivot_wider(a,names_from=var,values_from=Value)  %>%
    mutate(species=s, year=as.integer(y),age=parse_number(a)-data$off.age) %>% 
    mutate(y=NULL,a=NULL,Var3=NULL,s=NULL,Chat=exp(Chat))
} 
#  b<-outputAnnuToDF(obj)

outputToDF<-function(obj) {
  rep<-obj$report()
  a<-rbind(
    list_rbind(lapply(rep$logNq,array2DF),names_to='s') %>% mutate(var="N",Value=exp(Value)),
    list_rbind(lapply(rep$Zq,array2DF),names_to='s') %>% mutate(var="Z")
   )  
     pivot_wider(a,names_from=var,values_from=Value)  %>%
    mutate(species=s, year=as.integer(y),q=as.integer(q),age=parse_number(a)-data$off.age) %>% 
    mutate(y=NULL,a=NULL,Var3=NULL,s=NULL,deadZ=N*(1-exp(-Z)))
} 
#  b<-outputToDF(obj)


FToDF<-function(obj,sdrep,data) {
  if (missing(sdrep)) sdrep <- sdreport(obj)
  x<-as.list(sdrep, what="Est")
  ff<-exp(x$Uf)
  colnames(ff)<-data$years
  FF<-NULL
  for (s in (1:data$nSpecies)) {
    i<-data$nlogFfromTo[s,]
    key<-data$keyLogFsta[s,][data$keyLogFsta[s,]>0]
    faf<-data$info[s,'faf']; la<-data$info[s,'la']
    fff<-ff[i[1]:i[2],,drop=FALSE][key,]
    if (data$info[s,'fModel']==2) {
      fSepar<-data$info[s,'fSepar']; la<-data$info[s,'la']
      key2<-t(data$keyLogSeparF[[s]][,fSepar:la] )
      dims<-dim(key2)
      seaF<-exp(x$logSeparF[key2])
      seaF<-matrix(seaF,ncol=dims[[2]])
      fff[fSepar:la,]<-fff[fSepar:la,]*seaF
    }
    rownames(fff)<-(faf:la)-data$off.age
    fff<-array2DF(fff) %>% mutate(s=s,year=as.integer(Var2),age=as.integer(Var1),Var1=NULL,Var2=NULL) %>%rename(FF=Value)
    
    FF<-rbind(FF,fff)
  }
  FF
}

seasonalF<-function(obj,sdrep,data) {
  FF<-FToDF(obj,sdrep,data) 
  a<-inputToDF(data) %>% select(s,seasFprop, species, year, q, age)
  left_join(FF,a,by = join_by(s, year, age)) %>% mutate(SF=FF*seasFprop)
}
#  FF<-FToDF(obj,sdrep=sdrep)

#b<-left_join(inputToDF(data),outputToDF(obj),by = join_by(species, year, q, age))

plotF_age<-function(obj,sdrep,data,combineAges=FALSE,pLabel='') {
  ff<- FToDF(obj,sdrep,data) %>% mutate(species=data$spNames[s])
  
  a<-t(data$keyLogFsta)
  a<-cbind(a,Age=data$minAge-1 +(1:data$nAges))
  a<-data.frame(a) %>% pivot_longer(cols=1:data$nSpecies) %>%rename(species=name,aGroup=value,age=Age) %>%   mutate(species=gsub('\\.',' ',species))
  
  ff<-left_join(ff,a,by = join_by(age, species)) %>% as_tibble()
  
  if (data$zeroCatchYearExists)  {
    for (s in (1:data$nSpecies)) if (data$zeroCatchYearExistsSp[s]==1) {
      zcy<-data$zeroCatchYear[[s]]-data$off.year
      crit<- ff$year %in% zcy & ff$s==s
      ff[crit,'FF']<-0.0
    }
  }
  
  
  ag<-a%>%  group_by(species,aGroup) %>% summarise(mina=min(age),maxa=max(age)) %>% ungroup() %>%
    mutate(ages=paste(mina,maxa,sep='-'),mina=NULL,maxa=NULL)
  ag
  
  ff<-left_join(ff,ag,by = join_by(species, aGroup))%>% mutate(ages=factor(ages),age=factor(age))
  
  fplt<-by(ff,ff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=age,col=age))+
             geom_point(size=2)  + geom_line()+labs(title=paste(x[1,'species'],pLabel),ylab='F')+
             facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
  )

  if (combineAges) {
    fff<-ff %>%mutate(age=NULL) %>% unique()
    
    fplt<-by(fff,fff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=ages,col=ages))+
         geom_point(size=2)  + geom_line()+labs(title=paste(x[1,'species'],pLabel),ylab='F')+
         facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
    )
  }
  return(fplt)
}


tableNoObs<-function(runName,doPrint=TRUE) {
  load(file=file.path(data.path,paste0(runName,".Rdata")))
  noObs<-matrix(0,ncol=4,nrow=sms$data$nSpeciesAll ,dimnames=list(species=sms$data$allSpNames,type=c("Catch","Survey",'SSB/R','Stomachs')))
  
  noObs[,"Catch"]<-xtabs(~s,data=as.data.frame(data$keyCatch))
  noObs[,"Survey"]<-xtabs(~s,data=as.data.frame(data$keySurvey))
  
  if ("logSsbRsd" %in%   names(as.list(sms$sdrep, what="Estimate"))) {
    noObs[,"SSB/R"]<-apply(sms$data$recruitYears,1,sum)
  }
  if (doPrint) print(noObs) 
  invisible(noObs)
}
