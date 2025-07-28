tabsOut<-c('logYearEffectF')



inp=runName
verbose=TRUE
load(file=file.path(data.path,paste0(inp,".Rdata")),verbose=verbose)

str(sms,1)
str(sms$data,1)
sms$random

str(sms$sdrep,1)

a<-extractParameters(sms$sdrep,sms$map,sms$data)[[2]]
sort(unique(a$name))

b<-filter(a,name=="logYearEffectF") %>% transmute(s, Species=sms$data$spNames[s],Year=min_age,expEst=exp(estimate)) 



print(a,n=600)

incl_ICES_plot<-FALSE
t(t(opt$par));
exp(t(t(opt$par)))
rep<-obj$report()

overw<-as.data.frame(data$keySurvey.overview) %>%  rownames_to_column(var = "fleet") %>% select(f,fleet,s,type)

if (FALSE) {
  source(file.path(rsms.root.prog,"survey_residuals.R"))
  source(file.path(rsms.root.prog,"catch_residuals.R"))
  source(file.path(rsms.root.prog,"process_noise.R"))
}

#round(ftable(data$catchNumber[[1]][,,]/1000))

#lapply(rep$Zq,function(x) round(exp(x[,1,]),2))
#lapply(rep$Chat,function(x) round(exp(x)))

cbind(rep$nlls,all=rowSums(rep$nlls))

Est<-as.list(rep, "Est", report=TRUE)

if (!'sdrep' %in% ls()) sdrep <- sdreport(obj)
sdrep$pdHess
sdrep

# obj$fn()  #  value er den samme som sum af min nnls

x<-as.list(sdrep, what="Est")
x$Uf[,1:4]


#summary(sdrep, "random")                      ## Only random effects
#summary(sdrep, "fixed", p.value = TRUE)       ## Only non-random effects
#summary(sdrep, "fixed", p.value = FALSE)    ## Only non-random effects
#summary(sdrep, "report")                      ## Only report
# 
# ssb<-summary(sdrep, "report") 
# plot(ssb[,"Estimate"],type='l')
# lines(ssb[,"Estimate"]+2*ssb[,"Std. Error"],type='l',col='red')
# lines(ssb[,"Estimate"]-2*ssb[,"Std. Error"],type='l',col='red')

aa<-summary(sdrep, "fixed", p.value = FALSE) 
aa<-as.data.frame(aa) %>% mutate(n=1:dplyr::n(),param=rownames(aa))

parVal<-as.list(sdrep, what="Est")
parSd<-as.list(sdrep, what="Std")

make_tab<-function(d,key,printIt=TRUE,roundIt=2,addSd=TRUE,expIt=TRUE) {
  cat(d,'\n')
  dd<-parVal[[d]]
  if (expIt) dd<-exp(dd)
  tab<-key
  for (i in 1:length(dd)) tab[tab[,]==i]<- dd[i]
  if (addSd) {
    dd<-parSd[[d]]
    if (expIt) dd<-exp(dd)
    tabSd<-key
    for (i in 1:length(dd)) tabSd[tabSd[,]==i]<- dd[i]
    tabSd
  }
  
  if (printIt) {
    ptab<-tab
    ptab[ptab== -1] <-NA
    print(round(ptab,roundIt),na.print='.')
    if (addSd) {
      ptab<-tabSd
      ptab[ptab== -1] <-NA
      print(round(ptab,roundIt),na.print='.')
    }
    
  } 
  invisible(tab)
}

addSd<-TRUE;expIt<-T
make_tab(d="logSdLogN",key=data$keyVarLogN,roundIt=2,addSd=addSd,expIt=expIt)


#round(exp(x$Uf),3)
#round(exp(x$Un),0)

x$rec_loga; exp(x$rec_loga)
x$rec_logb; exp(x$rec_logb)

make_tab<-function(d,key,printIt=TRUE,roundIt=2) {
  d<-exp(d)
  tab<-key
  for (i in 1:length(d)) tab[tab[,]==i]<- d[i]
  tab
  if (printIt) {
    ptab<-tab
    ptab[ptab== -1] <-NA
    print(round(ptab,roundIt),na.print='.')
  } 
  invisible(tab)
}


make_tab3<-function(d,key,printIt=TRUE,roundIt=2) {
  for (s in (1:data$nSpecies)) {
    if (data$info[s,'fModel']==2) {
      key2<-key[[s]] 
      dims<-dim(key2)
      seaF<-exp(x$logSeparF[key2])
      seaF<-matrix(seaF,ncol=dims[[2]])
      rownames(seaF)<-data$years
      seaF<-seaF[which(!duplicated(seaF)),,drop=FALSE]
      faf<-data$info[s,'faf']; la<-data$info[s,'la']
      colnames(seaF)<-paste('age',(faf:la)-data$off.age)
      cat(data$spNames[s],'\n')
      print(round(seaF,5))
    }
  }
}
  


make_tab(d=x$logSdLogN,key=data$keyVarLogN,roundIt=2)
make_tab(d=x$logSdLogFsta,key=data$keyLogFstaSd,roundIt=2)
make_tab(d=x$logSdLogObsCatch,key=data$keyVarObsCatch)  
make_tab(d=x$logSdLogObsSurvey,key=data$keyVarObsSurvey) 
make_tab(d=x$logCatchability,key=data$keyCatchability) 

make_tab3(d=x$logSeparF,key=data$keyLogSeparF,roundIt=2)




summary_plot(obj,out=as.list(sdrep, "Est", report=TRUE),data,sdrep,incl_ICES_plot=F) 
plotF(obj,sdrep,data,combineAges=FALSE)  
  
ssb<-tail(summary(sdrep, "report"),1); cat('CV of SSB last year:',round(ssb[,2]/ssb[,1],3),'\n')


if (FALSE) for (s in (rsp)) {
  bb=filter(b,Species==s & source=='rsms')
  
  p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
    geom_line(aes(linetype=source,col=source))+
    geom_point(aes(shape=source,col=source))+
    facet_grid(variable ~ ., scales="free_y")+
    ggtitle(unlist(bb[1,'Species']))
  print(p)
  cat('press return to see the next plot:')
  readLines(n=1)
  cat('\n')
}



#############################
#Once obj has been constructed successfully, you should evaluate it

if (FALSE) obj$fn()

#and verify that it gives the same result as func(parameters).

if (FALSE) {  # give errors (as also seen in the SAM script from Anders)
  set.seed(1)
  chk <- checkConsistency(obj)
  chk
}

if (FALSE) {
  debug(func)
  obj <- MakeADFun(func, parameters, random=c("Un","Uf"),silent=FALSE)
}
