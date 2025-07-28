if (FALSE){ # run if you change SMS 
  fromDir<-"rsms_input_sms_keyRun_new_log2PI"
  file.copy(from=file.path(root,fromDir,'canum.in'),to=file.path(data.path,'canum.in'),overwrite=TRUE)
  file.copy(from=file.path(root,fromDir,'weca.in'),to=file.path(data.path,'weca.in'),overwrite=TRUE)
}

if (TRUE) {
  #compareSMS<-"rsms_input_sms_keyRun_new_log2PI"
  rsmsFile<-'rsms.dat'
  doMultiExtract<-FALSE
  runName<-'RSMS_SMS'

  rsms<-batch_SMS_like_configuration(outfile=rsmsFile,writeConrol=F)
  
  #                 COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  Nrr<-as.integer(c( 2,    2,    2,    2,    2,    2,    3,    2,    2,    3,    2,    2) )
  Nff<-as.integer(c( 2,    2,    2,    2,    2,    3,    3,    2,    3,    3,    2,    2) ) 
  
  rel.tol.mul<-1e-10  # Relative tolerance. Defaults to 1e-10. JEG KAN IKKE SE NOGEN EFFECT AF HØJERE VÆRDI
  xf.tol.mul<-2.2E-14   # false convergence tolerance. Defaults to 2.2e-14
  doMultiExtract<-FALSE
  
  fleetInfoFile<-"fleet_info.dat"
  fleetTemplate<-"fleet_info_template.dat"
  
  rsms<-batch_SMS_like_configuration(outfile=rsmsFile,writeConrol=F)
  rsms<-fineTune(control=rsms, writeControl=FALSE,outfile=rsmsFile,writeSurvey=TRUE, surveyIn=fleetTemplate, surveyOut=fleetInfoFile) 
  rsms<-batch_SMS_like_recState_Frandom_configuration(control=rsms,outfile=rsmsFile,writeConrol=T,Nr=Nrr,Nf=Nff) 
  
  
  inspect<-TRUE
}


### Extract data from SMS

if (TRUE) {  # transform  SMS data into RSMS format 
  switch(my.stock.dir,
   "rsms_input"       = inp_all<-make_rsms_data(dir=data.path,rsmsFile='rsms.dat',multi=doMultiExtract),
   "rsms_SAN-area-1r"       = inp_all<-make_rsms_data(dir=data.path,rsmsFile='rsms.dat',multi=FALSE),
    stop(paste0("Not valid stock dir: ",my.stock.dir))
  ) #end switch
  save(inp_all,file=file.path(data.path,"rsms_input_all.Rdata")) 
} else load(file=file.path(data.path,"rsms_input_all.Rdata"),verbose=TRUE)


smsConf<-0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated

# select a combination of species from the (full) data set, also including multi species information


#my.ps<-c(1,2,3,4,5,6,7,8,9,10)
my.ps<-c(1,2,3,4,5,6,7,8,  11, 12)  #virker
#my.ps<-c(11,12)
my.pso<-c(0L)
my.pso<-13L:27L


inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
#inp=inp_all

inp$data$sms.mode<-smsConf
inp$data$spNames

#### prepare to run
data<-inp[['data']]
data$silent<-FALSE
parameters<-inp[['parameters']]

if (FALSE) {
 load(paste0(runName,'.Rdata'),verbose=TRUE)
  parameters<-as.list(sms$sdrep, what="Est")  #parameters and random effects
  func(parameters)
 # osa <- oneStepPredict(obj, method="fullGaussian", discrete=FALSE)
}

#### Run rsms
cleanrun(silent=TRUE)

myMap<-map_param(data,parameters)

if (length(myMap$Un)>0) random<-c("Un") else random<-NULL
if (length(myMap$Uf)>0) random<-c(random,"Uf")

system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))
#checkConsistency(obj);

lu<-lowerUpper(obj,data,parameters )

system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=6000,eval.max=6000)))
announce(opt)

system.time(sdrep <- sdreport(obj, ignore.parm.uncertainty = F, getReportCovariance = FALSE))
cat('Hesssian:',sdrep$pdHess,'\n')

a<-extractParameters(sdrep,myMap,data,parameters,lu=lu,showNotUsed=FALSE)[[2]]
print(a,n=300)
print(filter(a,Var1=='NOP'),n=300)
print(filter(a,is.na(estimate.sd)),n=20)
print(filter(a,gradient==0),n=20)

print(filter(a,bounds),n=20)
print(a %>% arrange(desc(abs(gradient))),n=10)
# print(filter(a, s==1) ,n=1009)

myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a,1)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
#tableNoObs(runName)



transExternalSummary(inp='summary_table_raw_SMS_KeyRun_new_log2PI0_SS.out',outSet='SMS_SS',spNames=data$spNames) #, exSpeciesNames=data$spNames)
transExternalData(inp='summary_SMS_KeyRun_new_log2PI0_SS.out',outSet='SMS_details_SS',spNames=data$spNames)
transExternalSummary(inp='summary_table_raw_SMS_KeyRun_new_log2PI2_MS.out',outSet='SMS_MS',spNames=data$spNames) #, exSpeciesNames=data$spNames)
transExternalData(inp='summary_SMS_KeyRun_new_log2PI2_MS.out',outSet='SMS_details_MS',spNames=data$spNames)


  plotSurveyResiduals(runName,inclSp=1:100, outFormat=c('screen','pdf','png')[1],standardized=F) 
  plotCatchResiduals (runName,inclSp=1:100,fileLabel='catch_residuals',mult=3,outFormat=c('screen','pdf','png')[1],longSpNames=TRUE,standardized=F,scaleY=c('fixed',"free_y")[1]) 
 
  inpRdata=list(runName,'SMS_SS')
  labels=c(runName,'SMS_SS')
  
   plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[2],showSpecies=1:12,
                        inpRdata,
                        labels,
                        outFormat=c('screen','pdf','png')[1],
                        showAges=0:10, ncol=3,allInOne=T,
                        useF=c("simpleSumqF","fromCandZ")[2],
                        longSpNames=TRUE, fileLabel='single')




###################  multi species mode 1
   
runName<-'RSMS_SMS' 
load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

runName<-'RSMS_SMS_MS'
smsConf<-2 # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
doMultiExtract<-TRUE
cleanrun(silent=FALSE)

sdrep<-sms$sdrep # from single species run
cat('Hesssian,single:',sdrep$pdHess,'\n')

newPar<-as.list(sms$sdrep, what="Est")  #parameters and random effects

#my.ps=c(1:12)
my.pso<-c(13L:27L)
#my.pso<-c(0L)
my.ps; my.pso
inp_all<-make_rsms_data(dir=data.path,rsmsFile='rsms.dat',multi=doMultiExtract)
inp<-pick_species(ps=my.ps,pso=my.pso,inp=inp_all,smsConf) # to extract multi species data and parameters

inp$data$sms.mode<-smsConf
data<-inp[['data']]

# add MS parameters
newPar$vulnera<-inp$parameters$vulnera
newPar$overlapP<-inp$parameters$overlapP
newPar$logStomObsVar<-inp$parameters$logStomObsVar
parameters<-newPar

myMap<-map_param(data,parameters)

if (length(myMap$Un)>0) random<-c("Un") else random<-NULL
if (length(myMap$Uf)>0) random<-c(random,"Uf")

#lock single species parameters,if necessary
doLock<-FALSE
if (doLock) {
  multiPar<-c("vulnera","overlapP","logStomObsVar")
  lockP<-setdiff( names(myMap),multiPar)
  myMap<-lock_param(data,parameters,myMap,lockP) 
}

#func(parameters)

system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))

lu<-lowerUpper(obj,data,parameters)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=8000,eval.max=8000)))
announce(opt)

system.time(sdrep <- sdreport(obj,ignore.parm.uncertainty = F, getReportCovariance = FALSE)); 
cat('Hesssian:',sdrep$pdHess,'\n')
#sdrep

a<-extractParameters(sdrep,myMap,data,parameters,lu=lu,showNotUsed=FALSE)[[2]]
print(filter(a,is.na(estimate.sd)),n=20)
arrange(a,desc(abs(gradient)))
filter(a,name=="logStomObsVar")
print(filter(a,bounds),n=20)
# print(a,n=200)

myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
#load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[5],showSpecies=c(1:3,6:10),
                      inpRdata=runName,
                      labels=runName,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4, ncol=3,allInOne=F,
                      useF=c("simpleSumqF","fromCandZ")[2],
                      longSpNames=TRUE, fileLabel='single')

inpRdata=list('SMS_like_new','Multi1')
labels=c('Single','Multi1')

likelihoodCompare(inpRdata,labels,doPrint=TRUE)

plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[2],showSpecies=c(1:3,6:10),
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:7, ncol=3,allInOne=F,
                      useF=c("simpleSumqF","fromCandZ")[2],
                      longSpNames=TRUE, fileLabel='single')



###################  multi species mode 2
runName<-'Multi1'; load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

runName<-'Multi2'
smsConf<-2L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated

cleanrun(silent=TRUE)

sdrep<-sms$sdrep # from previous run
cat('Hesssian:',sdrep$pdHess,'\n')
newPar<-as.list(sms$sdrep, what="Est")  #parameters and random effects

inp$data$sms.mode<-smsConf
data<-inp[['data']]
parameters<-newPar
myMap<-map_param(data,parameters)

# LOCK !!!
#myMap<-lock_param(data,parameters,myMap,lockP="logFSeasonal") 

system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))

lu<-lowerUpper(obj,data,parameters )
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=7000,eval.max=7000)))
announce(opt)

system.time(sdrep <- sdreport(obj)) 
cat('Hesssian:',sdrep$pdHess,'\n')
#sdrep


a<-extractParameters(sdrep,myMap,data,parameters,lu=lu,showNotUsed=FALSE)[[2]]
print(filter(a,is.na(estimate.sd)),n=20)
arrange(a,desc(abs(gradient)))
print(filter(a,bounds),n=20)
print(a,n=200)

a<-extractParameters(sdrep,myMap,data,parameters,lu=lu,showNotUsed=TRUE)[[2]]
print(filter(a,name=='vulnera'),n=300)

myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
# load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

a<-sms$rep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)


inpRdata=list(runName,'SMS_new_MS')
labels=c(runName,'SMS_new_MS')

inpRdata=list('SMS_like_new','Multi1','Multi2')
labels=c('Single','Multi1','Multi2')

#AICCompare(inpRdata,labels, pval=TRUE)

plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[2],showSpecies=1:10,
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:10, ncol=3,allInOne=T,
                      useF=c("simpleSumqF","fromCandZ")[2],
                      longSpNames=TRUE, fileLabel='single')

plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[5],showSpecies=c(1:3,6:10),
                      inpRdata=runName,
                      labels=runName,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:7, ncol=3,allInOne=F,
                      useF=c("simpleSumqF","fromCandZ")[2],
                      longSpNames=TRUE, fileLabel='single')


inpRdata=list(runName,'SMS_new_details_MS')
labels=c(runName,'SMS_new_details')


plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[5],showSpecies=c(1:3,6:10),
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4, ncol=2,allInOne=F,
                      useF=c("simpleSumqF","fromCandZ")[2],
                      longSpNames=TRUE, fileLabel='single')

plotSeasonalData(inp=runName,Type=c("N","F","C","M","M1",'"M2','Z',"WEST","WECA","propMat","seasFprop","FiProp")[6],
                 CombineSeason=T,
                 showSpecies=c(1:3,6:10),
                 outFormat=c('screen','pdf','png')[1],
                 showAges=0:4,
                 multN=0.001,
                 ncols=2,
                 fileLabel='pl',
                 cummulate=FALSE,  
                 longSpNames=TRUE)


hertil




inp$data$stom[1,'data'][[1]][[1]]
rep1$stom[1,'data'][[1]][[1]]

inp$data$suitIdx
inp$data$suitIdx[1,'data'][[1]][[1]]
inp$data$suitIdx[1,'data'][[1]][[1]][['data']][[2]]


inp$data$stom[1,'data'][[1]][[1]][['data']][[2]]
rep1$stom[1,'data'][[1]][[1]][['data']][[2]]

inp$data$partM2Idx
inp$data$partM2Idx[1,'data'][[1]][[1]]
inp$data$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]

inp$data$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]
rep1$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]

stom2<-unnest(rep1$stom,cols = c(data))
stom2<-unnest(stom2,cols = c(data))

stom2 %>% select(pred,prey,vulneraIdx) %>% unique() %>% arrange(pred,prey)

rep1$res
M2 for 0-gruppe i Q1 ????
  
  
xx<-rep1$res %>% group_by(species,age) %>% mutate(sumM=sum(M2,na.rm=TRUE)) %>% filter(sumM>0) %>% 
  group_by(s,species,year,age) %>% summarize(M2=sum(M2,na.rm=TRUE))  %>% ungroup() %>% mutate(age=factor(age))
xx

  
  
ggplot(xx,aes(x=year,y=M2,shape=age,col=age))+
  geom_point(size=2)  + geom_line()+labs(ylab='M2',xlab='Year')+
  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")

rep1$nlls
rep1$availFood
# p<-extractParameters(sdrep1)
#view(p[[2]])

save(parameters,data,obj1,my.map1,newPar1,lu,random,file=file.path(data.path,"multi_sp1.Rdata"))
load(file=file.path(data.path,"multi_sp1.Rdata"),verbose=T)

##########################
smsConf<-2L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
load(file=file.path(data.path,"multi_sp1.Rdata"),verbose=T)

inp$data$sms.mode<-smsConf
data<-inp[['data']]

parameters<-newPar1

my.map<-map_param(data,parameters)
my.map2<-my.map
lockParameters<-FALSE
if (lockParameters) {
  #lockP<-c('logCatchability','logSdLogN','logSdLogFsta','logSdLogObsSurvey','logSdLogObsCatch','rho','rec_loga','rec_logb')
  lockP<-c('Uf','Un')
  my.map2<-lock_param(data,parameters,my.map,lockP)
}

system.time(obj2 <- MakeADFun(func, parameters, random,silent=F,map=my.map2))

lu<-lowerUpper(obj2,data,parameters)

if (exists('opt2')) rm(opt2);if (exists('sdrep2')) rm(sdrep2);
system.time(opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1300,eval.max=1300)))
announce(opt2)

sdrep2 <- sdreport(obj2,ignore.parm.uncertainty=FALSE); # ignore.parm.uncertainty=TRUE to speed up
cat('Hesssian:',sdrep2$pdHess,'\n')
sdrep2

a<-extractParameters(sdrep2)[[2]]
arrange(a,desc(abs(gradient)))


vulnerabilityTab<-function(sdrep,expIt=FALSE,asDF=FALSE) {
  vul<-sdrep[['par.fixed']]
  vul<-vul[names(vul)=='vulnera']
  if (expIt) vul<-exp(vul)
  vulnerability<-data$vulneraIdx
  vulnerability[match(1:length(vul),vulnerability)]<-vul
  vulnerability<-vulnerability[,data$preds]

  if (asDF) {
    vul.se<-sqrt(diag(sdrep$cov.fixed))
    vul.se<-vul.se[names(vul.se)=='vulnera']
    #if (expIt) vul.se<-exp(vul.se)
    vulnerability.se<-data$vulneraIdx
    vulnerability.se[match(1:length(vul.se),vulnerability.se)]<-vul.se
    vulnerability.se<-vulnerability.se[,data$preds]
    vulnerability.se<- array2DF(vulnerability.se)
    colnames(vulnerability.se)<-c('prey','pred','vulnerability.se')
    
    
    vulnerability<- array2DF(vulnerability)
    colnames(vulnerability)<-c('prey','pred','vulnerability')
    
    vulnerability<-left_join(vulnerability,vulnerability.se,by = join_by(prey, pred))
    vulnerability<-filter( vulnerability,!(vulnerability==0))
  }
  return(vulnerability)
}

vulnerabilityTab(sdrep1,expIt=TRUE,asDF=TRUE)
round(vulnerabilityTab(sdrep1,expIt=TRUE,asDF=FALSE),3) 

round(vulnerabilityTab(sdrep2,expIt=TRUE,asDF=FALSE),3) 

overlapTab<-function(sdrep,expIt=FALSE) {
  ov<-sdrep[['par.fixed']]
  ov<-ov[names(ov)=='overlapP']
  if (expIt) ov<-exp(ov)

  overlap<-data$overlapIdx
  overlap$overlap<-ov[overlap$overlap]

  return(overlap)
}

overlapTab(sdrep=sdrep1,expIt=TRUE) 
overlapTab(sdrep=sdrep2,expIt=TRUE) 


if (FALSE) {
   res<-obj2$report();  
   resul<-res$res
   resulAn<-res$resAnnual
   filter(resul,species=='WHG' & quarter==3)
   res$nlls; 
 
   M2Graph<-function(preys='COD',inp=resul,annual=TRUE,maxAgePlot=4) {
     res<-filter(inp,species %in% preys)
     res<-res %>% filter(age<=maxAgePlot) %>% group_by(species,year,age) %>% summarize(M2=sum(M2))  %>% mutate(age=factor(age))
    
     p<-ggplot(data=res, aes(x=year, y=M2, group=age)) +
       geom_line(aes(linetype=age,col=age))+
       geom_point(aes(shape=age,col=age))+
       facet_grid(species ~ ., scales="free_y")
       #ggtitle(unlist(data[1,'species']))
     print(p)
    ftable(round(xtabs(M2~species+year+age,data=filter(res,year %in% c(1974:1979,2020:2022))),6))
   }
   
   M2Graph(preys='COD') 
   M2Graph(preys=c('COD','WHG')) 
   M2Graph(preys=c('HER','NOP')) 
   M2Graph(preys=c('NSA','SSA')) 
 
   M2Graph(preys=c('HER')) 
   M2Graph(preys=c('NOP'))  # age 3 is zero ???????
}
