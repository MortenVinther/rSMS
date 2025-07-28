
extractDataAndRun<-function(runName='Single',
                      my.ps=1,my.pso=0,
                      doMultiExtract=FALSE,
                      dir=data.path,
                      silent=TRUE,
                      sms.dat='rsms.dat',
                      fleetIn="new_fleet_info.dat",
                      smsConf=0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
                      ) 
{
    
  inp_all<-make_rsms_data(dir,sms.dat,fleetIn=fleetIn,multi=doMultiExtract) 
  inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
  
  inp$data$sms.mode<-smsConf
  
   #### prepare to run
  data<-inp[['data']]
  data$silent<-FALSE
  parameters<-inp[['parameters']]

  myMap<-map_param(data,parameters)
  random=c("Un","Uf")
  environment(func) <- environment() # see https://github.com/kaskr/RTMB/issues/26
  system.time(obj <- MakeADFun(func, parameters, random,silent=silent,map=myMap))
  
  lu<-lowerUpper(obj,data,parameters )
  opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000))
  announce(opt)
  
  sdrep <- sdreport(obj)
  cat('Hesssian:',sdrep$pdHess,'\n')
 
  myRep<-obj$report()
  saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
  
  
  plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[2],showSpecies=1:12,
                        inpRdata=list(runName,'SandeelR1'),
                        labels=c(runName,'ICES 2025'),
                        outFormat=c('screen','pdf','png')[3],
                        longSpNames=FALSE, fileLabel=runName)

  
  
  return(TRUE)
}

# to read the results from the ICES assessment
transExternalSummary(inp='ICES_SAN_R1_2025.out',outSet='SandeelR1',spNames='SA1',exSpeciesNames=c('SA1'))
transExternalData(inp='ICES_SAN_R1_2025_detailed_twoBlocks.out',outSet='SandeelR1_detailed_twoBlocks',spNames=data$spNames, exSpeciesNames=data$spNames)

#### prepare the "baseline run"
ps<-1
pso<-0

runName<-'RUN1'
sms.dat.in<-'rsmsBatchTemplate.dat'
fleetIn<-"new_fleet_info_template.dat"

system.time(res<-extractDataAndRun(runName,my.ps=ps,my.pso=pso,doMultiExtract=FALSE,sms.dat=sms.dat.in,fleetIn=fleetIn,silent=TRUE) )


load(file=file.path(data.path,paste0(runName,'.Rdata')),verbose=T)
a<-extractParameters(sms$sdrep,sms$map,sms$data)[[2]]
print(a,n=300)

#Used later on
off.age<-sms$data$off.age
off.year<-sms$data$off.year
nSpecies<-sms$data$nSpecies
speciesNames<-sms$data$spNames

###


### here we begin a new run
runName<-'RUN2' # the same as RUN 1 but without effort fleet in second season
sms.dat.out<-paste0(runName,'.dat')


rsms<-read.RSMS.control(dir=data.path,file=sms.dat.in)
# do something to the control
write.RSMS.control(rsms,file=file.path(data.path,sms.dat.out),write.multi=F,nice=T)

fleetOut<-paste0("new_fleet_info_",runName,".dat")
key<-readNewFleetInfo(dir=data.path,of=fleetIn,off.age=off.age,off.year=off.year,nSpecies=1,verbose=FALSE) 

# change fleet
key[['k']][,'useFleet']<-0L
key[['k']][c(1,2,3),'useFleet']<-1L
writeNewFleetInfo(dir=dir,of=fleetOut,key,speciesNames)

system.time(res<-extractDataAndRun(runName,my.ps=ps,my.pso=ps0,doMultiExtract=FALSE,sms.dat=sms.dat.out,fleetIn=fleetOut,silent=TRUE) )

### here we begin a new run
runName<-'RUN3' # recruits from parameters
sms.dat.out<-paste0(runName,'.dat')

rsms<-read.RSMS.control(dir=data.path,file=sms.dat.in)
rsms@incl.process.noise<-3L
rsms@SSB.R<-2L
# do something to the control
write.RSMS.control(rsms,file=file.path(data.path,sms.dat.out),write.multi=F,nice=T)

fleetOut<-paste0("new_fleet_info_",runName,".dat")
key<-readNewFleetInfo(dir=data.path,of=fleetIn,off.age=off.age,off.year=off.year,nSpecies=1,verbose=FALSE) 

# change fleet
#key[['k']][,'useFleet']<-0L
#key[['k']][c(1,2,3),'useFleet']<-1L
writeNewFleetInfo(dir=dir,of=fleetOut,key,speciesNames)

system.time(res<-extractDataAndRun(runName,my.ps=ps,my.pso=ps0,doMultiExtract=FALSE,sms.dat=sms.dat.out,fleetIn=fleetOut,silent=TRUE) )

load(file=file.path(data.path,paste0(runName,'.Rdata')),verbose=T)
a<-extractParameters(sms$sdrep,sms$map,sms$data)[[2]]
print(a,n=300)

#############################
runName<-'RUN6' # process noise all ages, annual catch in likelihood, random walk F
sms.dat.out<-paste0(runName,'.dat')

rsms<-read.RSMS.control(dir=data.path,file=sms.dat.in)
rsms@combined.catches<-1L
rsms@incl.process.noise<-1L
rsms@fModel=1L
rsms@keyVarLogN[[1]]<-c(0L,1L,3L)
rsms@keyLogFsta[[1]]<-c(0L,1L,2L)
rsms@catch.sep.year[[1]]<-1983L
rsms@catch.sep.age[[1]]<-0
rsms@SSB.R<-1L
rsms@use_rho<-1L

data$seasFprop[[1]][,1,]<-0.90
data$seasFprop[[1]][,2,]<-0.10
data$seasFprop[[1]][,2,1]<-1.0
data$logSeasFprop[[1]][,,]<-log(data$seasFprop[[1]][,,])

# do something to the control
write.RSMS.control(rsms,file=file.path(data.path,sms.dat.out),write.multi=F,nice=T)

fleetOut<-paste0("new_fleet_info_",runName,".dat")
key<-readNewFleetInfo(dir=data.path,of=fleetIn,off.age=off.age,off.year=off.year,nSpecies=1,verbose=FALSE) 

# change fleet
#key[['k']][,'useFleet']<-0L
#key[['k']][c(1,2,3),'useFleet']<-1L
writeNewFleetInfo(dir=dir,of=fleetOut,key,speciesNames)

system.time(res<-extractDataAndRun(runName,my.ps=ps,my.pso=ps0,doMultiExtract=FALSE,sms.dat=sms.dat.out,fleetIn=fleetOut,silent=TRUE) )

load(file=file.path(data.path,paste0(runName,'.Rdata')),verbose=T)
a<-extractParameters(sms$sdrep,sms$map,sms$data)[[2]]
print(a,n=300)
#################


inpRdata<- list('RUN1','RUN2','RUN3','RUN6')
labels<- c('RUN1','RUN2','RUN3','RUN6')


AICCompare(inpRdata,labels) 
  
plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[4],showSpecies=1:12,
                      inpRdata=list(runName,'SandeelR1_detailed'),
                      labels=c(runName,'ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san1R_')

plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[5],showSpecies=1:12,
                      inpRdata=list(runName,'SandeelR1_detailed'),
                      labels=c(runName,'ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      multN=0.000001,
                      longSpNames=FALSE, fileLabel='san1R_')

