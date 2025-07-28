plotProcessNResiduals<-function(runName, inclSp=1:100,fileLabel='N_residuals',outDir=data.path,mult=3, standardized=TRUE,
                             outFormat=c('screen','pdf','png')[1],longSpNames=TRUE,scaleY=c('fixed',"free_y")[1]) {
  

  cleanup()
  load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=F)
  
  # N process 
  
  a<-sms$rep$resAnnual
  b<-sms$data$doProcessNoise; b<-b[b<3]
  b<-data.frame(species=names(b),pOpt=b)
  b<-left_join(b,a,by = join_by(species)) %>% 
    filter(pOpt==3 & quarter==sms$data$recSeason & age==sms$data$minAge  |
           pOpt==3 & quarter==1 & age>sms$data$minAge |
           pOpt==2 & age==sms$data$minAge & quarter==sms$data$recSeason) %>%
    transmute(species,year,quarter,age,predN)
  
  aa<-extractParameters(sdrep=sms$sdrep, myMap=sms$map,data=sms$data,parameters=sms$parameters,lu=sms$lu)[[1]] %>% 
    filter(name=="logSdLogN" ) %>% transmute(species=Var1,age,variance=expEst^2)
  b<-left_join(aa,b,by = join_by(species, age))
  Resid<-left_join(b,sms$rep$res,by = join_by(species, year, quarter, age)) %>%
       transmute(species,s=match(species,sms$data$spNames),  year,quarter=factor(paste('Q:',quarter)),Age=factor(age),predict=log(predN),observed=log(N),Residual=log(predN/N), variance)
  Resid<-Resid %>% mutate(cols=if_else(Residual<0,'negativ','positiv')) %>% filter(s %in% inclSp)  


  
  plt<-by(Resid,Resid$s,function(x){
    if (longSpNames) sp<-sms$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
    if (standardized) {x$Residual<- x$Residual/sqrt(x$variance) ; stan<-'standardized'} else stan<-''
    maxRes<-max((abs(x$Residual)))
    pl<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = (abs(Residual)))) +
      geom_point(shape = 21,alpha=0.75) +
      scale_color_manual(values = c("blue", "red")) +
      scale_size(range = c(0,maxRes)*mult) +
      guides(color = guide_legend(override.aes = list(size = maxRes*mult/2)))+
      labs(x = "", y = "Age",title=paste0('Process N', ' ,',stan,sp)) +
      facet_wrap(vars(quarter),ncol=1,scales=scaleY)+
      my_theme() + 
      theme(legend.title = element_blank(),    legend.position = "right")
    
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = file.path(outDir,paste0(fileLabel,'_',sp,'.',outFormat)),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
  
  #shapiro.test(my_data$len)

  plt<-by(Resid,Resid$s,function(x){
    if (longSpNames) sp<-sms$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
    pl<-x %>% ggplot(aes(sample=Residual)) + stat_qq() + stat_qq_line()+
      labs(x = "Theoretical", y = "Sample quantiles",title=paste('Process N, ',sp)) +
      facet_wrap(vars(quarter),ncol=1,scales=scaleY)+
      my_theme() + 
      theme(legend.title = element_blank(),    legend.position = "right")
    
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = file.path(outDir,paste0(fileLabel,'_',sp,'.',outFormat)),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })

  
   plt<-by(Resid,Resid$s,function(x){
    if (longSpNames) sp<-sms$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
      pl<-x %>% ggplot(aes(x=Residual)) + 
     # geom_histogram( binwidth=0.15, fill="#69b3a2", color="#e9ecef", alpha=0.9)+ 
     geom_histogram( bins=16,fill="#69b3a2", color="#e9ecef", alpha=0.9)+ 
      labs(x = "Residual", y = "Count",title=paste('Process N, ',sp)) +
      facet_wrap(vars(quarter),ncol=1,scales=scaleY)+
      my_theme() + 
      theme(legend.title = element_blank(),    legend.position = "right")
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = file.path(outDir,paste0(fileLabel,'_',sp,'.',outFormat)),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
  
  
  invisible(Resid)
}

#runName<-"SMS_like_MS_recState"
qq<-plotProcessNResiduals(runName, inclSp=1:100,fileLabel='N_residuals',mult=10,standardized=FALSE,
                                outFormat=c('screen','pdf','png')[1],longSpNames=TRUE,scaleY=c('fixed',"free_y")[1]) 
  


