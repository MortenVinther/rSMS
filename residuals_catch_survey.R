plotCatchResiduals<-function(runName, inclSp=1:100,fileLabel='catch_residuals',mult=3,
                              outFormat=c('screen','pdf','png')[1], outDir=data.path,longSpNames=TRUE,standardized=FALSE,scaleY=c('fixed',"free_y")[1]) {
  cleanup()
  load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=F)
  
  Resid<-sms$rep$residCatch  %>% 
    transmute(s,species,year,Age=factor(age),quarter=factor(paste('Q:',q)),predict=logPred,observed=logObs,Residual=resid,variance)
  if (standardized) Resid$Residual<- Resid$Residual/sqrt(Resid$variance)
  Resid<-Resid %>% mutate(cols=if_else(Residual<0,'negativ','positiv')) %>% filter(s %in% inclSp)  
  
  plt<-by(Resid,Resid$s,function(x){
    if (longSpNames) sp<-sms$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
    maxRes<-max((abs(x$Residual)))
    pl<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = (abs(Residual)))) +
      geom_point(shape = 21,alpha=0.75) +
      scale_color_manual(values = c("blue", "red")) +
      scale_size(range = c(0,maxRes)*mult) +
      guides(color = guide_legend(override.aes = list(size = maxRes*mult/2)))+
      labs(x = "", y = "Age",title=paste(if_else(standardized,'Standardised',''),'Catch residuals, ',sp)) +
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
  invisible(NULL)
}

plotSurveyResiduals<-function(runName,inclFleets=1:300, inclSp=1:100, outFormat=c('screen','pdf','png')[1],outDir=data.path,fileLabel='survey_residuals',
                              mult=3,longSpNames=TRUE,standardized=TRUE,scaleY=c('fixed',"free_y")[1])  {
  
  cleanup()
  load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=F)
  
  Resid<-sms$rep$residSurv  %>% 
    transmute(s,species,year,Age=factor(age),quarter=q,f,fleet=fleetName,predict=logPred,observed=logObs,Residual=resid,variance)
    if (standardized) Resid$Residual<- Resid$Residual/sqrt(Resid$variance)
    Resid<-Resid %>% mutate(cols=if_else(Residual<0,'negativ','positiv')) %>% filter(f %in% inclFleets & s %in% inclSp)  
  
  plt<-by(Resid,Resid$s,function(x){
    if (longSpNames) sp<-sms$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
    maxRes<-max((abs(x$Residual)))
    pl<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = (abs(Residual)))) +
      geom_point(shape = 21,alpha=0.75) +
      scale_color_manual(values = c("blue", "red")) +
      scale_size(range = c(0,maxRes)*mult) +
      guides(color = guide_legend(override.aes = list(size = maxRes*mult/2)))+
      labs(x = "", y = "Age",title=paste(if_else(standardized,'Standardised',''),'Survey residuals, ',sp)) +
      facet_wrap(vars(fleet),ncol=1,scales=scaleY)+
      my_theme() + 
      theme(legend.title = element_blank(),    legend.position = "right")

    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename =file.path(outDir, paste0(fileLabel,'_',sp,'.',outFormat)),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
  invisible(NULL)
}
  
  
