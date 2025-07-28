into_annual<-function(inp) {

  data<-inp[['data']]
  data$nSeasons<-1L
  data$recSeason<-1L
  
  data$propMat <-lapply(data$propMat,function(x) x[,1,,drop=FALSE])
  data$stockMeanWeight<- lapply(data$stockMeanWeight,function(x) {
    x<-apply(x,c(1,3),mean) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$catchMeanWeight<- lapply(data$catchMeanWeight,function(x) x[,1,,drop=FALSE])  # to be changed
  data$seasFprop<- lapply(data$seasFprop,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$propF<- lapply(data$propF,function(x) x[,1,,drop=FALSE])   
  data$propM<- lapply(data$propM,function(x) x[,1,,drop=FALSE])
  data$natMor<- lapply(data$natMor,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$seasFprop<- lapply(data$seasFprop,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  
  x<-cbind(sampleTimeWithin=data$sampleTimeWithinSurvey,q=data$keySurvey.overview[,'q'])
  x[,"sampleTimeWithin"]<-x[,"sampleTimeWithin"]/4+(x[,"q"]-1)*0.25
  data$sampleTimeWithinSurvey<-x[,1]
  data$keySurvey[,'q']<-1L

  inp[['data']]<-data
  return(inp)
}




