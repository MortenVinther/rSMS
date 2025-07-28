data$logSurveyObs
tail(data$keySurvey)

survey<- data.frame(cbind(data$keySurvey,logObs=data$logSurveyObs)) %>% 
  mutate(Year=y+abs(data$off.year), Age=factor(a-abs(data$off.age)),obs=exp(logObs),Fleet=paste(formatC(f,w=2,flag='0'),data$fleetNames[f]),Species=data$spNames[s])

out<-by(survey,survey$f,function(x){
  (xtabs(obs~Age+Year,data=droplevels(x)))
})

out[[1]]
out[[2]]
out[[3]]
out[[4]]

out<-by(survey,survey$s,function(x){
  x<-droplevels(x)
  #x<-x %>% group_by(Fleet) %>% mutate(obs=obs/mean(obs)) %>% ungroup()
  #print(x,n=300)
  ggplot(data=x, aes(x=Year, y=obs, shape=Age,col=Age)) +
    geom_line()+
    geom_point()+
    facet_wrap(vars(Fleet),ncol=2, scales="free_y")
})
out[[1]]
out[[2]]



# log obs
out<-by(survey,survey$s,function(x){
  x<-droplevels(x)
  ggplot(data=x, aes(x=Year, y=logObs, shape=Age,col=Age)) +
    geom_line()+
    geom_point()+
    facet_wrap(vars(Fleet),ncol=2, scales="free_y")
})
out[[1]]
out[[2]]
