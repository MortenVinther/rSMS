
d<-inputToDF(data)
d$Age<-paste('Age',formatC(d$age,width=2,flag="0"))
d$quarter<-as.character(d$q)

d<-filter(d,data$seasonalCatches[s]==1)
p<-by(d,d$species,function(x) {
 ggplot(x, aes(x=year, y=seasFprop, fill=quarter ))+
  geom_bar(stat="identity")+theme_minimal()+ggtitle(x[1,'species'])+
  scale_fill_manual(values=c("lightblue","#999999", "#E69F00", "#56B4E9"))+
  facet_wrap(vars(Age),ncol=1)
})

p


out<-outputToDF(obj)

out2<-left_join(d,out,by = join_by(species, year, q, age)) %>% filter(west>0) 

ff<-FToDF(obj,sdrep) %>% filter(data$seasonalCatches[s]==1) %>%as_tibble()
seasFprop<- select(out2,species,s,year,q,age,seasFprop)
seasFprop
ff

ff<-left_join(ff,seasFprop,relationship = "many-to-many",by = join_by(s,year, age)) %>% mutate(FFq=FF*seasFprop)
ff

d<-left_join(out2,ff,by = join_by(s, seasFprop, species, year, q, age)) %>% 
   mutate(chat=FFq/Z*(N-1)*exp(-Z)) %>% group_by(s,year,age) %>% mutate(chatP=chat/sum(chat)) %>% ungroup()



p<-by(d,d$species,function(x) {
  x$quarter<-factor(x$q)
  # Add regression lines
  ggplot(x, aes(x=seasFprop, y=chatP, color=quarter, shape=quarter)) +
    geom_point() + 
    geom_smooth(method=lm)+
    facet_wrap(vars(Age),ncol=3)
 
})
p

p<-by(d,d$species,function(x) {
  x$quarter<-factor(x$q)
  # Add regression lines
  ggplot(x, aes(x=seasFprop, y=chatP/seasFprop, color=quarter, shape=quarter)) +
    geom_point() + 
    geom_smooth(method=lm)+
    facet_wrap(vars(Age),ncol=3)
  
})

p
p<-by(d,d$species,function(x) {
  ggplot(x, aes(x=year, y=seasFprop, fill=quarter ))+
    geom_bar(stat="identity")+theme_minimal()+ggtitle(x[1,'species'])+
    scale_fill_manual(values=c("lightblue","#999999", "#E69F00", "#56B4E9"))+
    facet_wrap(vars(Age),ncol=1)
})

p

  

d


survResid<-as.data.frame(data$keySurvey) %>% as_tibble()  %>% 
  mutate(predict=rep$predSurveyObs,obs=data$logSurveyObs,Residual=predict-obs,year=y-data$off.year,Age=factor(a-data$off.age),species=data$spNames[s]) %>% 
  left_join(.,overw,by = join_by(f, s)) %>% select(s,species,f,fleet,type,year,q,Age,s,predict,obs,Residual) %>%
  mutate(cols=if_else(Residual<0,'negativ','positiv')) 

plt<-by(survResid,survResid$s,function(x){
  plt<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=x[1,'species']) +
    # facet_grid(rows =vars(fleet), scales="free_y")
    facet_wrap(vars(fleet),ncol=1)+
    theme_bw() + 
    theme(strip.text = element_text(size = 10,margin = margin(0.01,0,0.01,0, "cm")))
  print(plt)
})
plt

