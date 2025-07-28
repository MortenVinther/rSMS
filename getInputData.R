getInputData<-function(inp=inp_all,mode=1){
  d<-inp[['data']]

  toDF1<-function(x,vname=c('year','quarter','age','M')) {
     do.call(rbind,lapply(seq_along(x), function(i){ 
        a<-array2DF(x[[i]]) 
        colnames(a)<-vname;
        a<-a %>% mutate(year=as.integer(year)-d$off.year,quarter=as.integer(quarter),age=as.integer(age)-d$off.age, species= names(x)[[i]])
      }))
  }

 wsea<-toDF1(x=d$stockMeanWeight,vname=c('year','quarter','age','wsea')) 
 M<-toDF1(x=d$natMor,vname=c('year','quarter','age','M')) 
 propMat<-toDF1(x=d$propMat,vname=c('year','quarter','age','propMat')) 
 canum<-toDF1(x=d$catchNumber,vname=c('year','quarter','age','canum')) 
 weca<-toDF1(x=d$catchMeanWeight,vname=c('year','quarter','age','weca')) 
 propM<-toDF1(x=d$ propM,vname=c('year','quarter','age','propM')) 
 propF<-toDF1(x=d$ propF,vname=c('year','quarter','age','propF')) 

 a<-left_join(wsea,M,by = join_by(year, quarter, age, species))
 a<-left_join(a,propMat,by = join_by(year, quarter, age, species))
 a<-left_join(a,canum,by = join_by(year, quarter, age, species))
 a<-left_join(a,weca,by = join_by(year, quarter, age, species))
 a<-left_join(a,propM,by = join_by(year, quarter, age, species))
 a<-left_join(a,propF,by = join_by(year, quarter, age, species)) %>% as_tibble()
 

if (mode==1) { 
 meanL<-toDF1(x=d$meanL,vname=c('year','quarter','age','meanL')) 
 consum<-toDF1(x=d$consum,vname=c('year','quarter','age','consum')) 
 natMor1<-toDF1(x=d$natMor1,vname=c('year','quarter','age','M1')) 
 propM2<-toDF1(x=d$propM2,vname=c('year','quarter','age','propM2')) 
 otherN<-toDF1(x=d$otherN,vname=c('year','quarter','age','otherN')) 
 
 b<-left_join(meanL,consum,by = join_by(year, quarter, age, species))
 b<-left_join(b,otherN,by = join_by(year, quarter, age, species))
 b<-left_join(b,natMor1,by = join_by(year, quarter, age, species))
 b<-left_join(b,propM2,by = join_by(year, quarter, age, species)) %>% as_tibble()
 b 

 sort(unique(b$species))
 a<-full_join(a,b,by = join_by(year, quarter, age, species)) %>% as_tibble()
}
 
 return(a)
}
