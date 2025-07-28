ls()[grep('batch',ls())]

b<-c("batch_default_configuration" ,
      "batch_final_single_configuration",
      "batch_seasonal_catch_configuration")

 b1<-batch_default_configuration(writeConrol=FALSE)
 b2<-batch_final_single_configuration(writeConrol=FALSE)  
 b3<-batch_seasonal_catch_configuration(writeConrol = FALSE)
 
 b<-list(b1=b1,b2=b2,b3=b3)
 n.<-slotNames(b1)
 theSame<-'the same'
 compare<-matrix(theSame,ncol=length(b),nrow=length(n.),dimnames=list(n.,names(b)))
 compare[,1]<-'  '

 for (i in (1:(length(b)-1))) {
   for (x in n.) {
     #cat(x,'\n')
     cla<-class(slot(b[[i]],x))[1]
     #cat(cl,'\n')
     if (cla=='list') {
       c1<-all(mapply(function(a,b)(length(a)==length(b)),slot(b[[i]],x),slot(b[[i+1]],x)))
       if (c1) c1<-all(unlist(mapply(function(a,b)(a==b),slot(b[[i]],x),slot(b[[i+1]],x)))) 
     } else c1<-all(slot(b[[i]],x)==slot(b[[i+1]],x))
     if (!c1) compare[x,i+1]<-'FALSE'
   }
 }
 
 if (!all(compare==theSame)) compare
 
 for (i in (1:(length(b)-1))) {
   for (x in n.) {
     if (compare[x,i+1]=='FALSE') {
       cat("###############  ",names(b)[i],'\n')
       print(slot(b[[i]],x))
       cat("##  ",names(b)[i+1],'\n')
       print(slot(b[[i+1]],x))
     }
  }}
 
 
