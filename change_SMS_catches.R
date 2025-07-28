change_SMS_data<-function(fixCatch=c('NSA','SSA','NOP'),changes=NULL) {
  
  scanData<-function(name='canum_org.in',announce=FALSE){
    if (announce) cat('Reading ',name,'\n')
    x<-scan(file.path(data.path,name),comment.char='#',quiet=TRUE)
    stopifnot(tail(x,1)== -999)
    x<-head(x,-1)
  }
  CATCHN<-scanData('canum_org.in')
  WCATCH<-scanData('weca_org.in')
 
  b<-expand.grid(sub_area=1:SMS.control@no.areas,species.n=first.VPA:nsp,year=SMS.control@first.year:SMS.control@last.year,
                 quarter=1:SMS.control@last.season,age=SMS.control@first.age:SMS.control@max.age.all)
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  b<-data.frame(b,CATCHN=CATCHN,WCATCH=WCATCH) %>%  mutate(sp=SMS.control@species.names[species.n])
  
  fixCatchn <- match(fixCatch,SMS.control@species.names)
  
  if (length(fixCatchn)>0) {
    
    bb<-filter(b,species.n %in% fixCatchn)
    head(bb)
    
    crit1<-bb$sp=='NSA' & bb$age==0   # all 0 groups to quarter 3
    crit2<-bb$sp=='NSA' & bb$age>0   # all 1+ to quarter 2
    
    crit3<-bb$sp=='SSA' & bb$age==0   # all 0 groups to quarter 3
    crit4<-bb$sp=='SSA' & bb$age>0    # all 1+ to quarter 2
    
    crit5<-bb$sp=='NOP' & bb$age==0  & bb$quarter== 3  # all 0 groups in quarter 3 to quarter 4
    
    if (any(crit1)) bb[crit1,'quarter']<-3L
    if (any(crit2)) bb[crit2,'quarter']<-2L
    if (any(crit3)) bb[crit3,'quarter']<-3L
    if (any(crit4)) bb[crit4,'quarter']<-2L
    if (any(crit5)) bb[crit5,'quarter']<-4L
    
    bbb<-bb %>%   group_by(year, species.n, sp,quarter, sub_area, age) %>%
      summarize(WCATCH2=weighted.mean(WCATCH,w=CATCHN),CATCHN2=sum(CATCHN),.groups='drop') %>%
      mutate(WCATCH2=if_else(is.na(WCATCH2),0,WCATCH2)) 
    
    filter(bbb,year==1990 & CATCHN2>0 & species.n==22)
    filter(bbb,year==1990 & CATCHN2>0 & species.n==23)
    filter(bbb,year==1990 & CATCHN2>0 & species.n==24)
    
    bbb<-bbb %>% rename(CATCHN=CATCHN2,WCATCH=WCATCH2)
    b<-rbind(filter(b,!(species.n %in% fixCatchn)),bbb) 
  }
  
  
  if (!is.null(changes)) {
    a<-read.csv(changes) %>% mutate(remark=NULL)
    b<-left_join(b,a,by = join_by(sub_area, species.n, year, quarter, age)) %>%
      mutate(CATCHN=if_else(is.na(CATCHN_new),CATCHN,CATCHN_new),CATCHN_new=NULL)
  }
  
  CN<-tapply(b$CATCHN,list(b$year,b$quarter,b$species.n,b$age),sum)
  CW<-tapply(b$WCATCH,list(b$year,b$quarter,b$species.n,b$age),sum)
  
  CN[1,,'22',]
  CN[1,,'23',]
  CN[1,,'24',]
  
  transf<-function(item,file.name,dig) {
    out<-file.path(data.path,file.name)
    unlink(out)

    cat(paste("# Date ", date(), '\n'),file=out,append=TRUE)
    
    if (file.name=="canum.in") {
      cat(paste("#The file includes catch numbers at age",'\n'),file=out,append=TRUE)
      cat(paste("# if catch number=0  is applied for all ages in a year-season Zero fishing is assumed","\n"),file=out,append=TRUE)
      cat(paste("# if catch number<0  is catch numbers are assumened unknown (NA)","\n"),file=out,append=TRUE)
    }
    for (sp in dimnames(item)[[3]]) {
      out1<-item[,,sp,]
      for (y in dimnames(item)[[1]]) {
        out2<-out1[as.character(y),,]
        out2[is.na(out2)]<-0
        cat(paste("#",SMS.control@species.names[as.integer(sp)] ,"year:",y,"\n"),file=out,append=TRUE)
        write.table(format(out2,digits=dig),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
      }
    }
    cat(paste("-999 # checksum\n"),file=out,append=TRUE)
  }
  transf(CW,"weca.in",5)
  transf(CN,"canum.in",8)
}


change_SMS_data(fixCatch=c('NSA','SSA','NOP'),changes=file.path(data.path,'gem',"changes_CATCHN.csv")) 

  
