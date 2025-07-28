# the directory with our input data on list format
#list.data.path<-file.path(root,"data_northSea")


## See the end of this file for more options


##  input files  ###############



# ALK_stom_list.dat  optional wih addition of stomMark variable
# ALK_all_list.dat   optional wih addition of stomMark variable
# stomcon_list.dat   optional wih addition of stomMark variable

#########################################################################################

#trans.stomach        # transform stomach data from list input format to SMS format
#relevant if trans.stom = TRUE
#stom.first           # stomcon for inserted (first length class) values
#stom.mid             # stomcon for inserted (mid length class) values
#stom.last            # stomcon for inserted (last length class) values
#stom.exp            # stomcon for inserted ("expanded" length class) values
#stom.min.abs          # absolut minimum value for stom con for inserted values
#delete.tails         # Detelete tails (first and last "observations")
#inserted.haul.no.propor         # proportion of number hauls in case of invented values (first, mid and last values)
#max.other.food           # delete strata with more than max.other.food% oher food, use 100 for no deletion
#formatted.output              # make a nice table output (takes some minutes) or fast and unformatted

RSMS.stom.transform<-function( list.data.path=" ",stomMark=NULL,trans.bio=FALSE, trans.catch=FALSE,
                              trans.meanL=FALSE,  trans.meanL.from.weight=FALSE, trans.stockDist=FALSE,
                              trans.stomach=FALSE, trans.stomach.number=FALSE, use.stom.meanl.ALK=FALSE,
                              trans.ALK.stomach=FALSE, trans.ALK.all=FALSE, trans.other=FALSE, trans.Consum=FALSE,
                              stom.first=1E-06, stom.mid=1E-06, stom.last=1E-06, stom.exp=1E-06, stom.min.abs=1E-08, delete.tails=TRUE,
                              inserted.haul.no.propor=1.0,
                              max.other.food=100, sampling_effort=c("stom.no","haul.no")[2],
                              formatted.output=TRUE,selected.years=NULL,year.q=NULL,min.pred.length=0) {
  
  if (FALSE) {  #test
    list.data.path=file.path(data.path,'stom_input');
    stomMark="_Boots_0500_haul_as_observed";
    trans.bio=T; trans.catch=T;
    trans.meanL=T;  trans.meanL.from.weight=FALSE;trans.stockDist=F; trans.stomach.number<-FALSE
    trans.stomach=T; use.stom.meanl.ALK=FALSE; trans.ALK.stomach=T; trans.other=T; trans.Consum=T; trans.ALK.all=F;
    stom.first=1E-06; stom.mid=1E-06; stom.last=1E-06; stom.min.abs=1E-06; delete.tails=TRUE;
    inserted.haul.no.propor=1.0;
    formatted.output=T; selected.years=selected.years; year.q=year.q;min.pred.length=0; max.other.food=100
    stom.first=1E-06; stom.mid=1E-06; stom.last=1E-06; stom.exp=1E-06; stom.min.abs=1E-08; delete.tails=TRUE;
    inserted.haul.no.propor=1.0; sampling_effort=c("stom.no","haul.no")[1]
  
  }
  
  
  if (max.other.food!=100) stop("Sorry, a value different from 100 cannot be used for max.other.food")
  
  if (is.null(selected.years) & is.null(year.q)) stop("you have to specify selected.years or year.q.")
  if (!is.null(selected.years) & !is.null(year.q)) stop("you cannot specify both selected.years and year.q")
  
  species.no.prey<-data.frame(species=code.name,prey.no=0:(length(code.name)-1) )
  species.no.pred<-data.frame(species=code.name,pred.no=0:(length(code.name)-1) )
  
  maximum.age.all.species<-SMS.control@max.age.all
  years<-c(1,1)
  years[1]<-SMS.control@first.year
  years[2]<-SMS.control@last.year
  npr<-sum(SMS.control@species.info[,'predator'])
  no_areas<- SMS.control@no.areas
  area.names<-Read.area.names()
  #SMS_areas<-as.character(1:no_areas)
  SMS_areas<-(1:no_areas)
  
  checksum<-function(file='a'){
    cat("-999 # Checksum",file=file,append=TRUE)
  }
  
  select.ALK<-function(a){   #function for selection of ALK data
    if (!('ALK.MeanL' %in% names(a))) a$ALK.MeanL<-0
    if (is.null(year.q)) {
      b<-subset(a, (year %in% selected.years)   & (prey %in% code.name.prey),
                select=c(SMS_area,year,quarter,prey, prey.age,prey.size.class,prey.size,ALK.MeanL,ALK))
    } else {
      a$year.qq<-paste(a$year,'q',a$quarter,sep='')
      b<-subset(a, (year.qq %in% year.q)   & (prey %in% code.name.prey),
                select=c(SMS_area,year,quarter,prey, prey.age,prey.size.class,prey.size,ALK.MeanL,ALK))
    }
    merge(b,species.no.prey,by.x="prey",by.y="species")
  }
  
  
  write_stom_incl<-function(tot) {
    
    n<-subset(tot, first.year.quarter.area.pred.predL)%>%dplyr::select(year,quarter,pred,pred.size.class,pred.no,samp.eff,samp.eff.scaled)
    n$pred.no<-paste(formatC(n$pred.no,wid=2,flag='0'),n$pred,sep='_')
    
    aa<-tapply(n$samp.eff,list(n$pred.size.class,n$year,n$quarter,n$pred.no),sum)
    aa[is.na(aa)]<-0
    aadim<-dimnames(aa)
    ofile<-'incl_stom.in'
    cat('# file incl_stom.in. Date ', date(),' \n# A value >=1 indicate that the stomach sample is used in SMS. Default number is no. of hauls (a number > 0)
  \n##############\n',file=ofile)
    
    for (s in aadim[[4]] ){
      cat('##############\n',file=ofile,append=TRUE)
      for (q in aadim[[3]]) {
        cat("# Predator: ",s,'\n',file=ofile,append=TRUE)
        cat('# quarter',q,'\n',file=ofile,append=TRUE)
        cat('#',format(aadim[[2]],width=5),'\n',file=ofile,append=TRUE)
        write.table(format(round(aa[,,q,s],0),width=5),file=ofile,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
      }
    }
    cat('-999 # checksum\n',file=ofile,append=TRUE)
  }
  
  select.ALK.all<-function(a){   #function for selection of ALK data
    if (is.null(year.q)) {
      b<-subset(a,  (prey %in% code.name),
                select=c(SMS_area,year,quarter,prey, prey.age,prey.size.class,prey.size,ALK.MeanL,ALK))
    } else {
      a$year.qq<-paste(a$year,'q',a$quarter,sep='')
      b<-subset(a, (prey %in% code.name),
                select=c(SMS_area,year,quarter,prey, prey.age,prey.size.class,prey.size,ALK.MeanL,ALK))
    }
    merge(b,species.no.prey,by.x="prey",by.y="species")
  }
  
  #######
  
  select.stom<-function(a){   #function for selection of stomach data
    if (is.null(year.q)) {
      if (!trans.stomach.number)  { a$used.prey.number<-0; a$calc.prey.number<-0 }
      b<-subset(a,(year %in% selected.years)  & (pred %in% code.name.pred)
                & (haul.no>=min.stomach.sampled.in.stratum) & (pred.mean.length>=min.pred.length)
                & ((type %in% c('first','last','exp') & !delete.tails) | type %in% c('obs','mid')),
                select=c(SMS_area,year,quarter,pred,pred.size,pred.size.class,pred.mean.length,phi,
                         prey,prey.size,prey.size.class,prey.mean.length,type,
                         stomcon,mean.weight,samp.eff,samp.eff.scaled, calc.prey.number,used.prey.number))
    } else {
      a$year.qq<-paste(a$year,'q',a$quarter,sep='')
      b<-subset(a,(year.qq %in% year.q)  & (pred %in% code.name.pred)
                & (haul.no>=min.stomach.sampled.in.stratum)
                & ((type %in% c('first','last','exp') & !delete.tails) | type %in% c('obs','mid')),
                select=c(SMS_area,year,quarter,pred,pred.size,pred.size.class,pred.mean.length,phi,
                         prey,prey.size,prey.size.class,prey.mean.length,type,
                         stomcon,mean.weight,samp.eff,samp.eff.scaled,,calc.prey.number,used.prey.number))
    }
    
    #cat('\ntest1\n');print(head(b))
    
    non.prey<-subset(b,!(prey %in% code.name.prey),
                     select=c(SMS_area,year,quarter,pred,pred.size,prey,stomcon))
    
    if (dim(non.prey)[1]>0) {
      print(paste(unique(non.prey$prey)," to other food"))
      b<-subset(b,(prey %in% code.name.prey))
      non.prey[,'prey']<-'OTH'
      oth<-tapply(non.prey$stomcon,list(non.prey$SMS_area,non.prey$year,non.prey$quarter,non.prey$pred,non.prey$pred.size,non.prey$prey),sum)
      a<-subset(arr2df(oth),stomobs>0)
      names(a)<-list('SMS_area','year','quarter','pred','pred.size','prey','stomcon.new')
      c<-merge(b,a,all.x=TRUE)
      c[!is.na(c$stomcon.new),'stomcon']<-c[!is.na(c$stomcon.new),'stomcon']+c[!is.na(c$stomcon.new),'stomcon.new']
      b<-subset(c,select=-stomcon.new)
      
    }
    
    b<-merge(b,species.no.pred,by.x="pred",by.y="species")
    b<-merge(b,species.no.prey,by.x="prey",by.y="species")
    
    if (max.other.food<100) {
      
      # calculate relative stomach contents,
      a<-by(b,list(a$SMS_area,b$year,b$quarter,b$pred,b$pred.size.class),function(x) sum(x$stomcon))
      a<-subset(arr2df(a),stomobs>0)
      names(a)<-list('SMS_area','year','quarter','pred','pred.size.class','sum')
      
      b<-merge(b,a)
      b$stomcon<-b$stomcon/b$sum;
      
      
      # extract strata with more than x% other food
      a<-subset(b,prey.no==0 & stomcon>(max.other.food/100))
      a<-data.frame(year=a$year,quarter=a$quarter,pred.no=a$pred.no,pred.size.class=a$pred.size.class,del=TRUE)
      
      b<-merge(b,a,all=TRUE)
      cat("\nDeleted strata:\n")
      print(a)
      b<-b[is.na(b$del),]
    }
    #cat('\ntest2\n');print(head(b))
    b
  }
  
  
  #######
  
  
  
  #################################################################
 
  
  
  trans.4M.SMS.stomach<-function(){
    
    cat("trans.4M.SMS.stomach\n")
    if (is.null(stomMark)) filen<-'ALK_stom_list.dat' else  filen<-paste0('ALK_stom_list',stomMark,'.dat')
    file<-file.path(list.data.path,filen)
    cat(file,'\n')
    lak<-read.table(file,header=TRUE)
    used_ages<-data.frame(prey=tail(code.name,-1),plusage=SMS.control@species.info[,'last-age'])
    lak<-full_join(lak,used_ages,by = join_by(prey)) %>% filter(prey.age<=plusage)
    #filter(lak,prey=='WHG')
    lak<-select.ALK(lak)
    
    lak$ALK<-lak$ALK/100
    
    if (is.null(stomMark)) filen<-'stomcon_list.dat' else  filen<-paste0('stomcon_list',stomMark,'.dat')
    file<-file.path(list.data.path,filen)
    cat(file,'\n')
    stom<-read.table(file,header=TRUE,na.strings=".")
    if (any(is.na(stom))) stop(paste("program halted: file ", filen," includes missing data"))
    
    
    if (use.stom.meanl.ALK) stom$prey.mean.length<-stom$prey.mean.length.ALK
    
    if (is.null(var.groups.size)) stom$var.groups<-stom$pred.size else stom$var.groups<-var.groups[match(stom$pred.size,var.groups.size)]
    
    
    if (sampling_effort=="haul_no") stom$samp.eff<-stom$haul.no else stom$samp.eff<-stom$stom.no
    
    a<-aggregate(samp.eff~SMS_area+pred+var.groups,data=stom,FUN=max)
    b<-names(a)
    names(a)<- c("SMS_area","pred",    "var.groups",   "max.N" )
    stom<-merge(stom,a)
    stom$samp.eff.scaled<-round(stom$samp.eff/stom$max.N*100)
    stom[stom$samp.eff.scaled<min.stom.groups,"samp.eff.scaled"]<-min.stom.groups
    
    stom<-select.stom(stom)
    #print(summary(stom))
    # print(subset(stom,is.na(stomcon)))
    
    
    # adjust invented observations (first,mid and last)
    stom[stom$type=='first','stomcon']<-stom.first
    stom[stom$type=='mid','stomcon']<-stom.mid
    stom[stom$type=='last','stomcon']<-stom.last
    stom[stom$type=='exp','stomcon']<-stom.exp
    
    stom[stom$type!='obs' & stom$stomcon<stom.min.abs,'stomcon']<-stom.min.abs
    
    stom[stom$type!='obs','samp.eff']<-stom[stom$type!='obs','samp.eff']*inserted.haul.no.propor
    stom[stom$type!='obs','samp.eff.scaled']<-stom[stom$type!='obs','samp.eff.scaled']*inserted.haul.no.propor
    
    stom$type.no<-0
    stom[stom$type=='obs','type.no']<-1
    stom[stom$type=='mid','type.no']<-2
    stom[stom$type=='first','type.no']<-3
    stom[stom$type=='last','type.no']<-3
    stom[stom$type=='exp','type.no']<-4
    
    if (trans.stomach) {
      sizes<-rbind(
        select(stom,year,quarter,species=prey,species.no=prey.no,size=prey.size,size.class=prey.size.class),
        select(stom,year,quarter,species=pred,species.no=pred.no,size=pred.size,size.class=pred.size.class),
        select(lak ,year,quarter,species=prey,species.no=prey.no,size=prey.size,size.class=prey.size.class)
      ) %>% unique() %>% arrange(year,quarter,species.no,size,size.class) %>% as_tibble()
      out<-file.path(data.path,'stom_size_classes.in')
      unlink(out)
      write.table(sizes,file=out,append=FALSE, quote = TRUE, row.names = FALSE,col.names = TRUE,sep=',')
    }
    
    # recalculate to relative stomach contents, summing up to 1
    a<-by(stom,list(stom$SMS_area,stom$year,stom$quarter,stom$pred,stom$pred.size.class),function(x) sum(x$stomcon))
    a<-subset(arr2df(a),stomobs>0)
    names(a)<-list('SMS_area','year','quarter','pred','pred.size.class','sum')
    
    stom<-merge(stom,a)
    stom$stomcon<-stom$stomcon/stom$sum;
    #cat('\ntest3\n');print(head(stom))
    
    #check existence of ALK keys for all prey sizes
    a<-unique(data.frame(SMS_area=lak$SMS_area,year=lak$year,quarter=lak$quarter,species=lak$prey,size=lak$prey.size.class,a="a"))
    # cat('\nLAK\n');print(head(a,20)) # remove
    
    b2<-unique(data.frame(MS_area=stom$SMS_area,year=stom$year,quarter=stom$quarter,species=stom$prey,size=stom$prey.size.class))
    b3<-b2
    b3<-unique(subset(b3,species!='OTH'))
    b3$key<-paste(b3$SMS_area,b3$year,b3$quarter,b3$species,formatC(b3$size,wid=2,flag="0"))
    s<-order(b3$key)
    b4<-b3[s,]
    b4<-subset(b4,select=-key)
    b4$b<-"b"
    # cat('\nSTOM\n'); print(head(b4,20))  # remove
    check<-merge(a,b4,all=TRUE)
    
    mis.alk<-subset(check,is.na(a))
    if (dim(mis.alk)[1]>0) {
      print(mis.alk)
      print(check)
      stop('program stop')
    } else cat("Data are OK\n")
    
    
    ###############################################################################
    
    
    tot<-stom
    tot<- tot %>% arrange(SMS_area,year,quarter,pred.no,pred.size.class,prey.no,prey.size.class) %>% as_tibble() %>% rename(area=SMS_area)
    
    #test
    #tot<-select(tot, -samp.eff, -samp.eff.scaled, -calc.prey.number,-pred.size, -pred.mean.length,-prey.mean.length,-used.prey.number,-prey.size,-pred.size,-phi)
    nobs<-dim(tot)[[1]]
    tot$obs_no<- 1:nobs
  
    s<-    tot %>% group_by(area,year,quarter,pred.no,pred.size.class,prey.no) %>% mutate(f_prey.no=if_else(!duplicated(prey.no),obs_no,NA)) %>%
                group_by(area,year,quarter,pred.no,pred.size.class) %>% mutate(f_pred.size.class=if_else(!duplicated(pred.size.class),obs_no,NA)) %>%
                group_by(area,year,quarter,pred.no) %>% mutate(f_pred.no=if_else(!duplicated(pred.no),obs_no,NA)) %>%
                group_by(area,year,quarter) %>% mutate(f_quarter=if_else(!duplicated(quarter),obs_no,NA)) %>% 
                group_by(area,year) %>% mutate(f_year=if_else(!duplicated(year),obs_no,NA)) %>% 
                group_by(area) %>% mutate(f_area=if_else(!duplicated(area),obs_no,NA)) %>%
                ungroup()
 
   s %>%print(n=20)
 
   key<-filter(s, !is.na(f_area)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,area);  key[dim(key)[1],'last']<-nobs
   keySA<-key
   

   key<-filter(s, !is.na(f_year)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,year);  key[dim(key)[1],'last']<-nobs
   keySAY<-key
 
   key<-filter(s, !is.na(f_quarter)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,year);  key[dim(key)[1],'last']<-nobs
   keySAYQ<-key
 
   key<-filter(s, !is.na(f_pred.no)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,pred.no,pred);  key[dim(key)[1],'last']<-nobs
   keySAYQP<-key
 
   key<-filter(s, !is.na(f_pred.size.class)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,pred.size.class,pred.size,pred.mean.length,samp.eff,samp.eff.scaled,phi);  key[dim(key)[1],'last']<-nobs
   keySAYQPL<-key
   
   key<-filter(s, !is.na(f_prey.no)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,prey.no,prey);  key[dim(key)[1],'last']<-nobs
   keySAYQPLP<-key
   
   key<-filter(s, !is.na(f_prey.no)) %>% mutate(last=lead(obs_no)-1L) %>% transmute(first=obs_no,last,prey.no,prey);  key[dim(key)[1],'last']<-nobs
   keySAYQPLP<-key
   
   stomObs<- s %>% transmute(obs_no,stomcon,meanWeight=mean.weight,preySizeClass=prey.size.class,preyMeanL=prey.mean.length,type=type.no) 
  
   
   
   
    ############## Print ALK ##############
    out<-file.path(data.path,'alk_stom.in')
    unlink(out)
    cat(line,file=out,append=TRUE)
    cat(paste("# File for splitting age-group on length classess (ALK). Used for likelihood for stomach data\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    cat(paste("# n_ALK_y:  number of years with ALK\n", length(unique(tot$year))," \n"),file=out,append=TRUE)
    cat(paste("# n_ALK_yq:  number of year, quarter combinations with ALK \n",length(unique(tot$year.quarter)),"\n"),file=out,append=TRUE)
    cat(paste("# n_ALK_yqd:  number of year, quarter area combinations with ALK \n",length(unique(tot$year.quarter.area)),"\n"),file=out,append=TRUE)
    cat(paste("# n_ALK_yqds:  number of year, quarter area species combinations with ALK \n"),file=out,append=TRUE)
    cat(length(unique(tot$year.quarter.area.prey)),file=out,append=TRUE)
    cat(paste("\n# n_ALK_yqdsa:  number of year, quarter, area species, age combinations with ALK \n"),file=out,append=TRUE)
    cat(paste(length(unique(tot$year.quarter.area.prey.age)),"\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    
    cat(paste("# Year name, and first and last year-quarter index\n# ALK_y_name if_ALK_y il_ALK_y\n"),file=out,append=TRUE)
    
    aa<-tot[!duplicated(tot$year.quarter),]
    aa$index<-seq(1,dim(aa)[1])
    minn<-tapply(aa$index,list(aa$year),min)
    maxx<-tapply(aa$index,list(aa$year),max)
    cat(paste("# stl_y: year name, first and last year-quarter index\n"),file=out,append=TRUE)
    yy<-as.numeric(unlist(dimnames(minn)))
    for (i in (1:dim(minn))) cat(paste(yy[i],minn[i],maxx[i],"\n"),file=out,append=TRUE)
    
    aa<-tot[!duplicated(tot$year.quarter.area),]
    aa$index<-seq(1,dim(aa)[1])
    aa<-subset(aa,select=c(index,year,quarter,SMS_area))
    minn<-tapply(aa$index,list(aa$year,aa$quarter),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter),max)
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    cat(paste("# stl_yq: quarter name, first and last quarter-area index\n"),file=out,append=TRUE)
    for (j in (1:length(yy))) for(i in (1:length(qq))) if (!is.na(minn[j,i]) & !is.na(maxx[j,i])) cat(paste(qq[i],minn[j,i],maxx[j,i],"\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    
    cat(paste("# stl_yqd: area name, and first and last year-quarter-species-age index\n"),file=out,append=TRUE)
    aa<-tot[!duplicated(tot$year.quarter.area.prey),]
    aa$index<-seq(1,dim(aa)[1])
    minn<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area),max)
    
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    pp<-as.numeric(unlist(dimnames(minn)[3]))
    
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(pp))) {
      if (!is.na(minn[j,i,k]) & !is.na(maxx[j,i,k])) cat(paste(pp[k],minn[j,i,k],maxx[j,i,k],"\n"),file=out,append=TRUE)
    }
    
    cat(line,file=out,append=TRUE)
    cat(paste("# stl_yqdp: species name, and first and last area-species index\n"),file=out,append=TRUE)
    aa<-tot[!duplicated(tot$year.quarter.area.prey.age),]
    aa$index<-seq(1,dim(aa)[1])
    aa<-subset(aa,select=c(index,year,quarter,SMS_area,prey.no))
    minn<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area,aa$prey.no),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area,aa$prey.no),max)
    
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    dd<-as.numeric(unlist(dimnames(minn)[3]))
    pp<-as.numeric(unlist(dimnames(minn)[4]))
    
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(dd))) for(l in (1:length(pp))) {
      if (!is.na(minn[j,i,k,l]) & !is.na(maxx[j,i,k,l])) cat(paste(pp[l],minn[j,i,k,l],maxx[j,i,k,l],"\n"),file=out,append=TRUE)
    }
    
    cat(line,file=out,append=TRUE)
    cat(paste("# species age name, and first and last species-age index\n"),file=out,append=TRUE)
    cat(paste("# ALK_yqasa_name if_ALK_yqasa il_ALK_yqasa\n"),file=out,append=TRUE)
    minn<-tapply(tot$prey.size.class,list(tot$year,tot$quarter,tot$SMS_area,tot$prey.no,tot$prey.age),min)
    maxx<-tapply(tot$prey.size.class,list(tot$year,tot$quarter,tot$SMS_area,tot$prey.no,tot$prey.age),max)
    
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    dd<-as.numeric(unlist(dimnames(minn)[3]))
    pp<-as.numeric(unlist(dimnames(minn)[4]))
    pa<-as.numeric(unlist(dimnames(minn)[5]))
    
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(dd))) for(l in (1:length(pp))) for(m in (1:length(pa))){
      if (!is.na(minn[j,i,k,l,m]) & !is.na(maxx[j,i,k,l,m])) cat(paste(pa[m],minn[j,i,k,l,m],maxx[j,i,k,l,m],"\n"),file=out,append=TRUE)
    }
    
    
    cat(line,file=out,append=TRUE)
    cat(paste("# ALK by year, quarter, species, age and length\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    if (formatted.output) {
      for (i in (1:dim(tot)[1])) {
        if (tot$first.year.quarter[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# year:",tot$year[i],"Quarter:",tot$quarter[i],'\n'),file=out,append=TRUE)
          cat(paste('ALK output, ALKS year:',tot$year[i]," Quarter:",tot$quarter[i],'\n'))
        }
        if (tot$first.year.quarter.area[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# SMS_area:",tot$SMS_area[i],'\n'),file=out,append=TRUE)
          cat(line,file=out,append=TRUE)
        }
        
        if (tot$first.year.quarter.area.prey[i]) cat(paste("# species:",tot$prey[i],'\n'),file=out,append=TRUE)
        if (tot$first.year.quarter.area.prey.age[i]) cat(paste("# age:",tot$prey.age[i],"first size:",tot$prey.size.class[i],tot$prey.size[i],'\n'),file=out,append=TRUE)
        cat(formatC(tot$ALK[i],format="f",dig=5,width=8),file=out,append=TRUE)
        if (tot$last.age[i]) cat('\n',file=out,append=TRUE)
      }
    } else {
      cat(formatC(tot$ALK,format="f",dig=5,width=8),file=out,append=TRUE)
    }
    
    
    ####
    # print mean length at size class
    
    cat(line,file=out,append=TRUE)
    cat(paste("# mean length (mm) by year, quarter, species, age and length\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    if (formatted.output) {
      for (i in (1:dim(tot)[1])) {
        if (tot$first.year.quarter[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# year:",tot$year[i],"Quarter:",tot$quarter[i],'\n'),file=out,append=TRUE)
          cat(paste('ALK-mean length output, ALKS year:',tot$year[i]," Quarter:",tot$quarter[i],'\n'))
        }
        if (tot$first.year.quarter.area[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# SMS_area:",tot$SMS_area[i],'\n'),file=out,append=TRUE)
          cat(line,file=out,append=TRUE)
        }
        
        if (tot$first.year.quarter.area.prey[i]) cat(paste("# species:",tot$prey[i],'\n'),file=out,append=TRUE)
        if (tot$first.year.quarter.area.prey.age[i]) cat(paste("# age:",tot$prey.age[i],"first size:",tot$prey.size.class[i],tot$prey.size[i],'\n'),file=out,append=TRUE)
        cat(formatC(tot$ALK.MeanL[i],format="f",dig=0,width=6),file=out,append=TRUE)
        if (tot$last.age[i]) cat('\n',file=out,append=TRUE)
      }
    } else {
      cat(formatC(tot$ALK.MeanL,format="f",dig=0,width=6),file=out,append=TRUE)
    }
    cat(paste("#\n -999.0 # Check sum \n"),file=out,append=TRUE)
    
  }
  
  
  
  trans.4M.SMS.ALK.all<-function(){
    cat("trans.4M.SMS.ALK.all\n")
    if (is.null(stomMark)) filen<-'ALK_all_list.dat' else  filen<-paste0('ALK_all_list',stomMark,'.dat')
    file<-file.path(list.data.path,filen)
    lak<-read.table(file,header=TRUE)
    
    lak<-select.ALK.all(lak)
    
    lak$ALK=lak$ALK/100
    
    tot<-lak
    tot$year.quarter<-                   paste(tot$year,tot$quarter)
    tot$year.quarter.area<-              paste(tot$year.quarter,formatC(tot$SMS_area,wid = 2, flag = "0"))
    tot$year.quarter.area.prey<-         paste(tot$year.quarter.area,formatC(tot$prey.no,wid = 2, flag = "0"))
    tot$year.quarter.area.prey.age<-     paste(tot$year.quarter.area.prey,formatC(tot$prey.age,wid = 2, flag = "0"))
    tot$year.quarter.area.prey.age.size<-paste(tot$year.quarter.area.prey.age,formatC(tot$prey.size.class,wid = 2, flag = "0"))
    
    key<-order(tot$year.quarter.area.prey.age.size)
    tot<-tot[key,]
    #print(dim(tot))
    #print(tot)
    
    tot$first.year.quarter<-!duplicated(tot$year.quarter)
    tot$first.year.quarter.area<-!duplicated(tot$year.quarter.area)
    tot$first.year.quarter.area.prey<-!duplicated(tot$year.quarter.area.prey)
    tot$first.year.quarter.area.prey.age<-!duplicated(tot$year.quarter.area.prey.age)
    tot$last.age<-c(tot$first.year.quarter.area.prey.age[2:dim(tot)[1]],TRUE)
    
    print(sum(tot$first.year.quarter.area.prey.age))
    
    line<-'##############################################################\n'
    
    ############## Print ALKS ##############
    out<-file.path(data.path,'alk_all.in')
    unlink(out)
    cat(line,file=out,append=TRUE)
    cat(paste("# File for splitting age-group on length classess (ALK). Used for calculation of M2 values\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    cat(paste("# n_ALKS_y:  number of years with ALK\n", length(unique(tot$year))," \n"),file=out,append=TRUE)
    cat(paste("# n_ALKS_yq:  number of year, quarter combinations with ALK \n",length(unique(tot$year.quarter)),"\n"),file=out,append=TRUE)
    cat(paste("# n_ALKS_yqa:  number of year, quarter, area combinations with ALK \n",length(unique(tot$year.quarter.area)),"\n"),file=out,append=TRUE)
    cat(paste("# n_ALKS_yqas:  number of year, quarter, area, species combinations with ALK \n"),file=out,append=TRUE)
    cat(length(unique(tot$year.quarter.area.prey)),file=out,append=TRUE)
    cat(paste("\n# n_ALKS_yqasa:  number of year, quarter, area, species, age combinations with ALK \n"),file=out,append=TRUE)
    cat(paste(length(unique(tot$year.quarter.area.prey.age)),"\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    
    cat(paste("# Year name, and first and last year-quarter index\n# ALKS_y_name if_ALKS_y il_ALKS_y\n"),file=out,append=TRUE)
    
    aa<-tot[!duplicated(tot$year.quarter),]
    aa$index<-seq(1,dim(aa)[1])
    minn<-tapply(aa$index,list(aa$year),min)
    maxx<-tapply(aa$index,list(aa$year),max)
    cat(paste("# stl_y: year name, first and last year-quarter index\n"),file=out,append=TRUE)
    yy<-as.numeric(unlist(dimnames(minn)))
    for (i in (1:dim(minn))) cat(paste(yy[i],minn[i],maxx[i],"\n"),file=out,append=TRUE)
    
    
    aa<-tot[!duplicated(tot$year.quarter.area),]
    aa$index<-seq(1,dim(aa)[1])
    aa<-subset(aa,select=c(index,year,quarter,SMS_area))
    minn<-tapply(aa$index,list(aa$year,aa$quarter),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter),max)
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    cat(paste("# stl_yq: quarter name, first and last quarter-area index\n"),file=out,append=TRUE)
    for (j in (1:length(yy))) for(i in (1:length(qq))) if (!is.na(minn[j,i]) & !is.na(maxx[j,i])) cat(paste(qq[i],minn[j,i],maxx[j,i],"\n"),file=out,append=TRUE)
    
    
    aa<-tot[!duplicated(tot$year.quarter.area.prey),]
    aa$index<-seq(1,dim(aa)[1])
    aa<-subset(aa,select=c(index,year,quarter,SMS_area))
    minn<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area),max)
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    dd<-as.numeric(unlist(dimnames(minn)[3]))
    cat(paste("# stl_yqa: area name, first and last area-species index\n"),file=out,append=TRUE)
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(dd))) {
      if (!is.na(minn[j,i,k]) & !is.na(maxx[j,i,k])) cat(paste(dd[k],minn[j,i,k],maxx[j,i,k],"\n"),file=out,append=TRUE)
    }
    
    cat(line,file=out,append=TRUE)
    cat(paste("#  species name, and first and last area-species index\n"),file=out,append=TRUE)
    aa<-tot[!duplicated(tot$year.quarter.area.prey.age),]
    aa$index<-seq(1,dim(aa)[1])
    aa<-subset(aa,select=c(index,year,quarter,SMS_area,prey.no))
    minn<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area,aa$prey.no),min)
    maxx<-tapply(aa$index,list(aa$year,aa$quarter,aa$SMS_area,aa$prey.no),max)
    
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    dd<-as.numeric(unlist(dimnames(minn)[3]))
    pp<-as.numeric(unlist(dimnames(minn)[4]))
    
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(dd))) for(l in (1:length(pp))) {
      if (!is.na(minn[j,i,k,l]) & !is.na(maxx[j,i,k,l])) cat(paste(pp[l],minn[j,i,k,l],maxx[j,i,k,l],"\n"),file=out,append=TRUE)
    }
    
    cat(line,file=out,append=TRUE)
    cat(paste("# species age name, and first and last species-age index\n"),file=out,append=TRUE)
    cat(paste("# ALKS_yqasa_name if_ALKS_yqasa il_ALKS_yqasa\n"),file=out,append=TRUE)
    minn<-tapply(tot$prey.size.class,list(tot$year,tot$quarter,tot$SMS_area,tot$prey.no,tot$prey.age),min)
    maxx<-tapply(tot$prey.size.class,list(tot$year,tot$quarter,tot$SMS_area,tot$prey.no,tot$prey.age),max)
    
    yy<-as.numeric(unlist(dimnames(minn)[1]))
    qq<-as.numeric(unlist(dimnames(minn)[2]))
    dd<-as.numeric(unlist(dimnames(minn)[3]))
    pp<-as.numeric(unlist(dimnames(minn)[4]))
    pa<-as.numeric(unlist(dimnames(minn)[5]))
    
    for (j in (1:length(yy))) for(i in (1:length(qq))) for(k in (1:length(dd))) for(l in (1:length(pp))) for(m in (1:length(pa))){
      if (!is.na(minn[j,i,k,l,m]) & !is.na(maxx[j,i,k,l,m])) cat(paste(pa[m],minn[j,i,k,l,m],maxx[j,i,k,l,m],"\n"),file=out,append=TRUE)
    }
    
    
    cat(line,file=out,append=TRUE)
    cat(paste("# ALK by year, quarter, species, age and length\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    if (formatted.output) {
      for (i in (1:dim(tot)[1])) {
        if (tot$first.year.quarter[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# year:",tot$year[i],"Quarter:",tot$quarter[i],'\n'),file=out,append=TRUE)
          cat(paste('ALK output, ALKS year:',tot$year[i]," Quarter:",tot$quarter[i],'\n'))
        }
        if (tot$first.year.quarter.area[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# SMS_area:",tot$SMS_area[i],'\n'),file=out,append=TRUE)
          cat(line,file=out,append=TRUE)
        }
        
        if (tot$first.year.quarter.area.prey[i]) cat(paste("# species:",tot$prey[i],'\n'),file=out,append=TRUE)
        if (tot$first.year.quarter.area.prey.age[i]) cat(paste("# age:",tot$prey.age[i],"first size:",tot$prey.size.class[i],tot$prey.size[i],'\n'),file=out,append=TRUE)
        cat(formatC(tot$ALK[i],format="f",dig=5,width=8),file=out,append=TRUE)
        if (tot$last.age[i]) cat('\n',file=out,append=TRUE)
      }
    } else {
      cat(formatC(tot$ALK,format="f",dig=5,width=8),file=out,append=TRUE)
    }
    ####
    # print mean length at size class
    
    cat(line,file=out,append=TRUE)
    cat(paste("# Size (mm) by year, quarter, species, age and length\n"),file=out,append=TRUE)
    cat(line,file=out,append=TRUE)
    if (formatted.output) {
      for (i in (1:dim(tot)[1])) {
        if (tot$first.year.quarter[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# year:",tot$year[i],"Quarter:",tot$quarter[i],'\n'),file=out,append=TRUE)
          cat(paste('mean length output, ALKS year:',tot$year[i]," Quarter:",tot$quarter[i],'\n'))
        }
        if (tot$first.year.quarter.area[i]) {
          cat(line,file=out,append=TRUE)
          cat(paste("# SMS_area:",tot$SMS_area[i],'\n'),file=out,append=TRUE)
          cat(line,file=out,append=TRUE)
        }
        
        if (tot$first.year.quarter.area.prey[i]) cat(paste("# species:",tot$prey[i],'\n'),file=out,append=TRUE)
        if (tot$first.year.quarter.area.prey.age[i]) cat(paste("# age:",tot$prey.age[i],"first size:",tot$prey.size.class[i],tot$prey.size[i],'\n'),file=out,append=TRUE)
        cat(formatC(tot$ALK.MeanL[i],format="f",dig=0,width=6),file=out,append=TRUE)
        if (tot$last.age[i]) cat('\n',file=out,append=TRUE)
      }
    } else {
      cat(formatC(tot$ALK.MeanL,format="f",dig=0,width=6),file=out,append=TRUE)
    }
    cat(paste("#\n -999.0 # Check sum \n"),file=out,append=TRUE)
  }
  
  
  if (trans.ALK.all) trans.4M.SMS.ALK.all()
  if (trans.ALK.stomach) trans.4M.SMS.ALK()
  
}


# ICES/4M code for species in the same order as given variable name in species_names.in

########################################################################################


## North Sea all other pred (one Mackerel), Two sandeel stocks and hake 
if (TRUE) {
  min.stomach.sampled.in.stratum<-5
  code.name<-       c("OTH","FUL","GLT","HEG","KTW","GBG","GNT","PUF","RAZ","RAJ","GUR","W_H","N_H","GSE","HBP","HAK",'COD','WHG','HAD','POK','MAC','HER','NSA','SSA','NOP','SPR','PLE','SOL')
  code.name.pred<-        c("FUL","GLT","HEG","KTW","GBG","GNT","PUF","RAZ","RAJ","GUR","W_H","N_H","GSE","HBP","HAK",'COD','WHG','HAD','POK','MAC')
  code.name.prey<- c("OTH",'COD','WHG','HAD','HER','NSA','SSA','SPR','NOP')
  selected.years<-c(1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1995,2000,2002,2005,2013)
  
  year.q<-NULL
  var.groups.size<-NULL
  var.groups<-     NULL
  min.stom.groups<-5
  #
  
  RSMS.stom.transform(list.data.path=file.path(root,"SMS-input-key-run-2017"),
                     trans.bio=T, trans.catch=T,
                     trans.meanL=F,  trans.meanL.from.weight=FALSE,
                     trans.stomach=F, trans.ALK.stomach=F, trans.other=T, trans.Consum=F, trans.ALK.all=F,
                     stom.first=1E-06, stom.mid=1E-06, stom.last=1E-06, stom.min.abs=1E-06, delete.tails=TRUE,
                     inserted.haul.no.propor=1.0, 
                     formatted.output=T,selected.years=selected.years,year.q=year.q,min.pred.length=0)
  
}


