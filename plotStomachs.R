

#' Plot diet data.
#'
#' @param d Diet data set of class STOMdiet.
#' @param cut_pred_size From to in substring of predator size
#' @param show_plot Show the resulting graphs on screen (or save the results for later processing)
#' @param addTitle Add predator name on top of the plot.
#' @param tAngle Angle X-axis text.
#' @param addNstom Show number of stomachs sampled
#' @param nstomAngle Angle number of stomachs text
#' @param Ncol number of columns in facet plot
#' @param Nrow number of rows in facet plot
#' @param Strip.position strip.position: "top" | "bottom" | "left" | "right"
#' @param Colours vector of colours for preys.
#' @param otherCol Colour for "other prey"
#' @param refac_prey Reorder preys
#' @return nothing (if show_plot=TRUE) or a list of plots.
#' @importFrom ggplot2 ggplot facet_wrap geom_col labs geom_text theme theme_minimal scale_fill_manual aes element_text element_line
#' @importFrom rlang .data
#' @method  plot STOMdiet
#' @export
plot.STOMdiet<-function(d,show_plot=TRUE,cut_pred_size=c(1,10),addTitle=FALSE,tAngle=90,addNstom=FALSE,nstomAngle=45,Ncol=2,Nrow=NULL,
                    Strip.position = c("top", "bottom", "left", "right"),Colours,otherCol='grey',refac_prey=FALSE) {

  if (missing(Colours)) Colours<-c('red','green','plum','blue','cyan','yellow','coral','skyblue','purple','magenta','limegreen','pink' )
  Strip.position <- match.arg(Strip.position)

  if (refac_prey) d<-refac_prey(d)
  Colours[nlevels(d[['PREY']]$prey_name)]<-otherCol
  pn<-levels(d[['PREY']]$prey_name)
  allNames<- Colours[1:length(pn)]
  names(allNames)<-pn

 
  x<-x %>%  dplyr::group_by(key) %>% dplyr::mutate(prey_w=prey_w/sum(prey_w)*100) %>% dplyr::ungroup() %>%
    dplyr::mutate(pred_size=substr(pred_size,cut_pred_size[1],cut_pred_size[2])) %>%
    dplyr::group_by(stratum_time,pred_name, pred_size,n_tot,year,quarter,prey_name) %>% dplyr::summarise(prey_w=sum(prey_w)) %>%
    dplyr::ungroup() %>%
    dplyr::select(stratum_time,pred_name, pred_size,n_tot,year,quarter,prey_name,prey_w)


  out<-by(x,list(x$pred_name,x$year),function(x) {
    if (addTitle) tit<- as.character(x$pred_name)[1] else tit<-NULL
    a<-ggplot(x) +
       facet_wrap(~stratum_time,  ncol=Ncol, nrow=Nrow, strip.position = Strip.position)+
       scale_fill_manual(  values = allNames,name='Prey')+
       geom_col(aes(x=pred_size, y = prey_w, fill = prey_name ))+
       labs(x='Predator size',y='weight percentage',title=tit)+
      theme_minimal() +
      theme( panel.grid.major = element_line(linetype = "blank"),
             panel.grid.minor = element_line(linetype = "blank"),
             axis.text.x = element_text(angle = tAngle, vjust = 0.5),
             )

    if (addNstom) a<-a+geom_text(aes(pred_size, one, label = n_tot),
              vjust = 0.5, angle = nstomAngle, size=3,
              data = . %>% dplyr::select(pred_size,stratum_time,n_tot) %>% unique() %>%
                dplyr::mutate(one=15) %>% dplyr::arrange(stratum_time,pred_size))

    if (show_plot) print(a) else return(a)
  })
  if (show_plot) return() else return(out)
}

data$stom
#  stom2<-unnest(rep1$stom,cols = c(data))
stom2<-unnest(data$stom,cols = c(data))

stom2<-unnest(stom2,cols = c(data))
tst<-select(stom2,area,y,q,pred,predSizeClass,prey,preySizeClass)
dim(tst); dim(unique(tst))
tst<-stom2 %>% mutate(dup=duplicated(tst), n=1:dplyr::n())
filter(tst,dup)$n
view(tst)

d<-stom2 %>% mutate(Year=y-data$off.year, yq=paste0(Year,'-Q',q), Predator=data$predNames[pred],Prey=c(data$preyNames,data$otherFoodName)[prey]) %>%
  mutate(Predator=factor(Predator,data$predNames),Prey=factor(Prey,c(data$preyNames,data$otherFoodName))) %>%
  group_by(area,Year,q,yq,pred,Predator, predSizeClass,prey,Prey,noStom,noHaul) %>% summarize(stomcon=sum(stomcon),   Estom=sum(Estom)) %>% ungroup()



filter(stom2,is.na(Prey))
