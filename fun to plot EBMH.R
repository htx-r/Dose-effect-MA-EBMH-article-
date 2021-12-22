#===========================
require(dplyr)
library(wesanderson)
# 1. Table 
tab1 <- function(data=study_87){
  return(data.frame(
    level=0:(nrow(data)-1),
    dose=data$hayasaka_ddd,
    response=data$Responders,
    total=data$No_randomised,
    OR=round(exp(data$logOR),2),
    lb=round(exp(data$logOR - 1.96*data$selogOR),2),
    ub=round(exp(data$logOR + 1.96*data$selogOR),2),
    logOR=round(data$logOR,2),
    selogOR=round(data$selogOR,2)
  ))
}

# fun to Figure 1,2,3
plotdata.fun <- function(drma=rcs_1stage,
                         knots=knots,
                         data=study_87, 
                         p.eff=0.3771
                         ){
  xref = 0
  newd = data.frame(hayasaka_ddd = c(0, seq(0, max(data$hayasaka_ddd), length.out = 50)))
  pred <- predict(drma, newd, xref = xref, expo = T)
  
  # OR 
  OR <- pred$pred
  ubo <- pred$ci.ub
  lbo <- pred$ci.lb
  
  # prob 
  odds0 <- (p.eff/(1-p.eff))
  odds <- OR*odds0
  prob <- odds/(1+odds)
  ubp <- (ubo*odds0)/(1+(ubo*odds0))
  lbp <- (lbo*odds0)/(1+(lbo*odds0))
  
  # >> to return
  plotdata <- data.frame(dose=pred[,1],
                        # odds ratio 
                         OR=OR,
                         ubo=ubo,
                         lbo=lbo,
                         
                        # probabilities
                         prob=prob,
                         ubp=ubp,
                         lbp=lbp,
                         p.eff=p.eff)
  return(plotdata)
  }
  

#===========================
# plot
doseres.plot <- function(plotdata,
                          data,
                          ymax=5,
                          ymin=0,
                          y='OR',
                          ub='ubo',
                          lb='lbo',
                          col=c("darkred","mistyrose3"),#c('grey27','gray70'),
                          labs=c('Odds Ratio (OR)','Fluoxetine-equivalent dose'),
                          linetype='solid',
                          add2=NULL,
                          add3=NULL){
  
  # rename
  names(plotdata)[names(plotdata)==y] <- 'y'
  names(plotdata)[names(plotdata)==ub] <- 'ub'
  names(plotdata)[names(plotdata)==lb] <- 'lb'

  
    # plot OR vs dose with CI
    g1 <- ggplot(data = plotdata,aes(x=dose,y=y))+
      geom_smooth(                                    
        aes(x=dose,y=y,ymin=lb,ymax=ub),
        color=col[1],fill=col[2],
        data=plotdata, stat="identity",linetype=linetype,size=1.5)+
      coord_cartesian(ylim = c(ymin, ymax))
  # add rug of observed doses
  g2 <- g1+geom_rug(data=data,mapping = aes(x=hayasaka_ddd),
                   inherit.aes = F
                   ,col='royalblue4'
                   )
  
  # change labs
  g3 <- g2 + ggplot2::labs(y=labs[1], x=labs[2])
  
  
  # edit the text size
  g4 <- g3+#+theme(
    #panel.background = element_rect(fill = 'honeydew',colour = 'white'),
    # axis.text.x = element_text(face='bold',size=16),axis.text.y = element_text(face='bold',size=14),
    # axis.title.x=element_text(size=16,face = "bold"),axis.title.y=element_text(size=16,face = "bold"),
    # #strip.background =element_rect(fill="snow3"),
    # strip.text.x = element_text(size = 16),
    # legend.text = element_text(size=12))+
    theme_set(theme_minimal(
      base_size = 30
    ))
  
  if(!is.null(add2)){
    g.add2 <- g4+geom_smooth(                                    
      aes(x=dose,y=OR,ymin=lbo,ymax=ubo),
      #color='darkred',fill='lightcoral',
      color='grey17',fill='gray65',
      data=add2, stat="identity",linetype='dashed',size=1.5)
    g.add <- g.add2
    
    if(!is.null(add3)){
      g.add3 <- g.add2+geom_smooth(                                    
        aes(x=dose,y=OR,ymin=lbo,ymax=ubo),
        color='pink4',fill='pink',
        data=add3, stat="identity",linetype='dotted')
      g.add <- g.add3
    }else{
      g.add <- g.add2
    }
  }else{
    g.add <- g4 
  }
  g <- g.add
  if(y=='prob'){
    # add placebo effect
    g <- g+geom_hline(yintercept=plotdata$p.eff, color='grey25',linetype='dashed',size=1) # add placebo effect
  } else{
    g <- g
  }
  g
}
#



# appendix figure 1: dose distribution for flux_dose
dose_dist1 <- function(data=antidep){
  ggplot(data = data, aes(x='SSRI',y=hayasaka_ddd)) +
    xlab('')+
    ylab('')+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1,col='orange')+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_text(size = 16,family="serif"),
          axis.text.y=element_text(size=12))+theme_set(
            theme_minimal() 
          )
}


# appendix figure 2:  dose distribution per drug

dose_dist2 <- function(data=antidep){
  data2 <- data[!data$Drug%in%c('placebo'),]
  ggplot(data = data2, aes(x=Drug, y=Dose_delivered_mean)) +
    facet_wrap(~Drug, scales="free")+
    xlab('')+
    ylab('')+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=3,col='orange')+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_text(size = 16,family="serif"),
          axis.text.y=element_text(size=12))
}


# appendix table 1: characterstics
app.tab <- function(data=antidep){
  # names
  drug.name <- levels(factor(eval(substitute(Drug), data)))
  
  # create a list of drugs - row with 
  # n.events, n.studies, n.doses and n.patients
  data.drug <- sapply(drug.name, 
                      function(d){
                        data_per_drug <- data%>%filter(data$Drug==d)
                        n.events <- data_per_drug%>%select(Responders)%>%sum()
                        n.studies <-data_per_drug%>%select(Study_No)%>%unique()%>%nrow()
                        n.doses <-data_per_drug%>%nrow()
                        n.patients <- data_per_drug%>%select(No_randomised)%>%sum()
                        
                        return(c(n.events,
                                 n.patients,
                                 n.studies,
                                 n.doses
                        ))
                      }
                      ,simplify = F)
  
  # merge rows
  tbl <- do.call(rbind,data.drug)
  # to table format
  colnames(tbl) <- c(
    "Number of events",
    "Number of patients",
    "Number of studies",
    "Number of non-zero doses")
  add_tab <- as_tibble(tbl,rownames = NA)
  return(add_tab)
  
}





