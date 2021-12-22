source('analyze EBMH.R')
source('fun to plot EBMH.R')

# Table 1
tab1()

# Figure 1
plotdata1s = plotdata.fun(drma = rcs_1study,
                         data = study_87,
                         knots=knots)
plotdata2s = plotdata.fun(drma = lin_1study,
                          data = study_87,
                          knots=knots)

doseres.plot(plotdata = plotdata1s,
              data=study_87,
              ymax = 4.6,
              ymin=1,
              y='OR',
              ub='ubo',
              lb='lbo',
              add2=plotdata2s)


# Figure 2: OR vs dose - RCS: 2stage & 1stage + Linear$
plotdata1 = plotdata.fun(drma = rcs_pooled2,
                         data = antidep,
                         knots=knots,
                         p.eff = pl_eff)

plotdata2 = plotdata.fun(drma = rcs_pooled1,
                         data = antidep,
                         knots=knots,
                         p.eff = pl_eff)




doseres.plot(plotdata =plotdata1,
             data=antidep,
             ymax = 2,
             ymin=0.5,
             y='OR',
             ub='ubo',
             lb='lbo',
             add2=plotdata2,
             add3=NULL)

# Figure 3: prob vs dose - 1stage
doseres.plot(plotdata =plotdata2,
             data=antidep,
             ymax = 0.55,
             ymin=0.3,
             y='prob',
             ub='ubp',
             lb='lbp',
             labs = c('Predicted absolute effect','Fluoxetine-equivalent dose'))



# Figure 4: VPC vs dose

# VPC 
df <- antidep[!is.na(antidep$selogOR),]
df$vpc <- vpc(rcs_pooled1)
min(df$vpc[df$hayasaka_ddd==20])
max(df$vpc[df$hayasaka_ddd==20])
max(df$vpc)
ggplot(df, aes(hayasaka_ddd,vpc)) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.3, fill= "darkseagreen2", alpha=0.3)  +
  theme_light()+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.3, ymax = 0.6 , fill= "lightgoldenrod2", alpha=0.3) + 
  theme_light()+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = Inf, fill= "lightcoral", alpha=0.3) +
  theme_light()+
  geom_point() + 
  labs(y="variance partition component (VPC)", x="Fluoxetine-equivalent dose")+
  geom_smooth(method = "loess",se=FALSE)+
  # ylim(0,0.9)+
  coord_cartesian(clip="off", ylim=c(0,0.9))+
  theme(axis.title.x=element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"),
        plot.margin = unit(c(5,10,10,5), "mm"),
        axis.text.x = element_text(face='bold',size=14),axis.text.y = element_text(face='bold',size=14),)


# appendix figure 1
dose_dist1()

# appendix figure 2
dose_dist2()

# 
ggplot(df, aes(hayasaka_ddd,vpc)) + 
  coord_flip()+
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2),color="slateblue4")+
  coord_cartesian(clip="off", ylim=c(0,0.9))+
  theme(axis.title=element_blank(),
        plot.margin = unit(c(5,10,10,5), "mm"))
# appendix table
app.tab()

