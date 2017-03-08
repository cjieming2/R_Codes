setwd("/Users/jiemingchen/Documents/transplantation/a_donor/unos_optn")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(RMySQL)
library(ggplot2)
library(gridExtra)
library(gapminder)
library(reshape2)
library(reshape)
library(lattice)


## where the serialized data is locally
studies_dir = "/Users/jiemingchen/Documents/transplantation/a_donor/"



######################################################

## donors ####
data.donor = read.table("unos-LIVING_DONOR_DATA_DM_race_organ_relation-kidney-us-128407.txt", header = T, sep = "\t")
a.d = data.donor
b1.d=ddply(.data = a.d, .variables = c("AGE_DON","GENDER"), .fun=summarise, count = length(AGE_DON))
b3.d=ddply(.data = a.d, .variables = c("AGE_DON","GENDER","ETHCAT_DON"), .fun=summarise, count = length(AGE_DON))

## plot and save
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="GENDER"]  <- "Sex"
  names(data)[names(data)=="AGE_DON"]  <- "AGE"
  names(data)[names(data)=="ETHCAT_DON"]  <- "RACE"
  
  ## this is a histogram, age 40 is considered at bin 39; uses un-ddply a - nope
  # pmain = ggplot(data=data, aes(x=round(data$AGE), fill=Sex))
  # phisto = geom_histogram(breaks=seq(min(round(data$AGE)),round(max(data$AGE)),by=1), position = "dodge")
  
  ## this groups each age group as a category and then takes 15 min to plot!! uses ddply b -- nope
  # pmain = ggplot(b,aes(x=SEX,y=count,fill=Organ))
  # phisto = geom_bar(stat = "identity",color="white") + facet_wrap(~AGE,nrow=1)
  
  if(flag == 1)
  {
    ############
    # this one b1 
    pmain = ggplot(data=data, aes(x=AGE, y=count, fill=Sex))
    phisto = geom_bar(stat="identity", position = position_dodge())
    plabels = labs(x="Age",y="Count")
    pticks = scale_x_continuous(breaks=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    myplot = pmain + phisto + plabels + pticks + paxes + ptitle + plegend
    return(myplot)
  } else if(flag == 2)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=ORGAN))
    phisto = geom_bar(stat = "identity") 
    plabels = labs(x="Age",y="Count")
    # pticks = scale_x_continuous(breaks=range, labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    #pticks = scale_x_continuous(breaks=range,labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    # pticks = scale_x_continuous() ##debug
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend + pcolor
    return(myplot)
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots study with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=RACE))
    phisto = geom_bar(stat = "identity") 
    plabels = labs(x="Age",y="Count")
    # pticks = scale_x_continuous(breaks=range, labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    #pticks = scale_x_continuous(breaks=range,labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    # pticks = scale_x_continuous() ##debug
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  }
}

plot_b1.d = plot_age_sex(b1.d,1)

# ## to check the range, check the label
# ## min(round(b2.d$AGE); max(round(b2.d$AGE)
# range.d=seq(0.5,length(seq(min(round(b2.d$AGE)),max(round(b2.d$AGE)),by=1))*4+0.5-4,by=4)
# plot_b2.d = plot_age_sex(b2.d,2,range.d)
# 
# ## plot by race
plot_b3.d = plot_age_sex(b3.d,3,range.d)

x11(type="cairo")
grid.arrange(plot_b1.d, plot_b3.d, ncol=1) ## plot x axis is off
# savePlot("unosjun2016_age_sex_histogram_kidney_donors.png", type = "png")

## use ggsave for printing
## ggsave saves only the last plot though, so you can't use grid arrange
# ggsave("Plot2.png")

######
## group by race
pmain = ggplot(data=a.d, aes(x=factor(ETHCAT_DON)))
phisto = geom_bar(stat="count")
plabels = labs(x="Racial categories", y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

plot_race_ethnicity = pmain + phisto + plabels + paxes + ptitle + plegend


# plot by age and gender on xaxis
b.d.age = ddply(.data = a.d, .variables = c("AGE_DON","ETHCAT_DON"), .fun=summarise,
                count = length(AGE_DON))
pmain.age = ggplot(data=b.d.age, aes(x=AGE_DON, y=count, fill=ETHCAT_DON))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)
plot_race_age = pmain.age + phisto + plabels + paxes + ptitle + plegend

## combine plot
x11(type="cairo")
grid.arrange(plot_race_ethnicity, plot_race_age, ncol=1) ## plot x axis is off
# savePlot("unos-jun2016_age_race_sex_ethnicity_histogram_kidney_donors.png", type = "png")


#####
## split by race
b2.d=ddply(.data = a.d, .variables = c("AGE_DON","ETHCAT_DON"), .fun=summarise, count = length(AGE_DON))
pmain = ggplot(data=b2.d, aes(x = AGE_DON, y = count, fill = ETHCAT_DON))
phisto = geom_bar(stat="identity", position = "stack")
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b1.d$AGE_DON),max(b1.d$AGE_DON),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

jm = pmain + phisto + plabels + paxes + pticks + ptitle + plegend + scale_fill_brewer(palette = "Set1")

## facet
# jm + facet_grid(ETHCAT_DON ~ .) + theme(legend.position = "none")
# savePlot("unos-jun2016_age_splitbyrace_histogram_kidney_donors.png", type = "png")

## scaled facet to see those of very low coutns (check y axis) 
jmjm =jm + facet_grid(ETHCAT_DON ~ ., scale = "free_y") + theme(legend.position = "none")

## include the total
total = cast(b2.d, AGE_DON ~ ., sum, value="count")
total = rename(total, c(`(all)` = "count"))
total$ETHCAT_DON <- "zTotal"
b2.d.total <- rbind(total, b2.d) ## it auto rearranges the columns
jmjmjm = jmjm %+% b2.d.total
jmjmjm
# savePlot("unos-jun2016_age_splitbyrace_histogram_kidney_donors_scaled.png", type = "png")
  
#####
## split by race and gender using lattice
total3 = cast(b3.d, AGE_DON + GENDER ~ ., sum, value="count")
total3 = rename(total3, c(`(all)` = "count"))
total3$ETHCAT_DON <- "zTotal"
b3.d.total <- rbind(total3, b3.d) ## it auto rearranges the columns

barchart(count ~ AGE_DON | factor(ETHCAT_DON), b3.d.total,  groups = GENDER, type = "l", strip = FALSE, strip.left = TRUE,
       layout = c(1, 9), horizontal=FALSE, scales = list(y = "free"),
       index.cond=list(seq(9,1)), ## this gives order of panels
       ylab = "count",  auto.key = list(lines = TRUE, points = FALSE, columns = 2))

# savePlot("unos-jun2016_age_gender_splitbyrace_histogram_kidney_donors_lattice.png", type = "png")

######
## split by race and gender using ggplot2
pmain = ggplot(b3.d.total, aes(AGE_DON, count, fill = GENDER))
phisto = geom_bar(stat = "identity", position=position_dodge())
pfacet = facet_grid(ETHCAT_DON ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b1.d$AGE_DON),max(b1.d$AGE_DON),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend
# savePlot("unos-jun2016_age_gender_splitbyrace_histogram_kidney_donors_ggplot2.png", type = "png")

## stacked
pmain = ggplot(b3.d.total, aes(AGE_DON, count, fill = GENDER))
phisto = geom_bar(stat = "identity")
pfacet = facet_grid(ETHCAT_DON ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b1.d$AGE_DON),max(b1.d$AGE_DON),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend
# savePlot("unos-jun2016_age_gender_splitbyrace_histogram_kidney_donors_ggplot2_stacked.png", type = "png")

## overlapping bar plots
pmain = ggplot(b3.d.total, aes(AGE_DON, count, fill = GENDER))
phisto = geom_bar(stat="identity", position = "identity",  alpha=.3)
pfacet = facet_grid(ETHCAT_DON ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b1.d$AGE_DON),max(b1.d$AGE_DON),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend
savePlot("unos-jun2016_age_gender_splitbyrace_histogram_kidney_donors_ggplot2_transpplot.png", type = "png")


## line plot
pmain = ggplot(b3.d.total, aes(AGE_DON, count, color = GENDER))
phisto = geom_line()
pfacet = facet_grid(ETHCAT_DON ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b1.d$AGE_DON),max(b1.d$AGE_DON),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend
# savePlot("unos-jun2016_age_gender_splitbyrace_histogram_kidney_donors_ggplot2_lineplot.png", type = "png")

