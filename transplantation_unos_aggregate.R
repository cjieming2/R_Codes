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

## living donors
data.donor = read.table("unos_2011_jul2016_livingdonors.txt", header = T, sep = "\t", stringsAsFactors = F)

## plot
x11(type="cairo")
barplot(data.donor$Male, xlab="donation age", ylab="count", names.arg=data.donor$age)
barplot(rbind(data.donor$Male,data.donor$Female), xlab="Donation age", ylab="number", names.arg=data.donor$age, beside=TRUE, col=c("blue","red"))
legend(3,7000,c("Male","Female"),col=c("blue","red"),pch=15, cex = 2.5)
# savePlot("unos_2011_jul2016_livingdonors_kidney.png", type=c("png"))


## age, sex and ethnicity
myaggdata = read.table("unos_2011_jul2016_livingdonors_ethnicity_age_sex.txt", header = T, sep = "\t", stringsAsFactors = F)

## tell R that the columns are ordered already
myaggdata$age = factor(myaggdata$age, levels = myaggdata$age)

## plot 
plot_age_sex_here = function(data)
{
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = interaction(gender,age), y=count, fill=ethnicity))
    phisto = geom_bar(stat = "identity") 
    plabels = labs(x="Age",y="Count")
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 18, angle = 90, hjust = 1), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
}

# plot4me = plot_age_sex(myaggdata)
x11(type="cairo")
plot_age_sex_here(myaggdata)
# savePlot("unos_2011_jul2016_ldonors_ethnicity_age_sex_kidney.png", type=c("png"))


###########

## deceased donors
d.donor = read.table("unos_2011_jul2016_ddonors.txt", header = T, sep = "\t", stringsAsFactors = F)

## plot
x11(type="cairo")
barplot(d.donor$Male, xlab="donation age", ylab="count", names.arg=d.donor$age)
barplot(rbind(d.donor$Male,d.donor$Female), xlab="Donation age", ylab="number", names.arg=d.donor$age, beside=TRUE, col=c("blue","red"))
legend(3,7000,c("Male","Female"),col=c("blue","red"),pch=15, cex = 2.5)
# savePlot("unos_2011_jul2016_ddonors_kidney.png", type=c("png"))

###########
setwd("/Users/jiemingchen/Documents/transplantation/a_recipient")

## recipients
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="gender"]  <- "Sex"
  names(data)[names(data)=="age"]  <- "AGE"
  
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
    #pticks = scale_x_continuous(AGE)
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  } else if(flag == 2)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=Organ))
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
    
    pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend + pcolor
    return(myplot)
  } else if(flag == 3)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots study with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=STUDYID))
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
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
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
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  }
}

d.recip = read.table("unos_2008_jul2016_recipient_kidney.txt", header = T, sep = "\t", stringsAsFactors = F)

## tell R that the columns are ordered already
d.recip$age = factor(d.recip$age, levels = d.recip$age)

x11(type="cairo")
plot_age_sex(d.recip,1,"")
# savePlot("unos_2008_jul2016_recip_kidney.png", type=c("png"))


#----
## age, sex and ethnicity
myaggdata2 = read.table("unos_2008_jul2016_age_gender_ethnicity_recipient_kidney.txt", header = T, sep = "\t", stringsAsFactors = F)

## tell R that the columns are ordered already
myaggdata2$age = factor(myaggdata2$age, levels = myaggdata2$age)

## plot 
plot_age_sex = function(data)
{
  ########
  ## this plots each x axis as an interaction of age and sex and then adjust label for age
  ## order of Sex and AGE matter - uses b
  ## plots organ with this data
  pmain = ggplot(data, aes(x = interaction(gender,age), y=count, fill=ethnicity))
  phisto = geom_bar(stat = "identity") 
  plabels = labs(x="Age",y="Count")
  paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                axis.text.x = element_text(size = 18, angle = 90, hjust = 1), axis.text.y = element_text(size = 15))
  ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
  plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
  
  myplot = pmain + phisto + plabels + paxes + ptitle + plegend
  return(myplot)
}

# plot4me = plot_age_sex(myaggdata)
x11(type="cairo")
plot_age_sex(myaggdata2)
# savePlot("unos_2011_jul2016_recip_ethnicity_age_sex_kidney.png", type=c("png"))


#####
## waitlist
setwd("/Users/jiemingchen/Documents/transplantation/a_recipient")

## recipients
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="gender"]  <- "Sex"
  names(data)[names(data)=="age"]  <- "AGE"
  
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
    #pticks = scale_x_continuous(AGE)
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  } else if(flag == 2)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=Organ))
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
    
    pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend + pcolor
    return(myplot)
  } else if(flag == 3)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots study with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=STUDYID))
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
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
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
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for recips
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  }
}

d.recip = read.table("unos_2008_jul2016_recipient_waitlist_kidney.txt", header = T, sep = "\t", stringsAsFactors = F)

## tell R that the columns are ordered already
d.recip$age = factor(d.recip$age, levels = d.recip$age)

x11(type="cairo")
plot_age_sex(d.recip,1,"")
# savePlot("unos_2008_jul2016_recip_waitlist_kidney.png", type=c("png"))


#----
## age, sex and ethnicity
myaggdata2 = read.table("unos_2008_jul2016_age_gender_ethnicity_recipient_waitlist_kidney.txt", header = T, sep = "\t", stringsAsFactors = F)

## tell R that the columns are ordered already
myaggdata2$age = factor(myaggdata2$age, levels = myaggdata2$age)

## plot 
plot_age_sex2 = function(data)
{
  ########
  ## this plots each x axis as an interaction of age and sex and then adjust label for age
  ## order of Sex and AGE matter - uses b
  ## plots organ with this data
  pmain = ggplot(data, aes(x = interaction(gender,age), y=count, fill=ethnicity))
  phisto = geom_bar(stat = "identity") 
  plabels = labs(x="Age",y="Count")
  paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                axis.text.x = element_text(size = 18, angle = 90, hjust = 1), axis.text.y = element_text(size = 15))
  ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
  plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
  
  myplot = pmain + phisto + plabels + paxes + ptitle + plegend
  return(myplot)
}

# plot4me = plot_age_sex(myaggdata)
x11(type="cairo")
plot_age_sex2(myaggdata2)
# savePlot("unos_2011_jul2016_recip_waitlist_ethnicity_age_sex_kidney.png", type=c("png"))
