setwd("/Users/jiemingchen/Documents/transplantation/a_recipient")

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


######################################################

## donors ####
data.recipient = read.table("emr_ucsf_2016jun_kidney_recipients_with_inpatient_bill_2120subj.txt", header=T, sep="\t", stringsAsFactors = F)
a.d = data.recipient

b1.d=ddply(.data = a.d, .variables = c("Billing_Age2","Patient_Sex"), .fun=summarise, count = length(Billing_Age))
b3.d=ddply(.data = a.d, .variables = c("Billing_Age2","Patient_Sex","Patient_Race2"), .fun=summarise, count = length(Billing_Age))

## plot and save
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="Patient_Sex"]  <- "Sex"
  names(data)[names(data)=="Billing_Age2"]  <- "AGE"
  names(data)[names(data)=="Patient_Race2"]  <- "RACE"
  
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
# savePlot("emrjun2016_ucsf_age_sex_histogram_kidney_recipients.png", type = "png")

# use ggsave for printing
## ggsave saves only the last plot though, so you can't use grid arrange
# ggsave("Plot2.png")

######
## plot race and ethnicity
a.d$Patient_Race2 = factor(a.d$Patient_Race2, levels=c("American_Indian_or_Alaska_Native","Asian",
                                                       "Black_or_African_American",
                                                       "Native_Hawaiian_or_Other_Pacific_Islander",
                                                       "White_or_Caucasian","Other","Not_Specified",
                                                       "Unknown/Declined" ,"multiracial"))

b.d.race = ddply(.data = a.d, .variables = c("Patient_Race2","Patient_Ethnicity"), .fun=summarise,
                 count = length(Patient_Race2))
pmain = ggplot(data=b.d.race, aes(x=Patient_Race2, y=count, fill=Patient_Ethnicity))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

plot_race_ethnicity = pmain + phisto + plabels + paxes + ptitle + plegend

b.d.age = ddply(.data = a.d, .variables = c("Billing_Age2","Patient_Race2"), .fun=summarise,
                  count = length(Billing_Age2))
pmain.age = ggplot(data=b.d.age, aes(x=Billing_Age2, y=count, fill=Patient_Race2))
plot_race_age = pmain.age + phisto + plabels + paxes + ptitle + plegend

x11(type="cairo")
grid.arrange(plot_race_ethnicity, plot_race_age, ncol=1) ## plot x axis is off
# savePlot("emrjun2016_ucsf_age_race_sex_ethnicity_histogram_kidney_recipients.png", type = "png")


