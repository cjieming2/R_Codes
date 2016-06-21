setwd("/Users/jiemingchen/Documents/transplantation/")

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

## studies i want
study_ids=c('SDY131','SDY132','SDY133','SDY134','SDY289','SDY290','SDY291',
            'SDY292','SDY294','SDY352','SDY354','SDY355','SDY356','SDY357',
            'SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

donor_ids=c('SDY289','SDY290','SDY291','SDY292','SDY294','SDY356')

## where the serialized data is locally
studies_dir = "/Users/jiemingchen/Documents/transplantation/data"



######################################################

## donors ####
data.donor = read.table("DM_immport2tsv_donors_subjects_7studies.mergeDup.organ.tsv", header = T, sep = "\t")
a.d=data.donor[data.donor$AGEU != "Not_Specified",]
a.d=mutate(a.d, AGE=round(AGE))
b1.d=ddply(.data = a.d, .variables = c("AGE","SEX"), .fun=summarise, count = length(AGE))
b2.d=ddply(.data = a.d, .variables = c("AGE","SEX","Organ"), .fun=summarise, count = length(AGE))
b3.d=ddply(.data = a.d, .variables = c("AGE","SEX","STUDYID"), .fun=summarise, count = length(AGE))

## plot and save
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="SEX"]  <- "Sex"
  
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
    
    pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend + pcolor
    return(myplot)
  } else
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
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  }
}

plot_b1.d = plot_age_sex(b1.d,1)

## to check the range, check the label
## min(round(b2.d$AGE); max(round(b2.d$AGE)
range.d=seq(0.5,length(seq(min(round(b2.d$AGE)),max(round(b2.d$AGE)),by=1))*4+0.5-4,by=4)
plot_b2.d = plot_age_sex(b2.d,2,range.d)

## plot by studies
plot_b3.d = plot_age_sex(b3.d,3,range.d)

x11(type="cairo")
grid.arrange(plot_b1.d, plot_b2.d, plot_b3.d, ncol=1) ## plot x axis is off
# savePlot("DM_age_sex_organ_studies_combined_histogram_all_donors.png", type = "png")

## use ggsave for printing
## ggsave saves only the last plot though, so you can't use grid arrange
# ggsave("Plot2.png")

######
## plot race and ethnicity
a.d$RACE = factor(a.d$RACE, levels=c("American Indian or Alaska Native","Asian","Black or African American",
                                     "Native Hawaiian or Other Pacific Islander","White",
                                     "Other","Not_Specified","Unknown"))

b.d.race = ddply(.data = a.d, .variables = c("RACE","ETHNIC"), .fun=summarise,
                 count = length(RACE))
pmain = ggplot(data=b.d.race, aes(x=RACE, y=count, fill=ETHNIC))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

x11(type="cairo")
plot_race_ethnicity = pmain + phisto + plabels + paxes + ptitle + plegend

b.d.study = ddply(.data = a.d, .variables = c("RACE","STUDYID"), .fun=summarise,
                  count = length(RACE))
pmain.study = ggplot(data=b.d.study, aes(x=RACE, y=count, fill=STUDYID))
plot_race_study = pmain.study + phisto + plabels + paxes + ptitle + plegend

grid.arrange(plot_race_ethnicity, plot_race_study, ncol=1) ## plot x axis is off
# savePlot("DM_race_ethnicity_study_combined_histogram_all_donors.png", type = "png")

