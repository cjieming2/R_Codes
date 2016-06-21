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

recip_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355',
            'SDY356','SDY357','SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

recip_kidney_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355','SDY356',
                   'SDY357','SDY358','SDY546','SDY668','SDY671','SDY674')

## where the serialized data is locally
studies_dir = "/Users/jiemingchen/Documents/transplantation/data"


#### studies stats ####

# list.files(studies_dir)


#### subject info (Demographics) ####

## loop through the studies and compile
## checked that the total is just 1000 - small - safe to compile
#study_id = "SDY356"

## recipients ####
loadMultipleStudies = function(ids,studies_dir,domain="Demographics") {
  for (i in 1:length(ids))
  {
    ## load serialized study data into R
    data.study = as.data.frame(loadSerializedStudyData(studies_dir, ids[i], domain)[1])
  
    ## gender and samplesize stats
    if(!exists("data.recip"))
    {
      data.recip = data.study
    } else
    {
      data.recip = rbind(data.recip,data.study)
    }
  
  }
  
  return(data.recip)
}

## load serialized data and write to file
# data.recip = loadMultipleStudies(recip_ids,studies_dir,"Demographics")
# write.table(data.recip, "/Users/jiemingchen/Documents/transplantation/immport2tsv_recipients_1272subjects_16studies.tsv", sep="\t",quote = FALSE, row.names = FALSE)



## stats 
data.recip = read.table("immport2tsv_recipients_1204subjects_16studies_no_dups.tsv", header = T, sep = "\t")

## there are 6 Not_Specified and 5 NA in AGE with 1 '0'
## so AGEU is a better filter
a=data.recip[data.recip$AGEU != "Not_Specified",]
a=mutate(a, AGE = round(AGE))

## plotting
## you have remassage the df in each case
## option 1: takes parameters with counts from AGE and split by SEX
## option 2: takes parameters with counts from AGE, and spilt by SEX and Organ
## both returns handles for plots
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
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=Organ))
    phisto = geom_bar(stat = "identity", color = "white") 
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
  }
}

## plot and save
b1 = ddply(.data = a, .variables = c("AGE","SEX"), .fun=summarise,
                                     count = length(AGE))
# x11(type="cairo")
plot_b1 = plot_age_sex(b1,1)
# ggsave("age_sex_histogram_recip_all_fill_sex.png", dpi = 150)

## note that because there isnt age 72, 73, 74, 76
## the interaction of age and sex does not plot these values...
## includes organs in this plot
## massage dataframe into a new one
b2 = ddply(.data = a, .variables = c("AGE","SEX","Organ"), .fun = summarise, 
           count = length(AGE))
# x11(type="cairo")
range=seq(min(round(b2$AGE))+0.5,round(max(b2$AGE))*4+0.5,by=4)
plot_b2 = plot_age_sex(b2,2,range)
# ggsave("age_sex_histogram_recip_all_fill_organ.png", dpi = 150)


x11(type="cairo")
# multiplot(plot_b1, plot_b2, cols=1)
grid.arrange(plot_b1, plot_b2, ncol=1)
## plot x axis is somewhat off due to missing values in some of the age groups
# savePlot("age_sex_organ_combined_histogram_all_recipients.png", type = "png")


######
## plot race and ethnicity
a$RACE = factor(a$RACE, levels=c("American Indian or Alaska Native","Asian","Black or African American",
                                 "Native Hawaiian or Other Pacific Islander","White",
                                 "Other","Not_Specified","Unknown"))

# b.race = ddply(.data = a, .variables = c("RACE","ETHNIC"), .fun=summarise,
               # count = length(RACE), pos = cumsum(count) - 0.5*count)
b.race = ddply(.data = a, .variables = c("RACE","ETHNIC"), .fun=summarise,
               count = length(RACE))
data=b.race
pmain = ggplot(data=data, aes(x=RACE, y=count, fill=ETHNIC))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

x11(type="cairo")
pmain + phisto + plabels + paxes + ptitle + plegend
# savePlot("race_ethnicity_combined_histogram_all_recipients.png", type = "png")

######################################################


## recipients kidney ####
## load and save to file
# data.recip.kidney = loadMultipleStudies(recip_kidney_ids,studies_dir,"Demographics")
data.recip.kidney = read.table("immport2tsv_recipients_1145subjects_16studies_no_dups_kidney_only.tsv", header = T, sep = "\t")

#### stats
a.k=data.recip.kidney[data.recip.kidney$AGEU != "Not_Specified",]
a.k=mutate(a.k, AGE = round(AGE))
b1.k = ddply(.data = a.k, .variables = c("AGE","SEX"), .fun=summarise, count = length(AGE))

## plot and save
x11(type="cairo")
plot_age_sex(b1.k,1)
# savePlot("age_sex_organ_combined_histogram_kidney_recipients.png", type = "png")

######
## plot race and ethnicity
a.k$RACE = factor(a.k$RACE, levels=c("American Indian or Alaska Native","Asian","Black or African American",
                                 "Native Hawaiian or Other Pacific Islander","White",
                                 "Other","Not_Specified","Unknown"))

b.k.race = ddply(.data = a.k, .variables = c("RACE","ETHNIC"), .fun=summarise,
               count = length(RACE))
data.k=b.k.race
pmain = ggplot(data=data.k, aes(x=RACE, y=count, fill=ETHNIC))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

x11(type="cairo")
pmain + phisto + plabels + paxes + ptitle + plegend
# savePlot("race_ethnicity_combined_histogram_kidney_recipients.png", type = "png")

######################################################

## donors ####

## write to file 
# data.donor = loadMultipleStudies(donor_ids,studies_dir,"Demographics")
# write.table(data.donor, "/Users/jiemingchen/Documents/transplantation/immport2tsv_donors_10403subjects_6studies.tsv", sep="\t", quote = FALSE, row.names = FALSE)

data.donor = read.table("immport2tsv_donors_10364subjects_7studies_from_recip_cohorts.tsv", header = T, sep = "\t")
a.d=data.donor[data.donor$AGEU != "Not_Specified",]
a.d=mutate(a.d, AGE=round(AGE))
b1.d=ddply(.data = a.d, .variables = c("AGE","SEX"), .fun=summarise, count = length(AGE))
b2.d=ddply(.data = a.d, .variables = c("AGE","SEX","Organ"), .fun=summarise, count = length(AGE))

## plot and save
plot_b1.d = plot_age_sex(b1.d,1)

## to check the range, check the label
## min(round(b2.d$AGE); max(round(b2.d$AGE)
range.d=seq(0.5,length(seq(min(round(b2.d$AGE)),max(round(b2.d$AGE)),by=1))*4+0.5-4,by=4)
plot_b2.d = plot_age_sex(b2.d,2,range.d)

x11(type="cairo")
grid.arrange(plot_b1.d, plot_b2.d, ncol=1)## plot x axis is off
# savePlot("age_sex_organ_combined_histogram_all_donors.png", type = "png")

######
## plot race and ethnicity
a.d$RACE = factor(a.d$RACE, levels=c("American Indian or Alaska Native","Asian","Black or African American",
                                     "Native Hawaiian or Other Pacific Islander","White",
                                     "Other","Not_Specified","Unknown"))

b.d.race = ddply(.data = a.d, .variables = c("RACE","ETHNIC"), .fun=summarise,
                 count = length(RACE))
data.d=b.d.race
pmain = ggplot(data=data.d, aes(x=RACE, y=count, fill=ETHNIC))
phisto = geom_bar(stat="identity")
plabels = labs(x="Racial categories",y="Count")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(angle = 45, size = 15, hjust = 1), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
# ptext = geom_text(aes(label=count, y = pos, vjust=-0.5), size = 6)

x11(type="cairo")
pmain + phisto + plabels + paxes + ptitle + plegend
# savePlot("race_ethnicity_combined_histogram_all_donors.png", type = "png")

