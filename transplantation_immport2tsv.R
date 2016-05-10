setwd("/Users/jiemingchen/Documents/transplantation/")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(RMySQL)

## studies i want
study_ids=c('SDY131','SDY132','SDY133','SDY134','SDY289','SDY290','SDY291',
            'SDY292','SDY294','SDY352','SDY354','SDY355','SDY356','SDY357',
            'SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

recip_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355',
            'SDY356','SDY357','SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

donor_ids=c('SDY289','SDY290','SDY291','SDY292','SDY294','SDY356')

recip_kidney_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355','SDY356',
                   'SDY357','SDY358','SDY546','SDY668','SDY671','SDY674')

## where the serialized data is locally
studies_dir = "/Users/jiemingchen/Documents/transplantation/data"

#### starting up, download and serialize specific studies - do once ####
# ## set up tunnel
# system("/Users/jiemingchen/.start_tunnel.sh")
# 
# ## a script with db connect 
# source("/Users/jiemingchen/R_codes/immport_connect.R")
# 
# ## serialize specific studies on local directory
# serialzeStudyData(study_ids, data_dir=studies_dir)
# 
# ## close tunnel
# system("/Users/jiemingchen/.stop_tunnel.sh")


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
data.recip = loadMultipleStudies(recip_ids,studies_dir,"Demographics")
# write.table(data.recip, "/Users/jiemingchen/Documents/transplantation/immport2tsv_recipients_1272subjects_16studies.tsv", sep="\t",quote = FALSE, row.names = FALSE)



## stats 
data.recip = read.table("immport2tsv_recipients_1272subjects_16studies.tsv", stringsAsFactors = FALSE, header = T, sep = "\t")
library(ggplot2)
a=data.recip[data.recip$AGEU != "Not_Specified",]
a=mutate(a, AGE = round(AGE))



## massage dataframe into a new one
b = ddply(.data = a, .variables = c("AGE","SEX","Organ"), .fun = summarise, 
      count = length(AGE))

## plotting
plot_age_sex = function(data,flag) {
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
    # this one uses ddply  and b
    pmain = ggplot(data=data, aes(x=AGE, y=count, fill=Sex))
    phisto = geom_bar(stat="identity", position = position_dodge())
    plabels = labs(title="Age groups", x="Age",y="Count")
    pticks = scale_x_continuous(breaks=seq(0,80,by=1))
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

    pmain + phisto + plabels + pticks + paxes + ptitle + plegend
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=Organ))
    phisto = geom_bar(stat = "identity", color = "white") 
    plabels = labs(title="Age groups", x="Age",y="Count")
    pticks = scale_x_continuous(breaks=seq(0.5,320.5,by=4), labels=seq(0,80,by=1))
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
  
    pmain + phisto + plabels + pticks + paxes + ptitle + plegend
  }
}

## plot and save
x11(type="cairo")
plot_age_sex(b,1)
# ggsave("age_sex_histogram_recip_all_fill_sex.png", dpi = 150)

## note that because there isnt age 72, 73, 74, 76
## the interaction of age and sex does not plot these values...
x11(type="cairo")
plot_age_sex(b,2)
ggsave("age_sex_histogram_recip_all_fill_organ.png", dpi = 150)

######################################################

## donors ####

## write to file 
data.donor = loadMultipleStudies(donor_ids,studies_dir,"Demographics")
# write.table(data.donor, "/Users/jiemingchen/Documents/transplantation/immport2tsv_donors_10403subjects_6studies.tsv", sep="\t", quote = FALSE, row.names = FALSE)



######################################################


## recipients kidney ####
## load and save to file
data.recip.kidney = loadMultipleStudies(recip_kidney_ids,studies_dir,"Demographics")

#### stats
library(ggplot2)
a=data.recip.kidney[data.recip.kidney$AGEU != "Not_Specified",]
x11(type="cairo")

## plotting
## plot and save
x11(type="cairo")
plot_age_sex(b,1)
# ggsave("age_sex_histogram_recip_all_fill_sex.png", dpi = 150)

## note that because there isnt age 72, 73, 74, 76
## the interaction of age and sex does not plot these values...
x11(type="cairo")
plot_age_sex(b,2)
ggsave("age_sex_histogram_recip_all_fill_organ.png", dpi = 150)

