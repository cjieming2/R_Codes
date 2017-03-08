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
## left out non-clinical-trials e.g. SDY788 
# study_ids=c('SDY131','SDY132','SDY133','SDY134','SDY289','SDY290','SDY291',
#             'SDY292','SDY293','SDY294','SDY352','SDY354','SDY355','SDY356',
#             'SDY357','SDY358','SDY479','SDY546','SDY557','SDY567','SDY571',
#             'SDY662','SDY668','SDY670','SDY671','SDY674','SDY689')

recip_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355',
            'SDY356','SDY357','SDY358','SDY479','SDY546','SDY557','SDY567',
            'SDY571','SDY662','SDY668','SDY670','SDY671','SDY674','SDY689')

# recip_kidney_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355','SDY356',
#                    'SDY357','SDY358','SDY546','SDY668','SDY671','SDY674')

## SDY546 and SDY674 are a mix - remember to remove non-donors
## be sure to check for duplicated lines in subjects! (they can be different line due to differences in ARMs)
# donor_ids=c('SDY289','SDY290','SDY291','SDY292','SDY293','SDY294','SDY546','SDY674')


## where the serialized data is locally
# studies_dir = "/Users/jiemingchen/Documents/transplantation/data"
# 
# #### starting up, download and serialize specific studies - do once ####
# ## set up tunnel
system("/Users/jiemingchen/.start_tunnel.sh")
# 
# ## a script with db connect
source("/Users/jiemingchen/R_codes/immport_connect.R")
# 
# ## serialize specific studies on local directory
# serialzeStudyData(study_ids, data_dir=studies_dir)
# 
# ## close tunnel
system("/Users/jiemingchen/.stop_tunnel.sh")
# 
# 
# #### studies stats ####
# 
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

#### load serialized data and write to file ####
#### loop all domains from RImmPort
domains = getListOfDomains()

## donors 
# for (i in 1:nrow(domains))
# {
#   data.donor = loadMultipleStudies(donor_ids,studies_dir,domains$`Domain Name`[i])
#   write.table(data.donor, paste("/Users/jiemingchen/Documents/transplantation/",
#                                 domains$`Domain Code`[i],"_immport2tsv_donors_subjects_8studies.tsv",
#                                 sep=""), 
#               sep="\t", quote = FALSE, row.names = FALSE)
# }

## recipients
for (i in 1:nrow(domains))
{
  data.donor = loadMultipleStudies(recip_ids,studies_dir,domains$`Domain Name`[i])
  write.table(data.donor, paste("/Users/jiemingchen/Documents/transplantation/",
                                domains$`Domain Code`[i],"_immport2tsv_recip_subjects_21studies.tsv",
                                sep=""), 
              sep="\t", quote = FALSE, row.names = FALSE)
}