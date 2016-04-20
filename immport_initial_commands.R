## set up tunnel
system("/Users/jiemingchen/.start_tunnel.sh")

## provide appropriate MySQL database connection parameters
## requires a tunnel that opens a port 3307 for redirect to buttelab server port 3306, which connects to AWS 
## run the following on bash to create this tunnel in the background
## ssh chenj@buttelab-s01.ucsf.edu -fNg -L 3307:buttelab-aws-rds01.cd8zgucpvgtu.us-west-2.rds.amazonaws.com:3306
source('/Users/jiemingchen/R_codes/immport_connect.R')

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

## take in the study-id
study_id <- 'SDY514'
sdy514_data <- getStudy(study_id)

## access demographics data
dm_df <- sdy514_data$special_purpose$dm_l$dm_df
head(dm_df)

## access concomitant medications data
cm_df <- sdy514_data$interventions$cm_l$cm_df
head(cm_df)

## get Trial Title from Trial Summary
ts_df <- sdy514_data$trial_design$ts_l$ts_df
title <- ts_df$TSVAL[ts_df$TSPARMCD== "TITLE"]
title


## close tunnel
system("/Users/jiemingchen/.stop_tunnel.sh")
