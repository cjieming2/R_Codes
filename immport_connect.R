source("/Users/jiemingchen/R_codes/jmRlib.R")

## try http:// if https:// URLs are not supported
## import ImmPort ; doesnt work on 4/12/16
#source("https://bioconductor.org/biocLite.R")
#biocLite("RImmPort")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(RMySQL)

# get the directory of a sample SQLite database that has been 
# bundled into the RImmPort package
# db_dir <- system.file("extdata", "ImmPortStudies", "Db", package = "RImmPort")

# connect to the database
# sqlite_conn <- dbConnect(SQLite(), dbname=file.path(db_dir, "ImmPort.sqlite"))

# set the ImmPort SQLite database as the ImmPort data source
# setImmPortDataSource(sqlite_conn)

## provide appropriate MySQL database connection parameters
## requires a tunnel that opens a port 3307 for redirect to buttelab server port 3306, which connects to AWS 
## run the following on bash to create this tunnel in the background
## ssh chenj@buttelab-s01.ucsf.edu -fNg -L 3307:buttelab-aws-rds01.cd8zgucpvgtu.us-west-2.rds.amazonaws.com:3306
mysql_conn <- dbConnect(MySQL(), user="chenj", password="3VrTh60IlfiHjLATiVkKn8orM",
                        dbname="proj_study_ALLSTUDIES",
                        host="127.0.0.1", port=3307)

## Chethan used this; the host is also on AWS but it's public (doesnt require credentials) 
## hence we are able to connect from here
## not the latest release
# mysql_conn  = dbConnect(MySQL(), user='lab_user', 
#                         password='big_data', 
#                         dbname='proj_study_ALLSTUDIES',
#                         host='exxport-mysql.celiaramktyo.us-west-2.rds.amazonaws.com')

# set the data source as the ImmPort MySQL database.
setImmPortDataSource(mysql_conn)
