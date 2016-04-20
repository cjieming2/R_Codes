## set up tunnel
system("/Users/jiemingchen/.start_tunnel.sh")

## provide appropriate MySQL database connection parameters
## requires a tunnel that opens a port 3307 for redirect to buttelab server port 3306, which connects to AWS 
## run the following on bash to create this tunnel in the background
## ssh chenj@buttelab-s01.ucsf.edu -fNg -L 3307:buttelab-aws-rds01.cd8zgucpvgtu.us-west-2.rds.amazonaws.com:3306
source('/Users/jiemingchen/R_codes/immport_connect.R')

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")



## close tunnel
system("/Users/jiemingchen/.stop_tunnel.sh")
