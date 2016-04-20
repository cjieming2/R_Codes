#!/bin/bash -eu
 
# start ssh tunnel to buttelab-s01.ucsf.edu server
# get access to aws mysql server from local system:
# mysql -h 127.0.0.1 -u MYSQL_USER_NAME -P 3307 -p  
 
ssh chenj@buttelab-s01.ucsf.edu -Ng -L 3307:buttelab-aws-rds01.cd8zgucpvgtu.us-west-2.rds.amazonaws.com:3306 &  echo $! > /tmp/tunnel_pid.txt
echo "tunnel started"
