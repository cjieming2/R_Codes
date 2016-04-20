# stop ssh tunnel to buttelab-s01.ucsf.edu server
# clean up
 
cat /tmp/tunnel_pid.txt | xargs kill
rm /tmp/tunnel_pid.txt
echo "tunnel stopped"
