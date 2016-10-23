#!/bin/bash

##for emailing message at completion:
LOGNAME="eDNA_Pipeline"
MAILTO="invertdna@gmail.com"
SUBJECT="Pipeline completed!"
ATTFILE=$(echo $(ls -t /Users/rpk/Desktop/*.csv | head -1))

echo -e "From: $LOGNAME\nTo: $MAILTO\nSubject: $SUBJECT\n\
Mime-Version: 1.0\nContent-Type: text/plain\n\nSample_Directory_contents:\n\n\n$(ls -lh $(echo $DIRECTORIES |cut -d' ' -f1))\n\n\nTop_Ten_Annotations\n" > /tmp/file
cat "/tmp/temp.txt" $(echo -e "\n\n") $ATTFILE>> /tmp/file
/usr/sbin/sendmail -t -oi < /tmp/file
