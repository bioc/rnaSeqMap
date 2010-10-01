
#!/bin/bash

echo ""
echo "############################## README ####################################"
echo "This Batch file should import extensions for the xmapcore_homo_sapiens_58 database"
echo " into your MySQL database. "
echo ""
echo "Then, change the following settings to match those used by your system,"
echo "and run this batch file."
echo "##########################################################################"
echo ""

# Which server is the database on?
SVR=localhost

# What username should I use to connect?
USR=CHANGEME

# User password (uncomment the one you require):
# --Ask for password each time
PWD=-p
# --Use this password
# PWD=--password=PASSWORD
# --No password
# PWD=

# Locations of the mysql and mysqlimport commands (these should work as is if
# your PATH is set correctly)
MYS=mysql

if [ $USR = "CHANGEME" ]; then
  echo "##ERROR## You need to edit this script, and define some variables before running it"
  exit
fi

$MYS    -h $SVR -u $USR $PWD xmapcore_homo_sapiens_58 < rnaSeqMap.sql
$MYS    -h $SVR -u $USR $PWD xmapcore_homo_sapiens_58 < bio_sample.sql
$MYS    -h $SVR -u $USR $PWD xmapcore_homo_sapiens_58 < seq_read.sql
$MYS    -h $SVR -u $USR $PWD xmapcore_homo_sapiens_58 < procedures.sql
$MYS    -h $SVR -u $USR $PWD xmapcore_homo_sapiens_58 < indeces.sql

echo ""
echo "############################## README ####################################"
echo "The data has now been added to your xmapcore_homo_sapiens_58 database"
echo "running in MySQL on $SVR."
echo ""
echo "Please check that you have no errors, warnings or skipped records in the"
echo "output above, as there is no way to test from this script"
echo ""
echo "Have fun with the database!"
echo "##########################################################################"

