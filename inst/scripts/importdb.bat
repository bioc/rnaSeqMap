
@ECHO OFF
ECHO.
ECHO ****************************** README ************************************
ECHO This Batch file should import extensions for the xmapcore_homo_sapiens_58 database
ECHO into your MySQL database.  
ECHO.
ECHO Then, change the following settings to match those used by your system,
ECHO and run this batch file.
ECHO **************************************************************************
ECHO.

REM ************************* USER DEFINED VARIABLES *************************

REM ** Path to mysql program (if not on %PATH% already)
SET MYS=mysql


REM ** The server your MySQL database is on
SET SVR=localhost

REM ** The username you connect to the db using
SET USR=CHANGEME

REM ** The password you use to connect to the db:
REM **     Use -p (default) to ask for password with each step
REM **     Use a blank value for no password
REM **     Use --password=YOURPASSWORD to automatically send your password
SET PWD=-p

REM ********************* END OF USER DEFINED VARIABLES **********************

IF NOT %USR%==CHANGEME GOTO OK
ECHO **ERROR** You need to edit this script, and define some variables before running it
ECHO.
GOTO END

:OK

%MYS%    -h %SVR% -u %USR% %PWD% xmapcore_homo_sapiens_58 < rna_seq_map.sql
%MYS%    -h %SVR% -u %USR% %PWD% xmapcore_homo_sapiens_58 < bio_sample.sql
%MYS%    -h %SVR% -u %USR% %PWD% xmapcore_homo_sapiens_58 < seq_read.sql
%MYS%    -h %SVR% -u %USR% %PWD% xmapcore_homo_sapiens_58 < procedures.sql
%MYS%    -h %SVR% -u %USR% %PWD% xmapcore_homo_sapiens_58 < indeces.sql


ECHO ****************************** README ************************************
ECHO The data has now been added to your xmapcore_homo_sapiens_58 database
ECHO running in MySQL on %SVR%.
ECHO.
ECHO Please check that you have no errors, warnings or skipped records in the
ECHO output above, as there is no way to test from this script
ECHO.
ECHO Have fun with the database!
ECHO **************************************************************************

:END


