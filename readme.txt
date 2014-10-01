Repository is handled with git

After copying file, some may be in at different line endings. In order to handle issue,
clear cached files and restore to head by which git copies files to working folder in correct line endings.
        git rm --cached -r .
        git reset --hard

To compile in Visual Studio
- install visual studio 2010 express
- in project directory type
        make

To compile in Debian
- install make, g++ and libpcap (wireshark)
        make clean : to clean old data
        make		

To compile in FreeBSD
- install gcc, and libpcap
        gmake clean : to clean old data
        gmake

To compile in CentOs
- install git, and libpcap-devel  (yum install libpcap-devel)
- install developer package (yum groupinstall "Development Tools")
        make clean : to clean old data
        make
