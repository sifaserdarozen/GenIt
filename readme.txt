Repository is handled with git

If file-copied instead of using repository, some of the files may be at different line endings. In order to handle the issue, clear cached files and restore to head by which git corrects line endings.
        git rm --cached -r .
        git reset --hard

To compile in Visual Studio
- install visual studio 2010 express and cuda tools
- in project directory type
        make

To compile in Debian
- install make, g++ and cuda tools
        make clean : to clean old data
        make		

To compile in FreeBSD
- install gcc, and cuda tools (use pkg install)
        gmake clean : to clean old data
        gmake

To compile in CentOs
- install cuda tools  (use yum)
- install developer package (yum groupinstall "Development Tools")
        make clean : to clean old data
        make
