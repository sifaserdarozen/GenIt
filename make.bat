@echo off
cls

echo ************************************************************
echo                   STARTING BUILD PROCESS
echo ************************************************************
echo .

REM go to library folder
echo Trying to create directory structure ...
mkdir lib
mkdir bin
cd lib
echo .

REM clear old content
del /Q *.*

REM compile cuda_coders.cu
nvcc --include-path ..\include -c -m32 --cl-version 2010 ..\src\cuda_coders.cu

REM compile cpu_coders.cpp
cl /Zi /I..\include /EHsc -c ..\src\cpu_coders.cpp 2>>null || echo COMPILE ERROR && exit /b

REM compile g722coder.cpp
cl /Zi /I..\include /EHsc -c ..\src\g722coder.cpp 2>>null || echo COMPILE ERROR && exit /b

REM building sip library
LIB *.obj /OUT:codecs.lib  1>>2>>null || echo BUILD BUILD ERROR && exit /b
REM clear object files
REM del /Q *.obj

echo .
echo .............Compilation of codes library is successful.............
echo .

REM copy trird party lib files to library folder
copy ..\3rdParty\lib\*.lib . 1>>null
REM build genit
cl /Zi /I..\include /EHsc codecs.lib cudart.lib ..\src\genit.cpp 2>>null || echo COMPILE ERROR && exit /b
move genit.exe ..\bin\genit_x86.exe 1>>null

echo .
echo .............Build of genit is successful.............
echo .

echo ************************************************************
echo                   BUILD PROCESS IS SUCCESSFULL
echo ************************************************************
echo .

cd ..








