# Microsoft Developer Studio Project File - Name="WindShape" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=WindShape - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "WindShape.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "WindShape.mak" CFG="WindShape - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "WindShape - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "WindShape - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""$/WindWax/source", RIBAAAAA"
# PROP Scc_LocalPath "."
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "WindShape - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release\Temp"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "_AFXDLL" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x40c /d "NDEBUG"
# ADD RSC /l 0x40c /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 /nologo /subsystem:console /debug /machine:I386

!ELSEIF  "$(CFG)" == "WindShape - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug\Temp"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "_AFXDLL" /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x40c /d "_DEBUG"
# ADD RSC /l 0x40c /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 glui32d.lib glut32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /nodefaultlib:"msvcrt.lib" /out:"Release/WindShape.exe" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "WindShape - Win32 Release"
# Name "WindShape - Win32 Debug"
# Begin Source File

SOURCE=..\commun\fichier.cpp

!IF  "$(CFG)" == "WindShape - Win32 Release"

!ELSEIF  "$(CFG)" == "WindShape - Win32 Debug"

# ADD CPP /Gm- /GR-

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\commun\fichier.h
# End Source File
# Begin Source File

SOURCE=..\commun\geom.cpp
# End Source File
# Begin Source File

SOURCE=..\commun\geom.h
# End Source File
# Begin Source File

SOURCE=..\commun\matrice.cpp
# End Source File
# Begin Source File

SOURCE=..\commun\matrice.h
# End Source File
# Begin Source File

SOURCE=..\commun\pince.cpp
# End Source File
# Begin Source File

SOURCE=..\commun\plot.cpp
# End Source File
# Begin Source File

SOURCE=..\commun\plot.h
# End Source File
# Begin Source File

SOURCE=..\commun\profil.cpp
# End Source File
# Begin Source File

SOURCE=..\commun\profil.h
# End Source File
# Begin Source File

SOURCE=..\commun\rasklad.cpp
# End Source File
# Begin Source File

SOURCE=.\WindShape.cpp
# End Source File
# End Target
# End Project