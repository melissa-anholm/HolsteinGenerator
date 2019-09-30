cmake_minimum_required (VERSION 2.6)
project (question_07 CXX) # This is where it looks for the compiler in env. vars CC and/or CXX.  Anything set afterwards doesn't matter.

#SET(CMAKE_CXX_COMPILER $ENV{GCC_DIR}/g++)  # apparently this has to go before "project(...)" otherwise it doesn't work.
#SET(CMAKE_C_COMPILER $ENV{GCC_DIR}/gcc)
message("* This is a message.  CMAKE_CXX_COMPILER = " ${CMAKE_CXX_COMPILER})

#ROOT includes =========================
#the FindROOT.cmake script does all of the magic for us to let us link to all of the root libraries
#INCLUDE(${CMAKE_SOURCE_DIR}/FindROOT.cmake)  # This will look for findROOT.cmake right in the current directory.
INCLUDE($ENV{MISC_INCLUDE_DIR}/FindROOT.cmake)  # because findroot.cmake now lives in *this* directory.

INCLUDE($ENV{G4_DIR}/lib/Geant4-10.5.1/Geant4Config.cmake)  # G4 Includes Part 1
INCLUDE($ENV{G4_DIR}/lib/Geant4-10.5.1/UseGeant4.cmake)     # G4 Includes Part 2.  This is the ${Geant4_USE_FILE}.
message("* This is a message.  Geant4_LIBRARIES = " ${Geant4_LIBRARIES} )
message("* This is a message.  Geant4_USE_FILE = " ${Geant4_USE_FILE} )


INCLUDE_DIRECTORIES (${ROOT_INCLUDE_DIR})
#INCLUDE_DIRECTORIES ($ENV{G4_INCLUDE_DIR})  # The 'ENV' *must* be included, and it still only finds files in "" style quotes.  To do it with <> style quotes, this directory must *also* be in the path.
INCLUDE_DIRECTORIES ($ENV{MISC_INCLUDE_DIR})  # Do I need this for anything?  maybe not right now.

message("* This is a message.  ROOT_INCLUDE_DIR = " ${ROOT_INCLUDE_DIR})
#message("* This is a message.  ROOT_INCLUDE_DIR = " $ENV{ROOT_INCLUDE_DIR})

LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})
LINK_DIRECTORIES (/usr/local/lib)  # for trinatdaq

#Compiler Flags=========================
#the -pthread -m64 are to make root happy
#the -g is for debugging symbols. -O2 is for

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -pthread -m64 -std=c++11")
#SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -Wall -pthread -m64 -stdlib=libc++ -std=c++11") # -stdlib=libc++ is for llvm not g++.
#SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -Wall -pthread -m64")
#SET(CMAKE_BUILD_TYPE ReLWithDebInfo)  # this makes the .dSYM file.  idk why this happens.  

#Linking the Executable ================
add_executable (holstein TheHolstein.cpp)
#add_executable (holsteinview holsteinview.cpp)
#add_executable (holsteinview roc_compare.cpp)

#link this library with everything that is needed
#to pop up a windows with the plots at the end
target_link_libraries (holstein
                      ${Geant4_LIBRARIES}
                      Gui
                      Core
                     # Cint #doesn't work with rootv6.  I can comment it out, but then does root still work?!
                      RIO
                      Net
                      Hist
                      Graf
                      Graf3d
                      Gpad
                      Tree
                      Rint
                      Postscript
                      Matrix
                      Physics
                      MathCore
                      Thread
                      pthread
                      m
                      dl
                      )   ## does this even work?!

