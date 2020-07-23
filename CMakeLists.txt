CMAKE_MINIMUM_REQUIRED(VERSION 2.7)
PROJECT(holstein_generator CXX) # This is where it looks for the compiler in env. vars CC and/or CXX.  Anything set afterwards doesn't matter.

#SET(CMAKE_CXX_COMPILER $ENV{GCC_DIR}/g++)  # apparently this has to go before "project(...)" otherwise it doesn't work.
#SET(CMAKE_C_COMPILER $ENV{GCC_DIR}/gcc)
#message("* This is a message.  CMAKE_CXX_COMPILER = " ${CMAKE_CXX_COMPILER})

#ROOT includes =========================
#the FindROOT.cmake script does all of the magic for us to let us link to all of the root libraries
#INCLUDE(${CMAKE_SOURCE_DIR}/FindROOT.cmake)  # This will look for findROOT.cmake right in the current directory.
#INCLUDE($ENV{MISC_INCLUDE_DIR}/FindROOT.cmake)  # because findroot.cmake now lives in *this* directory.

#INCLUDE($ENV{G4_DIR}/lib/Geant4-10.5.1/Geant4Config.cmake)  # G4 Includes Part 1
#INCLUDE($ENV{G4_DIR}/lib/Geant4-10.5.1/UseGeant4.cmake)     # G4 Includes Part 2.  This is the ${Geant4_USE_FILE}.
#message("* This is a message.  Geant4_LIBRARIES = " ${Geant4_LIBRARIES} )
#message("* This is a message.  Geant4_USE_FILE = " ${Geant4_USE_FILE} )

FIND_PACKAGE(Geant4 REQUIRED ui_all vis_all)
#INCLUDE(FindROOT.cmake)
INCLUDE($ENV{MISC_INCLUDE_DIR}/FindROOT.cmake)  # because findroot.cmake now lives in *this* directory.
INCLUDE(${Geant4_USE_FILE})


INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(SYSTEM ${ROOT_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
INCLUDE_DIRECTORIES ($ENV{MISC_INCLUDE_DIR})  # Do I need this for anything?  maybe not right now.


#INCLUDE_DIRECTORIES (K37_include_dir)

SET(K37_include_dir /Users/anholm/Desktop/Trinat_Geant/include)
SET( shared_hh 
	 Holstein52Generator.hh
	 Holstein52Isotope.hh
	 IsotopeValues.hh
	 K37AtomicSetup.hh
	 K37Cloud.hh
	 K37ComplexGammaFunction.hh
	 K37FermiFunction.hh
	 K37SublevelPopulations.hh
	)
foreach(_file ${shared_hh})
	configure_file(
		${K37_include_dir}/${_file}
		${CMAKE_SOURCE_DIR}/include/${_file}
		COPYONLY
		)
endforeach()

SET(K37_src_dir /Users/anholm/Desktop/Trinat_Geant/src)
SET( shared_cc 
	 Holstein52Generator.cc
	 Holstein52Isotope.cc
	 IsotopeValues.cc
	 K37AtomicSetup.cc
	 K37Cloud.cc
	 K37ComplexGammaFunction.cc
	 K37FermiFunction.cc
	 K37SublevelPopulations.cc
	)
foreach(_file ${shared_cc})
	configure_file(
		${K37_src_dir}/${_file}
		${CMAKE_SOURCE_DIR}/src/${_file}
		COPYONLY
		)
endforeach()



SET(a_bunch_of_root_libraries
	Gui Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread pthread m dl )

#INCLUDE_DIRECTORIES ($ENV{G4_INCLUDE_DIR})  # The 'ENV' *must* be included, and it still only finds files in "" style quotes.  To do it with <> style quotes, this directory must *also* be in the path.
#message("* This is a message.  ROOT_INCLUDE_DIR = " ${ROOT_INCLUDE_DIR})
#message("* This is a message.  ROOT_INCLUDE_DIR = " $ENV{ROOT_INCLUDE_DIR})

LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})
#LINK_DIRECTORIES(K37_src_dir)

#LINK_DIRECTORIES (/usr/local/lib)  # for trinatdaq

#Compiler Flags=========================
#the -pthread -m64 are to make root happy
#the -g is for debugging symbols. -O2 is for

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -pthread -m64 -std=c++11")
#SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -Wall -pthread -m64 -stdlib=libc++ -std=c++11") # -stdlib=libc++ is for llvm not g++.
#SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -Wall -pthread -m64")
#SET(CMAKE_BUILD_TYPE ReLWithDebInfo)  # this makes the .dSYM file.  idk why this happens.  

#INCLUDE_DIRECTORIES(include)

SET (CONFIGURATION_DIRECTORY "${PROJECT_BINARY_DIR}/")  
SET (OUTDIR "\"${PROJECT_SOURCE_DIR}/Output\"")  # no trailing slash?
configure_file (
	"${PROJECT_SOURCE_DIR}/include/K37Config.hh.in"
	"${PROJECT_SOURCE_DIR}/include/K37Config.hh"
	)


ADD_SUBDIRECTORY(src)
ADD_EXECUTABLE(holstein TheHolstein.cpp)
  
target_link_libraries(holstein
	holstein_lib ${Geant4_LIBRARIES} ${a_bunch_of_root_libraries} )

#add_executable (holsteinview holsteinview.cpp)
#target_link_libraries(holsteinview
#	holstein_lib ${Geant4_LIBRARIES} ${a_bunch_of_root_libraries} )

#add_executable (holsteinview roc_compare.cpp)
#target_link_libraries(holsteinview
#	holstein_lib ${Geant4_LIBRARIES} ${a_bunch_of_root_libraries} )


###  Everything is compiled, so let's copy some files over.  ###
###  ===========================================================
SET(K37_shared_setup /Users/anholm/Desktop/Trinat_Geant/ConfigurationFiles)
SET(setup_files 
	K_37_INPUT.txt
	K_37_POPULATIONS_INPUT.txt
	)

foreach(_file ${setup_files})
	configure_file(
		${K37_shared_setup}/${_file}
		${CMAKE_BINARY_DIR}/${_file}
		COPYONLY
		)
endforeach()



#Linking the Executable ================
#add_executable (holstein TheHolstein.cpp)
#add_executable (holsteinview holsteinview.cpp)
#add_executable (holsteinview roc_compare.cpp)

#link this library with everything that is needed
#to pop up a windows with the plots at the end
#target_link_libraries (holstein
#                      holstein_lib
#                      ${Geant4_LIBRARIES}
#                      Gui
#                      Core
#                      RIO
#                      Net
#                      Hist
#                      Graf
#                      Graf3d
#                      Gpad
#                      Tree
#                      Rint
#                      Postscript
#                      Matrix
#                      Physics
#                      MathCore
#                      Thread
#                      pthread
#                      m
#                      dl
#                      )   ## 


#target_link_libraries (holsteinview
#                      ${Geant4_LIBRARIES}
#                      holstein_lib  ## needed??
#                      Gui
#                      Core
#                      RIO
#                      Net
#                      Hist
#                      Graf
#                      Graf3d
#                      Gpad
#                      Tree
#                      Rint
#                      Postscript
#                      Matrix
#                      Physics
#                      MathCore
#                      Thread
#                      pthread
#                      m
#                      dl
#                      )   ## 

