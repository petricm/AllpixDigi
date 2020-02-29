cmake_minimum_required(VERSION 3.14)
# project name
PROJECT( AllpixDigi )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )
INCLUDE( CheckIncludeFileCXX )

FIND_PACKAGE(ROOT 6.08 REQUIRED)
LINK_LIBRARIES( ${ROOT_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} )
SET( CMAKE_REQUIRED_INCLUDES ${ROOT_INCLUDE_DIRS} )
CHECK_INCLUDE_FILE_CXX( "${ROOT_INCLUDE_DIRS}/vdt/atan.h"  ROOT_VDT_FOUND )

FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES( SYSTEM ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE(MarlinUtil 1.4 REQUIRED)
FIND_PACKAGE(Marlin 1.0 REQUIRED)
FIND_PACKAGE(MarlinTrk REQUIRED)
FIND_PACKAGE(AIDA REQUIRED)


### LIBRARY #################################################################

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# Prevent O3 flag from being overwritten by something exporting O2 or worse
ADD_DEFINITIONS(" -O3 ")

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# optional package
FOREACH( pkg Marlin MarlinTrk MarlinUtil AIDA)
	IF( ${pkg}_FOUND )
		INCLUDE_DIRECTORIES(SYSTEM ${${pkg}_INCLUDE_DIRS})
		TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${${pkg}_LIBRARIES})
		ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
		MESSAGE( STATUS "${pkg} -- found" )
	ELSE()
    MESSAGE( STATUS "${pkg} -- not found" )
	ENDIF()
ENDFOREACH()


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()