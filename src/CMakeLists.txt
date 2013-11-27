#================================================================================
# Fluid Dynamics Engine CMake file
#================================================================================

# cmake_minimum_required(VERSION 2.6)
# set(CMAKE_BUILD_TYPE Debug)
if("${CMAKE_BUILD_TYPE}" EQUAL "Debug")
  set(CMAKE_VERBOSE_MAKEFILE TRUE) # this does not even work...
endif()

set(BOOST_LIBS thread date_time system program_options)
find_package(Boost COMPONENTS ${BOOST_LIBS} REQUIRED)

# find_package(GSL REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)

if(APPLE AND UNIX)
  set(LIB_PATH "/usr/local/lib")
  set(INC_PATH "/usr/local/include")
elseif(UNIX AND NOT APPLE)
  set(LIB_PATH "/usr/local/lib")
  set(INC_PATH "/usr/local/include")
elseif(WIN32)
  set(LIB_PATH "C:/msys/1.0/local/lib")
  set(INC_PATH "C:/msys/1.0/local/include")
endif(APPLE AND UNIX)

#--------------------------------------------------------------------------------
# This is the list of source files that need to be compiled
#--------------------------------------------------------------------------------
set( fdl_SRCS
  core/main.cpp
  core/fluidsolver.cpp
  core/particlesystem.cpp
  io/exporterbase.cpp
  io/pngexporter.cpp
  io/df3exporter.cpp
  io/pbrtexporter.cpp
  logger/logger.cpp
  logger/logwriter.cpp
  logger/stdiowriter.cpp
  logger/syslogwriter.cpp
  render/glutapp.cpp
  render/camera.cpp
)

#--------------------------------------------------------------------------------
# Another list, this time it includes all header files that should be treated with moc
#--------------------------------------------------------------------------------
set( fdl_MOC_HDRS
  # application.h
)

#--------------------------------------------------------------------------------
# .ui files that are used in the project (none)
#--------------------------------------------------------------------------------
set( fdl_UIS
  # none yet.
)

#--------------------------------------------------------------------------------
# and finally any resource files (none)
#--------------------------------------------------------------------------------
set( fdl_RCS
  # none yet.
)

#--------------------------------------------------------------------------------
# This sets the windows build that will need the special winmain@16 call. Qt provides
#  this for us in the qtmain.lib file. Using this cmake code will ensure we have it
#  linked into our build. Not needed on Unix/OS X/Linux which is why we have the 
#  IF(WIN32) conditional.
if (WIN32)
  set (QT_USE_QTMAIN TRUE)
endif (WIN32)

#--------------------------------------------------------------------------------
# This command finds Qt4 libraries and sets all required variables
# find_package( Qt4 REQUIRED )

#--------------------------------------------------------------------------------
# QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that contains a path to
# a CMake script that sets up this directory for building with Qt4
# include( ${QT_USE_FILE} )

#--------------------------------------------------------------------------------
# This command will generate rules that will run rcc on all files from fdl_RCS
# in result fdl_GEN_RCS variable will contain paths to files produced by rcc
# qt4_add_resources( fdl_GEN_RCS ${fdl_RCS} )

#--------------------------------------------------------------------------------
# This will run uic on .ui files:
# qt4_wrap_ui( fdl_GEN_UIS ${fdl_UIS} )

#--------------------------------------------------------------------------------
# And finally this will run moc on the headers
# qt4_wrap_cpp( fdl_GEN_MOC_HDRS ${fdl_MOC_HDRS} )


#--------------------------------------------------------------------------------
# Use whatever QT executables were used for building VTK (not totally necessary)
# set(QT_QMAKE_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
# set(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
# set(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
# set(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")



#--------------------------------------------------------------------------------
# We need this to be able to include headers produced by uic in our code
# CMAKE_BINARY_DIR holds a path to the build directory, while INCLUDE_DIRECTORIES()
#  works just like INCLUDEPATH from qmake
include_directories (
  # Look in the cmake build directory (some generated headers could be there)
  ${INC_PATH}
  ${ZLIB_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR} 
  ${PNG_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/extern
  ${PROJECT_SOURCE_DIR}/include
)

#--------------------------------------------------------------------------------
# For Apple set the icns file containing icons
if(APPLE)
  # set how it shows up in the Info.plist file
  # set(MACOSX_BUNDLE_ICON_FILE fdl.icns)
  # set where in the bundle to put the icns file
  # set_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/fdl.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  # include the icns file in the target
  # set(fdl_SRCS ${fdl_SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/fdl.icns)
endif(APPLE)

#--------------------------------------------------------------------------------
#  CMake's way of creating an executable
# set( EXECUTABLE_OUTPUT_PATH ../bin )

# how can I setup an automated versioning that can be reflected in source code?
set( TARGET_VERSION_MAJOR 0 )
set( TARGET_VERSION_MINOR 2 )
add_definitions( -DCMAKE_TARGET_VERSION=1 -DTARGET_VERSION_MAJOR=${TARGET_VERSION_MAJOR} -DTARGET_VERSION_MINOR=${TARGET_VERSION_MINOR})

# ADD_EXECUTABLE( fdl MACOSX_BUNDLE WIN32
add_executable( fdl
  ${fdl_SRCS}
  ${fdl_GEN_MOC_HDRS}
  ${fdl_GEN_UIS}
  ${fdl_GEN_RCS}
)


#--------------------------------------------------------------------------------
# Tell CMake which libraries we need to link our executable against.
target_link_libraries ( fdl
  ${Boost_LIBRARIES}
  ${PNG_LIBRARY}
  ${ZLIB_LIBRARIES}
  ${GLUT_LIBRARY}
  ${OPENGL_LIBRARY}
  ${Ilmbase_HALF_LIBRARY}
  ${Ilmbase_IEX_LIBRARY}
  # Ilmbase_ILMTHREAD_LIBRARY
  ${Ilmbase_IMATH_LIBRARY}
)