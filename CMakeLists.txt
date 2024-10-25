cmake_minimum_required(VERSION 3.16)

# Name the project
project("PATTYMELT")

# Set up all root6 related macros and paths
find_package(ROOT REQUIRED)
# include(${ROOT_USE_FILE})
# include_directories(${ROOT_INCLUDE_DIRS})

# Define directories
set(INCLUDE_DIRECTORY $ENV{HOME}/include)
set(LIB_DIRECTORY $ENV{HOME}/lib)
set(SHARE_DIRECTORY $ENV{HOME}/share)

# Define experiment and calibration libraries for linking against
set(LIB_EXPT ${LIB_DIRECTORY}/lib063123.so)
# set(LIB_EXPTCALIB ${LIB_DIRECTORY}/lib112723Calib.so)

# Setting up experiment global include (to both targets and sources)
set(INCLUDE_EXPT ${INCLUDE_DIRECTORY}/063123)
include_directories(${INCLUDE_EXPT})

# Setting up CycApps global include (to both targets and sources)
exec_program(cycapp-config ARGS --incdir OUTPUT_VARIABLE CYCAPP_INCLUDE_DIRECTORY)
include_directories(${CYCAPP_INCLUDE_DIRECTORY})
exec_program(cycapp-config ARGS --libs OUTPUT_VARIABLE CYCAPP_LIBS)
set(CYCAPP_LIBS "${CYCAPP_LIBS} -lCycSrim")

exec_program(root-config ARGS --glibs OUTPUT_VARIABLE ROOT_GLIBS)

include("${ROOT_DIR}/modules/RootNewMacros.cmake")

# Define executable directory
set (INSTALL_EXEC_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Define project library
set(LIB_NAME PATTYMELT)
message("Entering src to create ${LIB_NAME} library")
add_subdirectory("src")