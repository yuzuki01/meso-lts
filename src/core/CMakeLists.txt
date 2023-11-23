######################
# CMakeFile for core #
######################
cmake_minimum_required(VERSION 3.17)
project(core)
set(CMAKE_BUILD_TYPE MinSizeRel)

# OS type
if(WIN32)
    # if Windows
    message("OS is Windows")
    set(SHARED_LIBRARY_SUFFIX ".dll")
elseif(UNIX AND NOT APPLE)
    # if Linux
    message("OS is Linux")
    set(SHARED_LIBRARY_SUFFIX ".so")
else()
    message(FATAL_ERROR "unsupported OS.")
endif()

set(PROJECT_HOME_PATH ${PROJECT_SOURCE_DIR}/../..)

file(GLOB SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mesomath/*.cpp
        )
include_directories(${PROJECT_HOME_PATH}/src/include)

set(LIBRARY_OUTPUT_PATH ${PROJECT_HOME_PATH}/build)     # output .dll file to ~/build
link_directories(${LIBRARY_OUTPUT_PATH})

set(CMAKE_SHARED_LIBRARY_PREFIX "")                     # remove prefix
add_library(${PROJECT_NAME}                             # make dll
        SHARED ${SRC})
