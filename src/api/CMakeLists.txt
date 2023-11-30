########################
# CMakeFile for solver #
########################
cmake_minimum_required(VERSION 3.17)
project(api)
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
        )
include_directories(${PROJECT_HOME_PATH}/src/include)

set(LIBRARY_OUTPUT_PATH ${PROJECT_HOME_PATH}/build)     # output .dll file to ~/build
link_directories(${LIBRARY_OUTPUT_PATH})

set(CMAKE_SHARED_LIBRARY_PREFIX "")                     # remove prefix
add_library(${PROJECT_NAME}                             # make dll
        SHARED ${SRC})
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/core${CMAKE_SHARED_LIBRARY_SUFFIX})    # link to core.dll
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/mesh${CMAKE_SHARED_LIBRARY_SUFFIX})    # link to mesh.dll
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/solver${CMAKE_SHARED_LIBRARY_SUFFIX})    # link to solver.dll

# check openMp
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
