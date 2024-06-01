cmake_minimum_required(VERSION 3.25)

project(meso)

set(CMAKE_CXX_STANDARD 17)
set(PROJECT_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)

# OS type
if (WIN32)
    # if Windows
    message("OS is Windows")
    set(OS_INDICATOR "win")
    set(SHARED_LIBRARY_SUFFIX ".dll")
elseif (UNIX AND NOT APPLE)
    # if Linux
    message("OS is Linux")
    set(OS_INDICATOR "linux")
    set(SHARED_LIBRARY_SUFFIX ".so")
else ()
    message(FATAL_ERROR "unsupported OS.")
endif ()

set(LIBRARY_OUTPUT_PATH ${PROJECT_HOME_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_HOME_DIR}/bin)

include_directories(${PROJECT_HOME_DIR}/include)

file(GLOB SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
)

# set(CMAKE_SHARED_LIBRARY_PREFIX "")                     # remove prefix

# 根据不同的目标执行不同的指令
add_executable(${PROJECT_NAME}
        main.cpp
        ${SRC}
)

target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_HOME_DIR}/bin/libcore${SHARED_LIBRARY_SUFFIX})
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_HOME_DIR}/bin/libmesh${SHARED_LIBRARY_SUFFIX})
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_HOME_DIR}/bin/libsolver${SHARED_LIBRARY_SUFFIX})

# check openMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# check MPI
find_package(MPI REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_LIBRARIES})
