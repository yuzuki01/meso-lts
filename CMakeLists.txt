cmake_minimum_required(VERSION 3.25)

project(meso-mpi)

set(CMAKE_CXX_STANDARD 17)
set(PROJECT_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})

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

set(LIBRARY_OUTPUT_PATH "build")
set(EXECUTABLE_OUTPUT_PATH "build")

include_directories("src/include")

file(GLOB SRC
        "src/*/*.cpp"
        "src/*/*/*.cpp"
)

# set(CMAKE_SHARED_LIBRARY_PREFIX "")                     # remove prefix

# 根据不同的目标执行不同的指令
add_executable(${PROJECT_NAME}
        ${SRC}
)

# METIS
target_link_libraries(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib/libmetis.a")

# check MPI
find_package(MPI REQUIRED)
if (MPI_FOUND)
    target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_INCLUDE_PATH})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${MPI_LIBRARIES})
endif ()
