cmake_minimum_required(VERSION 3.25)

project(meso)

set(PROJECT_VERSION 1.1.0)
set(CMAKE_CXX_STANDARD 17)
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

# head file
include_directories(./src/include)

# .exe output dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)

add_executable(${PROJECT_NAME}
        main.cpp
        main.h
        )

# .dll
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build) # output .dll to ~/build
link_directories(${LIBRARY_OUTPUT_PATH})
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/core${SHARED_LIBRARY_SUFFIX})  # link core.dll
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/mesh${SHARED_LIBRARY_SUFFIX})  # link mesh.dll
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/solver${SHARED_LIBRARY_SUFFIX})  # link solver.dll
target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/api${SHARED_LIBRARY_SUFFIX})  # link api.dll

# check openMp
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP Found.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
