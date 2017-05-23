cmake_minimum_required(VERSION 2.8)
project(mpihydro)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require GSL and openMP for this project:
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
set(CMAKE_CXX_FLAGS "-O3 -Wall ${OpenMP_CXX_FLAGS}")

add_subdirectory (src)  
