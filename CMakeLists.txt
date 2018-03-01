cmake_minimum_required(VERSION 2.8)
project(mpihydro)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require GSL and openMP for this project:
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIR})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
set(CMAKE_CXX_FLAGS "-g -xmic-avx512 -ipo -qopt-report=5 -qopt-report-phase:vec -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set(CMAKE_CXX_FLAGS "-g -flto -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++14 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif()

add_subdirectory (src)  
#set_property(TARGET mpihydro PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

#set(CMAKE_CXX_FLAGS "-g -qopt-report=5 -qopt-report-phase:vec -ipo -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -xmic-avx512 -fma -align -finline-functions -DDOCTEST_CONFIG_DISA
