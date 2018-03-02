cmake_minimum_required(VERSION 2.8)
project(mpihydro)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require GSL and openMP for this project:
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIR})
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})

option (KNL "Build executable on KNL" ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (KNL)
        set(CMAKE_CXX_FLAGS "-g -dynamic  -debug inline-debug-info -xmic-avx512 -qopt-report=5 -qopt-report-phase:vec -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    else ()
        set(CMAKE_CXX_FLAGS "-g -ipo -qopt-report=5 -qopt-report-phase:vec -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -fimf-use-svml -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-g -flto -fsanitize=address -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif()

add_subdirectory (src)  
#set_property(TARGET mpihydro PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

#set(CMAKE_CXX_FLAGS "-g -qopt-report=5 -qopt-report-phase:vec -ipo -Ofast -Wall ${OpenMP_CXX_FLAGS} -std=c++11 -xmic-avx512 -fma -align -finline-functions -DDOCTEST_CONFIG_DISA
