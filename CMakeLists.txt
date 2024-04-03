cmake_minimum_required(VERSION 3.5)
project(MUSIChydro VERSION 2.0 LANGUAGES CXX)

set(libname "music")
set(exename "MUSIChydro")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set(CMAKE_INSTALL_RPATH ${CMAKE_HOME_DIRECTORY})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Find GSL and openMP for this project:
find_package(GSL)
if (GSL_FOUND)
    message("Found GSL library ${GSL_INCLUDE_DIR}, will build MUSIC with GSL")
    include_directories(${GSL_INCLUDE_DIR})
elseif()
    message("Did not find GSL library, will build MUSIC alone")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler ... ")

    find_library(LIBOMP_LIBRARY NAMES libomp.dylib HINTS /opt/homebrew/opt/libomp/lib/ /opt/Homebrew/opt/libomp/lib/ /usr/local/opt/libomp/lib/ /usr/local/homebrew/opt/libomp/lib/ /usr/local/Homebrew/opt/libomp/lib/ NO_DEFAULT_PATH)

    if (LIBOMP_LIBRARY)

      set(OpenMP_CXX_FOUND TRUE)

      message("Found OpenMP with Clang ..." ${OpenMP_CXX_FOUND})
      find_path(OpenMP_incDir omp.h PATHS /opt/homebrew/opt/libomp/include /usr/local/homebrew/opt/libomp/include /opt/Homebrew/opt/libomp/include /usr/local/opt/libomp/include /usr/local/Homebrew/opt/libomp/include NO_DEFAULT_PATH)

      set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${OpenMP_incDir} -Wno-unused-command-line-argument")
      set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})

      message(STATUS "OpenMP_incDir: " ${OpenMP_incDir})
      message(STATUS "Found OpenMP library: " ${LIBOMP_LIBRARY})
      message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")

    endif ()
else ()
    find_package(OpenMP COMPONENTS CXX)
    if (OpenMP_CXX_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
    endif ()
endif()


# default build type: Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}/"
        CACHE PATH "Install path prefix, prepended onto install directories."
        FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

option (KNL "Build executable on KNL" OFF)
option (unittest "Build Unit tests" OFF)
option (link_with_lib "Link executable with the libarary" ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (KNL)
        set(CMAKE_CXX_FLAGS "-vecabi=cmdtarget -debug inline-debug-info -fimf-use-svml -xmic-avx512 -qopt-report=5 -qopt-report-phase:vec ${OpenMP_CXX_FLAGS} -std=c++11 -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    else ()
        set(CMAKE_CXX_FLAGS "-ipo -qopt-report=5 -qopt-report-phase:vec ${OpenMP_CXX_FLAGS} -std=c++11 -fimf-use-svml -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11  -stdlib=libc++ -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif()

if (unittest)
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")

find_package(Git)
if(Git_FOUND)
    # Get the current working branch
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest abbreviated commit hash of the working branch
    execute_process(
        COMMAND git log -1 --format=%H
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    add_definitions(-DGIT_COMMIT_HASH="${GIT_COMMIT_HASH}")
    add_definitions(-DGIT_BRANCH="${GIT_BRANCH}")
else ()
    add_definitions(-DGIT_COMMIT_HASH="unknown")
    add_definitions(-DGIT_BRANCH="unknown")
endif()

configure_file(
  src/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
)

include_directories(${CMAKE_BINARY_DIR}/generated)
add_subdirectory (src)
