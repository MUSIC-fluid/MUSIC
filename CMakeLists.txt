cmake_minimum_required(VERSION 3.4)
project(MUSIChydro VERSION 2.0 LANGUAGES CXX)

set(libname "music")
set(exename "MUSIChydro")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#Version number
add_custom_target(
    version_target ALL #Always runs
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/__version.h
)

set(tag_def "\\#define TAG ")
set(url_def "\\#define URL ")
set(branch_def "\\#define BRANCH ")
add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/__version.h # fake! ensure we run!
    COMMAND echo ${tag_def} \\\"`git describe --tag`\\\" > ${PROJECT_SOURCE_DIR}/src/version.h
    COMMAND echo ${url_def} \\\"`git remote get-url origin`\\\" >> ${PROJECT_SOURCE_DIR}/src/version.h
    COMMAND echo ${branch_def} \\\"`git rev-parse --abbrev-ref HEAD`\\\" >> ${PROJECT_SOURCE_DIR}/src/version.h
)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

set(CMAKE_INSTALL_RPATH ${CMAKE_HOME_DIRECTORY})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Find ROOT for this project - Not mandatory
find_package(ROOT 6)
if (ROOT_FOUND)
    message("Found ROOT. Will enable ROOT integration.")
    set(CMAKE_CXX_STANDARD 17) # ROOT requires c++17
    add_definitions( -DROOT_FOUND )
endif()

# Find GSL and openMP for this project:
find_package(GSL)
if (GSL_FOUND)
    message("Found GSL library ${GSL_INCLUDE_DIR}, will build MUSIC with GSL")
    include_directories(${GSL_INCLUDE_DIR})
    set (ENV{CPATH} ${GSL_INCLUDE_DIR})
elseif()
    message("Did not find GSL library, will build MUSIC alone")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using Clang compiler without OpenMP parallelization... ")
else ()
    find_package(OpenMP)
    if (OpenMP_FOUND)
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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    if (KNL)
        set(CMAKE_CXX_FLAGS "-vecabi=cmdtarget -debug inline-debug-info -fimf-use-svml -xmic-avx512 -qopt-report=5 -qopt-report-phase:vec ${OpenMP_CXX_FLAGS} -std=c++11 -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    else ()
        set(CMAKE_CXX_FLAGS "-ipo -qopt-report=5 -qopt-report-phase:vec ${OpenMP_CXX_FLAGS} -std=c++11 -fimf-use-svml -fma -align -finline-functions -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "-std=c++11  -stdlib=libc++ -DDOCTEST_CONFIG_DISABLE -DNDEBUG")
endif()

if (unittest)
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11 -DDOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")
endif()

string(APPEND CMAKE_CXX_FLAGS " -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")

add_subdirectory (src)
