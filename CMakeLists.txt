# CMake config for osmium-history-splitter
#
# by Matt Amos <zerebubuth@gmail.com>
# based on osmium-tool's CMakeLists.txt by Jochen Topf.

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# version

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;MinSizeRel;Dev"
  CACHE STRING
  "List of available configuration types"
  FORCE)

project(osmium-history-splitter)

set(OSMIUM_HISTORY_SPLITTER_VERSION_MAJOR 0)
set(OSMIUM_HISTORY_SPLITTER_VERSION_MINOR 0)
set(OSMIUM_HISTORY_SPLITTER_VERSION_PATCH 1)

set(OSMIUM_HISTORY_SPLITTER_VERSION "${OSMIUM_HISTORY_SPLITTER_VERSION_MAJOR}.${OSMIUM_HISTORY_SPLITTER_VERSION_MINOR}.${OSMIUM_HISTORY_SPLITTER_VERSION_PATCH}")

add_definitions(-DOSMIUM_HISTORY_SPLITTER_VERSION="${OSMIUM_HISTORY_SPLITTER_VERSION}")
set(AUTHOR "Matt Amos <zerebubuth@gmail.com>")

# external dependencies

find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

find_package(Osmium REQUIRED COMPONENTS io)
include_directories(${OSMIUM_INCLUDE_DIRS})

# optional cppcheck

message(STATUS "Looking for cppcheck")
find_program(CPPCHECK cppcheck)

if (CPPCHECK)
  message(STATUS "Looking for cppcheck - found")
  set(CPPCHECK_OPTIONS --enable=warning,style,performance,portability,information,missingInclude)

  # ???
  set(CPPCHECK_OPTIONS ${CPPCHECK_OPTIONS} --suppress=missingIncludeSystem)

  add_custom_target(cppcheck ${CPPCHECK} --std=c++11 ${CPPCHECK_OPTIONS} ${CMAKE_SOURCE_DIR}/src/*pp)
else()
  message(STATUS "Looking for cppcheck - not found")
  message(STATUS "  Target cppcheck not available.")
endif(CPPCHECK)

# use C++11, or user-supplied version

if(NOT USE_CPP_VERSION)
  set(USE_CPP_VERSION c++11)
endif()
add_definitions(-std=${USE_CPP_VERSION})

# compiler & linker flags

set(USUAL_COMPILE_OPTIONS "-g -ggdb -O0")

set(CMAKE_CXX_FLAGS_DEV "${USUAL_COMPILE_OPTIONS}"
  CACHE STRING "Flags used by the compiler during developer builds."
  FORCE)

set(CMAKE_EXE_LINKER_FLAGS_DEV ""
  CACHE STRING "Flags used by the linker during developer builds."
  FORCE)
mark_as_advanced(
  CMAKE_CXX_FLAGS_DEV
  CMAKE_EXE_LINKER_FLAGS_DEV
)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${USUAL_COMPILE_OPTIONS}"
  CACHE STRING "Flags used by the compiler during RELWITHDEBINFO builds."
  FORCE)

# build type options

if(CMAKE_BUILD_TYPE STREQUAL "Dev")
  add_definitions(-Werror)
  add_definitions(${OSMIUM_WARNING_OPTIONS})
endif()

if(CMAKE_BUILD_TYPE)
  set(build_type ${CMAKE_BUILD_TYPE})
else()
  set(build_type "RelWithDebInfo")
endif()

set(CMAKE_BUILD_TYPE ${build_type}
  CACHE STRING
  "Choose the type of build, options are: ${CMAKE_CONFIGURATION_TYPES}."
  FORCE)

# tests

enable_testing()
add_subdirectory(test)

# header files

include_directories(include)

# source files

add_subdirectory(src)

