cmake_minimum_required(VERSION 3.9)

project(compress)

set(CMAKE_CXX_STANDARD 11)

option(NORANDOM "use deterministic rounding for intml" OFF)
option(DEBUG "print some debug messages" OFF)

if(NORANDOM)
    message(STATUS "use deterministic rounding for intml")
    add_definitions(-DNORANDOM)
endif()

add_executable(compress main.cpp)
if(TARGET compress)
    target_include_directories(compress PUBLIC include)
endif()

if(DEBUG)
    add_definitions(-DDEBUG)
    message(STATUS "Debug build")
    string(APPEND CMAKE_CXX_FLAGS " -g -O0")
else()
    string(APPEND CMAKE_CXX_FLAGS " -O3")
endif()
