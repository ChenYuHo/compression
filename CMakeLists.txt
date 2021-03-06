cmake_minimum_required(VERSION 3.9)

project(compress)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(NORANDOM "use deterministic rounding for intml" OFF)
option(DEBUG "print some debug messages" OFF)

if(NORANDOM)
    message(STATUS "use deterministic rounding for intml")
    add_definitions(-DNORANDOM)
endif()

#add_executable(test hls/compress_test.cpp)
add_executable(compress main.cpp)
if(TARGET compress)
    find_package(OpenMP)
    target_include_directories(compress PUBLIC include)
    include(CMakeHostSystemInformationExtra)
    include(GetCPUSIMDFeatures)
    include(CMakeCompilerMachineOption)
    set(vendor)
    set(architecture)
    detect_x64_micro_architecture(vendor architecture)
    message(STATUS "Found ${vendor}:${architecture}")
    if(OpenMP_CXX_FOUND)
        target_link_libraries(compress PUBLIC OpenMP::OpenMP_CXX)
    endif()
    cmake_compiler_machine_option(ARCHITECTURE_FLAG ${architecture})
    message(STATUS "Use compiler option: ${ARCHITECTURE_FLAG}")
    target_compile_options(compress PRIVATE ${ARCHITECTURE_FLAG})
endif()

if(DEBUG)
    add_definitions(-DDEBUG)
    message(STATUS "Debug build")
    string(APPEND CMAKE_CXX_FLAGS " -g -O0")
else()
    string(APPEND CMAKE_CXX_FLAGS " -O3")
endif()
