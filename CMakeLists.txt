cmake_minimum_required(VERSION 3.5)
project(Marlin VERSION 0.01)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing=1 -Wswitch-enum -Wundef -pedantic  -Wfatal-errors -Wshadow")

################################
# Marlin library

add_library(marlin src/configuration.cc src/dictionary.cc src/marlin.cc src/compress.cc src/decompress.cc)
set_target_properties(marlin PROPERTIES
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER inc/marlin.h)
target_include_directories(marlin PRIVATE inc)

################################
# Samples

#add_executable(CompressImages images.cpp)
#target_link_libraries(project1 project1_lib)

################################
# Testing
enable_testing()

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cc)

## from list of files we'll create tests test_name.cpp -> test_name
foreach(_test_file ${TEST_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    target_include_directories(${_test_name} PRIVATE inc)
    target_link_libraries(${_test_name} marlin)
    add_test(${_test_name} ${_test_name})
endforeach()
