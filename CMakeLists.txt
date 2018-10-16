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

option(WITH_EXAMPLES "Build Examples" OFF)
if(WITH_EXAMPLES)

    file(GLOB EXAMPLE_SRC_FILES ${PROJECT_SOURCE_DIR}/examples/*.cc)
    find_package( OpenCV REQUIRED )
    ## from list of files we'll create examples example_name.cpp -> example_name
    foreach(_example_file ${EXAMPLE_SRC_FILES})
        get_filename_component(_example_name ${_example_file} NAME_WE)
        add_executable(${_example_name} ${_example_file})
        target_include_directories(${_example_name} PRIVATE inc)
        target_link_libraries(${_example_name} marlin  ${OpenCV_LIBS})
    endforeach()
endif()

################################
# Testing

option(WITH_TESTS "Build Tests" OFF)
if(WITH_TESTS)
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
endif()

