cmake_minimum_required(VERSION 3.5)
project(Marlin VERSION 0.01)

set(CMAKE_CXX_STANDARD 14)

################################
# Marlin library

add_library(marlin src/marlin.cc src/dictionary.cc)
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

# from list of files we'll create tests test_name.cpp -> test_name
foreach(_test_file ${TEST_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    add_dependencies(${_test_name} "googletest.git")
    target_link_libraries(${_test_name} gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
    add_test(${_test_name} ${_test_name})
    set_tests_properties(${_test_name} PROPERTIES TIMEOUT 5)
endforeach()
