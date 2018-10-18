cmake_minimum_required(VERSION 3.5)
project(Marlin VERSION 0.01)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing=1 -Wswitch-enum -Wundef -pedantic  -Wfatal-errors -Wshadow")

################################
# Marlin library

file(GLOB MAIN_SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_library(marlin ${MAIN_SRC_FILES})
set_target_properties(marlin PROPERTIES
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER inc/marlin.h)
target_include_directories(marlin PRIVATE inc)

################################
# Samples

option(WITH_UTILS "Build Utilities" ON)
if(WITH_UTILS)

    file(GLOB UTILS_SRC_FILES ${PROJECT_SOURCE_DIR}/utils/*.cc)
    find_package( OpenCV REQUIRED )
#    find_package( Boost 1.40 COMPONENTS program_options REQUIRED )
## from list of files we'll create utils util_name.cpp -> util_name
    foreach(_util_file ${UTILS_SRC_FILES})
        get_filename_component(_util_name ${_util_file} NAME_WE)
        add_executable(${_util_name} ${_util_file})
        target_include_directories(${_util_name} PRIVATE inc utils)
#        target_link_libraries(${_util_name} marlin  ${OpenCV_LIBS} ${Boost_LIBRARIES})
        target_link_libraries(${_util_name} marlin  ${OpenCV_LIBS})
    endforeach()
endif()

################################
# Testing

option(WITH_TESTS "Build Tests" ON)
if(WITH_TESTS)
    enable_testing()
    file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cc)

    ## from list of files we'll create tests test_name.cpp -> test_name
    foreach(_test_file ${TEST_SRC_FILES})
        get_filename_component(_test_name ${_test_file} NAME_WE)
        add_executable(${_test_name} ${_test_file})
        target_include_directories(${_test_name} PRIVATE inc utils)
        target_link_libraries(${_test_name} marlin)
        add_test(${_test_name} ${_test_name})
    endforeach()
endif()

