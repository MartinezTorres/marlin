cmake_minimum_required(VERSION 3.5)
project(Marlin VERSION 0.01)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing=1 -Wswitch-enum -Wundef -pedantic -Wfatal-errors -Wshadow")

################################
# Marlin library (entropy codec only)
file(GLOB MAIN_SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB IMAGE_ONLY_SRC_FILES ${PROJECT_SOURCE_DIR}/src/image*.cc)
list(REMOVE_ITEM MAIN_SRC_FILES ${IMAGE_ONLY_SRC_FILES})
add_library(marlin ${MAIN_SRC_FILES})
set_target_properties(marlin PROPERTIES
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER inc/marlin.h)
target_include_directories(marlin PRIVATE inc)

# ImageMarlin library (image codec + entropy codec)
file(GLOB MAIN_SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_library(imarlin ${MAIN_SRC_FILES} inc/imageMarlin.hpp)
set_target_properties(imarlin PROPERTIES
        VERSION ${PROJECT_VERSION}
        PUBLIC_HEADER inc/imageMarlin.hpp)
target_include_directories(imarlin PRIVATE inc)

################################
# Samples

option(WITH_UTILS "Build Utilities" ON)
if(WITH_UTILS)
    file(GLOB UTILS_SRC_FILES ${PROJECT_SOURCE_DIR}/utils/*.cc)
    find_package( OpenCV REQUIRED )

    foreach(_util_file ${UTILS_SRC_FILES})
        get_filename_component(_util_name ${_util_file} NAME_WE)
        add_executable(${_util_name} ${_util_file})
        target_include_directories(${_util_name} PRIVATE inc utils)
#        target_link_libraries(${_util_name} marlin  ${OpenCV_LIBS} ${Boost_LIBRARIES})
        target_link_libraries(${_util_name} imarlin ${OpenCV_LIBS})
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

