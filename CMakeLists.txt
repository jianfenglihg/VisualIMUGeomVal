cmake_minimum_required(VERSION 3.10)
project(Test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)

include_directories(
    ${EIGEN3_INCLUDE_DIR}
)
add_library(TestLib SHARED
    TestGeometry.cpp
    )
add_executable(test test.cpp)

target_link_libraries(test TestLib)
