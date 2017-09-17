cmake_minimum_required(VERSION 3.6)
project(Kavosh)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -std=c++11 -ftree-vectorize -pthread")

set(SOURCE_FILES
        include/GetTimeMs64.hpp
        include/GraphData.hpp
        src/GetTimeMs64.cpp
        src/GraphData.cpp
        src/main.cpp)

include_directories(include
        snap/snap-core
        snap/glib-core
        nauty26r7)

add_executable(Kavosh ${SOURCE_FILES})