cmake_minimum_required(VERSION 3.6)
project(Kavosh)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftree-vectorize -pthread")

set(SOURCE_FILES
    include/GetTimeMs64.hpp
    include/GraphData.hpp
    src/GetTimeMs64.cpp
    src/GraphData.cpp
    src/main.cpp
    Snap-2.4/snap-core/Snap.o)

include_directories(include
        Snap-2.4/snap-core
        Snap-2.4/glib-core
        nauty25r9)

link_directories(nauty25r9/)

link_libraries(nauty)

add_executable(Kavosh ${SOURCE_FILES})