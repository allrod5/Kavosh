cmake_minimum_required(VERSION 3.6)
project(Kavosh)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -std=c++11 -ftree-vectorize -pthread -DDEBUG")

set(SOURCE_FILES
        include/GetTimeMs64.hpp
        include/GraphData.hpp
        src/GetTimeMs64.cpp
        src/GraphData.cpp
        src/main.cpp
        Snap-2.4/snap-core/Snap.o
        nauty26r7/nauty.o)

include_directories(include
        Snap-2.4/snap-core
        Snap-2.4/glib-core
        nauty26r7)

link_directories(nauty26r7/)

add_library(nauty STATIC IMPORTED)

set_target_properties(nauty PROPERTIES
        LINKER_LANGUAGE CXX
        IMPORTED_LOCATION /home/rodrigo/Documentos/UFABC/Iniciação\ Científica/PIBIC\ -\ 01.2015/Software/Kavosh/nauty26r7/nauty.a)

link_libraries(nauty)

add_executable(Kavosh ${SOURCE_FILES})