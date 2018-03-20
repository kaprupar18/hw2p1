cmake_minimum_required (VERSION 2.8.0)
project (hw2p1)

set(hw2p1_VERSION "v0.1")
set(hw2p1_VERSION_QUALIFIERS "")
string(TIMESTAMP hw2p1_VERSION_DATE %Y-%m-%d\ %H%M UTC)
set(hw2p1_COMPLETE_VERSION_STRING "${hw2p1_VERSION} ${hw2p1_VERSION_QUALIFIERS} (${hw2p1_VERSION_DATE})")
message(STATUS "Using version string: ${hw2p1_COMPLETE_VERSION_STRING}")

set(CMAKE_VERBOSE_MAKEFILE off)
if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Werror -rdynamic -O2 -g")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
endif(CMAKE_COMPILER_IS_GNUCC)



set(HEADERS
    common.h
)

set(SOURCES
    common.cpp
)

add_executable(hw2p1_engine ${HEADERS} main.cpp ${SOURCES})
target_compile_definitions( hw2p1_engine
                            PUBLIC -Dhw2p1_VERSION_STRING="${hw2p1_COMPLETE_VERSION_STRING}" )
target_link_libraries( hw2p1_engine)
