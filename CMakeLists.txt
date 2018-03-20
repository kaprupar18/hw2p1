cmake_minimum_required (VERSION 2.8.0)
project (hw2p1)


set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Werror -rdynamic -O3 -g -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")


find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(HEADERS
    common.h
    
)

set(SOURCES
    common.cpp
)

add_executable(serial_h2p1 ${HEADERS} serial ${SOURCES})
target_link_libraries( serial_h2p1)

add_executable(autograder_h2p1 ${HEADERS} autograder.cpp ${SOURCES})
target_link_libraries( autograder_h2p1)

add_executable(openmp_h2p1 ${HEADERS} openmp.cpp ${SOURCES})
target_link_libraries( openmp_h2p1)

add_executable(mpi_h2p1 ${HEADERS} mpi.cpp ${SOURCES})
target_link_libraries( mpi_h2p1 ${MPI_C_LIBRARIES})
