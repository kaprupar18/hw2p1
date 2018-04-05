cmake_minimum_required (VERSION 2.8.0)
project (hw2p1)


set(CMAKE_CXX_FLAGS "-std=c++11 -rdynamic -O3 -g -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")

configure_file(auto-bridges-mpi16 auto-bridges-mpi16 COPYONLY)
configure_file(auto-bridges-pthreads16 auto-bridges-pthreads16 COPYONLY)
configure_file(auto-bridges-openmp16 auto-bridges-openmp16 COPYONLY)
configure_file(auto-bridges-serial auto-bridges-serial COPYONLY)
configure_file(job-bridges-mpi16 job-bridges-mpi16 COPYONLY)
configure_file(job-bridges-openmp16 job-bridges-openmp16 COPYONLY)
configure_file(job-bridges-pthreads16 job-bridges-pthreads16 COPYONLY)
configure_file(job-bridges-serial job-bridges-serial COPYONLY)
configure_file(full-auto full-auto COPYONLY)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(HEADERS
    common.h
    
)

set(SOURCES
    #common.cpp
)

add_executable(serial ${HEADERS} serial common.cpp)
target_link_libraries( serial)

add_executable(autograder ${HEADERS} autograder.cpp common.cpp)
target_link_libraries( autograder)

add_executable(openmp ${HEADERS} openmp.cpp ${SOURCES})
target_link_libraries( openmp)

add_executable(mpi ${HEADERS} mpi.cpp ${SOURCES})
target_link_libraries( mpi ${MPI_C_LIBRARIES})
