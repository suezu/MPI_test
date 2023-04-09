cmake_minimum_required(VERSION 3.10)

project(Nparticle)

#set(CMAKE_CXX_STANDARD 11)

include_directories(include)

add_executable(excute src/2D_vector.cpp)


find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${MPI_INCLUDE_PATH})
