cmake_minimum_required(VERSION 3.16)
project(rans2d CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -O2 -march=native)

include_directories(include)

add_executable(rans2d
    src/main.cpp
    src/Mesh.cpp
    src/Solver.cpp
    src/Flux.cpp
    src/BC.cpp
    src/IO.cpp
)
