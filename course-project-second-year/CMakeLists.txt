cmake_minimum_required(VERSION 3.21.2)
project(course_project CXX)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(course_project main.cpp types/complex.cpp types/parser.cpp)
