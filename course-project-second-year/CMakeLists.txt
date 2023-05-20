cmake_minimum_required(VERSION 3.21.2)
project(course_project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(course_project_lib types/complex.cpp utils/parser.cpp)

add_executable(course_project main.cpp)
target_link_libraries(course_project course_project_lib)

add_executable(course_project_test tests/test.cpp)
target_link_libraries(course_project_test course_project_lib)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

target_link_libraries(course_project_test gtest_main)
add_test(NAME test_all COMMAND course_project_test)
