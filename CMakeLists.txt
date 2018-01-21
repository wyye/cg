cmake_minimum_required(VERSION 2.8)

project(line-intersection)

configure_file(test.txt test.txt COPYONLY)

add_library(fixed_real "fixedreal.h" "fixedreal.cpp")

add_executable(${PROJECT_NAME} "main.cpp")
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)
target_link_libraries(${PROJECT_NAME} fixed_real)
