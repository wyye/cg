cmake_minimum_required(VERSION 2.8)

project(line-intersection)

configure_file(test.txt test.txt COPYONLY)

set(CXX_STANDARD 17)

add_executable(${PROJECT_NAME} "main.cpp")
