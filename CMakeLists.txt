cmake_minimum_required(VERSION 3.4)

project(line-intersection)

configure_file(test.txt test.txt COPYONLY)
configure_file(test1.txt test1.txt COPYONLY)
configure_file(test2.txt test2.txt COPYONLY)
configure_file(test3.txt test3.txt COPYONLY)
configure_file(test4.txt test4.txt COPYONLY)

add_library(fixed_real "fixedreal.h" "fixedreal.cpp")

add_definitions(--std=c++14 -g -DLOG_LEVEL=DEBUG)

add_executable(${PROJECT_NAME} "main.cpp" "logger.h" "ts_printf.h")
target_link_libraries(${PROJECT_NAME} fixed_real)
target_link_libraries(${PROJECT_NAME} gmpxx gmp)
