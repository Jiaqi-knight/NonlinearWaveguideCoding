cmake_minimum_required(VERSION 3.6)
project(cpp)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        interp.cpp
        interp.h
        main.cpp
        makefile)

add_executable(cpp ${SOURCE_FILES})