cmake_minimum_required(VERSION 3.5)


## Project

project(rtiow)

## Project targets

add_executable(rtiow main.cpp)
set_property(TARGET rtiow PROPERTY CXX_STANDARD 17)

