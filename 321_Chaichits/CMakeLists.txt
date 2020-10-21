cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(rt)

add_executable(rt main.cpp Bitmap.cpp)

find_package(OpenMP)

if (OpenMP_CXX_FOUND)
	target_link_libraries(rt OpenMP::OpenMP_CXX)
endif()
