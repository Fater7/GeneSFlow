cmake_minimum_required(VERSION 3.1.0) 
project(Demo)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# 目录源文件
aux_source_directory(. DIR_SRCS)
# 头文件目录
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/GeneSFlow)
# 子文件夹
add_subdirectory(GeneSFlow)
# Target
add_executable(Demo ${DIR_SRCS})
target_link_libraries(Demo GeneSFlow)