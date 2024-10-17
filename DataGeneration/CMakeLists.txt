cmake_minimum_required(VERSION 3.16)

project(example LANGUAGES C CXX)

if (APPLE)
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
endif()

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp -O2")

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../CommonLibraries/Eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/libsim
    ${CMAKE_CURRENT_SOURCE_DIR}/libsim/Objects
    ${CMAKE_CURRENT_SOURCE_DIR}/libsim/Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/libsim/DataGen
    ${CMAKE_CURRENT_SOURCE_DIR}/libsim/Materials
    ${CMAKE_CURRENT_SOURCE_DIR}/../CommonLibraries/libigl/include
)

file(GLOB SRC_FILES 
    *.cpp
    libsim/*
    libsim/Objects/*
    libsim/Elements/*
    libsim/DataGen/*
    libsim/Materials/*
)

add_executable(${PROJECT_NAME} ${SRC_FILES})
