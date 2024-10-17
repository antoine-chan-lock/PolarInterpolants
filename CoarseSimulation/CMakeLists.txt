cmake_minimum_required(VERSION 3.16)

if (APPLE)
  set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
  set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
endif()

project(example)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp -O2")

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(CHOLMOD)
find_package(SPQR)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../CommonLibraries/Eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libsim)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/symStrains)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../CommonLibraries/libigl/include)

include_directories(${CHOLMOD_INCLUDE_DIR})
include_directories(${SPQR_INCLUDE_DIR})

file(GLOB SRC_FILES *.cpp
                    libsim/*
                    symStrains/*.inc)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME}  ${CHOLMOD_LIBRARIES})
target_link_libraries(${PROJECT_NAME}  ${SPQR_LIBRARIES})