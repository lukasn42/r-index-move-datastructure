cmake_minimum_required(VERSION "3.16.3")
set (CMAKE_CXX_STANDARD 20)

CMAKE_POLICY(SET CMP0081 NEW)
CMAKE_POLICY(SET CMP0015 NEW)
CMAKE_POLICY(SET CMP0048 NEW)

SET(IPS4O_USE_OPENMP ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt-mds)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
include_directories(sdsl/include) #SDSL headers are here

LINK_DIRECTORIES(sdsl/lib) #SDSL lib are here

find_package(OpenMP)

link_libraries(${CMAKE_DL_LIBS})
add_subdirectory(external/move-datastructure)

message("Building in ${CMAKE_BUILD_TYPE} mode")

#set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -fstrict-aliasing -funroll-loops -finline-functions -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -OFast -fstrict-aliasing -march=native")

add_executable(ri-mds-build ri-mds-build.cpp)
TARGET_LINK_LIBRARIES(ri-mds-build sdsl)
TARGET_LINK_LIBRARIES(ri-mds-build divsufsort)
TARGET_LINK_LIBRARIES(ri-mds-build divsufsort64)
TARGET_LINK_LIBRARIES(ri-mds-build OpenMP::OpenMP_CXX mds)

add_executable(ri-mds-locate ri-mds-locate.cpp)
TARGET_LINK_LIBRARIES(ri-mds-locate sdsl)
TARGET_LINK_LIBRARIES(ri-mds-locate divsufsort)
TARGET_LINK_LIBRARIES(ri-mds-locate divsufsort64)
TARGET_LINK_LIBRARIES(ri-mds-locate OpenMP::OpenMP_CXX mds)

add_executable(ri-mds-count ri-mds-count.cpp)
TARGET_LINK_LIBRARIES(ri-mds-count sdsl)
TARGET_LINK_LIBRARIES(ri-mds-count divsufsort)
TARGET_LINK_LIBRARIES(ri-mds-count divsufsort64)
TARGET_LINK_LIBRARIES(ri-mds-count OpenMP::OpenMP_CXX mds)

add_executable(ri-mds-revert ri-mds-revert.cpp)
TARGET_LINK_LIBRARIES(ri-mds-revert sdsl)
TARGET_LINK_LIBRARIES(ri-mds-revert divsufsort)
TARGET_LINK_LIBRARIES(ri-mds-revert divsufsort64)
TARGET_LINK_LIBRARIES(ri-mds-revert OpenMP::OpenMP_CXX mds)

#add_executable(ri-mds-space ri-mds-space.cpp)
#TARGET_LINK_LIBRARIES(ri-mds-space sdsl)
#TARGET_LINK_LIBRARIES(ri-mds-space divsufsort)
#TARGET_LINK_LIBRARIES(ri-mds-space divsufsort64)
