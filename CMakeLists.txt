cmake_minimum_required(VERSION "3.16.3")
set (CMAKE_CXX_STANDARD 20)

CMAKE_POLICY(SET CMP0081 NEW)
CMAKE_POLICY(SET CMP0015 NEW)
CMAKE_POLICY(SET CMP0048 NEW)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (s-rlbwt-mds)

include_directories(/)
include_directories(internal/)

include_directories(sdsl/include) #SDSL headers are here
link_directories(sdsl/lib) #SDSL lib are here

find_package(OpenMP)

# ips4o
include_directories(external/move-datastructure/extern/ips4o)

# concurrentqueue
include_directories(external/move-datastructure/extern/concurrentqueue)

# malloc_count
include_directories(external/move-datastructure/extern/malloc_count)

# libsais
include_directories(external/move-datastructure/extern/libsais/src)

# move datastructure
include_directories(
    external/move-datastructure/src/ external/move-datastructure/src/mdsb/ external/move-datastructure/src/misc/
    external/move-datastructure/include/ external/move-datastructure/include/mdsb/ external/move-datastructure/include/misc/
)
add_subdirectory(external/move-datastructure)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -fstrict-aliasing -ftree-vectorize -funroll-loops -finline-functions -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-ggdb3 -Wall -Wextra -march=native -fstrict-aliasing")

add_executable(ri-mds-build ri-mds-build.cpp)
target_link_libraries(ri-mds-build dl OpenMP::OpenMP_CXX sdsl malloc_count libsais mds)

add_executable(ri-mds-locate ri-mds-locate.cpp)
target_link_libraries(ri-mds-locate dl OpenMP::OpenMP_CXX sdsl malloc_count libsais mds)

add_executable(ri-mds-count ri-mds-count.cpp)
target_link_libraries(ri-mds-count dl OpenMP::OpenMP_CXX sdsl malloc_count libsais mds)

add_executable(ri-mds-revert ri-mds-revert.cpp)
target_link_libraries(ri-mds-revert dl OpenMP::OpenMP_CXX sdsl malloc_count libsais mds)