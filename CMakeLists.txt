cmake_minimum_required (VERSION 2.8.12)
project (CPlantBox)


# set default build type to release Release
set(CMAKE_BUILD_TYPE Release)

# add source directory to the include path
include_directories(${PROJECT_SOURCE_DIR}/src)

# add subdirectories
add_subdirectory(src)
add_subdirectory(tutorial)
add_subdirectory(debug_test)
