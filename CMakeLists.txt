cmake_minimum_required (VERSION 3.15) 
project (CHNBR
        VERSION 0.1
        DESCRIPTION "Chromatic Numbers library"
        LANGUAGES CXX) 

set(MODERN_CMAKE_BUILD_TESTING ON)

include_directories(include)
add_subdirectory (examples) 
