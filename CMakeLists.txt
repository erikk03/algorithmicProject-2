# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(TriangulationProject)

# CGAL and its components
find_package(CGAL QUIET COMPONENTS Qt5)

if (NOT CGAL_FOUND)
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()

# Boost and its components
find_package(Boost 1.74 REQUIRED COMPONENTS system filesystem thread)

if (NOT Boost_FOUND)
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

# Include directories for local headers
include_directories(${CGAL_INCLUDE_DIRS} includes)  # Add 'includes' directory

# Create the main executable
add_executable(main main.cpp)  # main.cpp is your main source file

# Link the executable to CGAL and Boost libraries
target_link_libraries(main PRIVATE CGAL::CGAL ${Boost_LIBRARIES})

# Check if the Qt5 components for CGAL are found
if (CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER)  # Enable the Basic Viewer
    target_link_libraries(main PRIVATE CGAL::CGAL_Qt5)
endif()

# Set the build type to Release for better performance
set(CMAKE_BUILD_TYPE Release)

