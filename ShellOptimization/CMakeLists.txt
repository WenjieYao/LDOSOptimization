cmake_minimum_required(VERSION 3.0)
project(Shell)

set( CMAKE_POSITION_INDEPENDENT_CODE ON )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( SRC_DIRECTORY ${CMAKE_BINARY_DIR}/src )

if (EXISTS $ENV{CMAKE_INCLUDE_PATH})
    include_directories($ENV{CMAKE_INCLUDE_PATH})
endif()

include_directories(/Users/jayyao/scuff-em-installation/include/scuff-em)
link_directories(/Users/jayyao/scuff-em-installation/lib)


# Find NLOPT
find_package(NLopt)
if ( NOT NLopt_FOUND)
    message(STATUS "The main.cpp component requires NLopt, and will not be compiled.")
    return()
endif()

# Creating entries for target: SHOpt
add_executable( Shell  ${SRC_DIRECTORY}/main.cpp ${SRC_DIRECTORY}/scuff-scatter.cpp ${SRC_DIRECTORY}/LDOSmath.cpp)
target_link_libraries( Shell scuff ${NLOPT_LIBRARIES})
