cmake_minimum_required(VERSION 3.0)
project(SHOpt)

set( CMAKE_POSITION_INDEPENDENT_CODE ON )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( SRC_DIRECTORY ${CMAKE_BINARY_DIR}/src )

if (EXISTS $ENV{CMAKE_INCLUDE_PATH})
    include_directories($ENV{CMAKE_INCLUDE_PATH})
endif()

include_directories(/Users/jayyao/scuff-em-installation/include/scuff-em)
link_directories(/Users/jayyao/scuff-em-installation/lib)

# Boost and its components
find_package(Boost COMPONENTS system thread)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_library( basisfunctions SHARED ${SRC_DIRECTORY}/BaseGeo.cpp ${SRC_DIRECTORY}/SphericalHarmonics.cpp ${SRC_DIRECTORY}/genRandomSpherePoints.cpp )
else()
    message(STATUS "The basisfunctions library requires Boost, and will not be compiled.")
    message(STATUS "The SHMesher libraries require the basisfunctions library, and will not be compiled.")
    return()
endif()

# Find Eigen3 (requires 3.1.0 or greater)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if ( NOT EIGEN3_FOUND)
    message(STATUS "The main.cpp component requires Eigen3 >=3.1.0, and will not be compiled.")
    return()  
endif()
include( ${EIGEN3_USE_FILE} )

# Find NLOPT
find_package(NLopt)
if ( NOT NLopt_FOUND)
    message(STATUS "The main.cpp component requires NLopt, and will not be compiled.")
    return()
endif()

# Creating entries for target: SHOpt
add_executable( SHOpt  ${SRC_DIRECTORY}/main.cpp ${SRC_DIRECTORY}/scuff-scatter.cpp ${SRC_DIRECTORY}/LDOSmath.cpp)
target_link_libraries( SHOpt scuff basisfunctions ${BOOST_LIBRARIES} Eigen3::Eigen ${NLOPT_LIBRARIES})
