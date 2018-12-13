cmake_minimum_required(VERSION 3.0)
project( meshutils )

set( CMAKE_POSITION_INDEPENDENT_CODE ON )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( SRC_DIRECTORY ${CMAKE_BINARY_DIR}/src )

if (EXISTS $ENV{CMAKE_INCLUDE_PATH})
    include_directories($ENV{CMAKE_INCLUDE_PATH})
endif()

# meshutils library
add_library( meshutils SHARED ${SRC_DIRECTORY}/Mesh.cpp ${SRC_DIRECTORY}/MeshUtils.cpp )

# Boost and its components
find_package(Boost COMPONENTS system thread)

if (Boost_FOUND)
    add_library( basisfunctions SHARED ${SRC_DIRECTORY}/BaseGeo.cpp ${SRC_DIRECTORY}/SphericalHarmonics.cpp ${SRC_DIRECTORY}/genRandomSpherePoints.cpp )
else()
    message(STATUS "The basisfunctions library requires Boost, and will not be compiled.")
    message(STATUS "The SHMesher libraries require the basisfunctions library, and will not be compiled.")
    return()
endif()

# CGAL and its components
find_package(CGAL)

if( NOT CGAL_FOUND )
    message(STATUS "The cgal_mesher component requires the CGAL library, and will not be compiled.")
    return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

MESSAGE(STATUS "CGAL_LIBRARIES: " ${CGAL_LIBRARIES})
MESSAGE(STATUS "CGAL_3RD_LIBRARIES: " ${CGAL_3RD_PARTY_LIBRARIES})
# Creating entries for target: SHMesher
add_executable( SHMesher  ${SRC_DIRECTORY}/runMesher.cpp ${SRC_DIRECTORY}/cgalDistMesh.cpp )
target_link_libraries( SHMesher meshutils basisfunctions ${BOOST_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
