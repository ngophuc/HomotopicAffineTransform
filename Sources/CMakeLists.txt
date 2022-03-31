# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( HomotopicAffineTransform )

cmake_minimum_required(VERSION 3.5)
set (CMAKE_CXX_STANDARD 11)

aux_source_directory(. SRC_LIST)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

#Required in Gmp
#FIND_PACKAGE(GMP REQUIRED)
#INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIRS})
#LINK_DIRECTORIES(${GMP_LIBRARY_DIRS})

set(GMP_DIR "/usr/local/Cellar/gmp/6.2.1")
include_directories(${PROJECT_NAME} ${GMP_DIR}/include)
link_directories(${PROJECT_NAME} ${GMP_DIR}/lib)
find_library(GMP_LIBRARY NAMES gmp libgmp PATHS "${GMP_DIR}")

add_definitions(-g)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} gmp gmpxx)
