cmake_minimum_required(VERSION 2.8)
project(stella)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(cxx11)
include(cxx_config)

enable_language(C CXX Fortran)

include_directories(include)
install(DIRECTORY include/ DESTINATION include)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)


FIND_PACKAGE(Boost 1.44 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIR})

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})
include_directories(${MPI_Fortran_INCLUDE_PATH})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(cedar-deps ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} -lgfortran)

add_subdirectory(src)
add_library(stella ${stella_src})


set(cedar-dir /home/andrew/src/cedar)
set(cedar_LIBRARIES -lcedar -L${cedar-dir}/build)
set(cedar_INCLUDE_DIR ${cedar-dir}/include)

include_directories(${cedar_INCLUDE_DIR})


add_subdirectory(app)
install(TARGETS stella LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
