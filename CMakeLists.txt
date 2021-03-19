cmake_minimum_required(VERSION 3.17)
project(MCMC)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -fopenmp")

# Flags of sources, headers and libraries
set_directory_properties(PROPERTIES COMPILE_DEFINITIONS_DEBUG "_DEBUG")
include_directories(/opt/OpenBLAS/include)
link_directories(/opt/OpenBLAS/lib)
set(OpenGlLinkers -lSDL2 -lGL -lGLEW)
set(SOURCE_FILES main.cpp src/Matrix.h src/Matrix.cpp src/Helper.h src/Helper.cpp src/LAPACK.h src/BLAS.h src/Vector.h src/Vector.cpp src/CIsing.cpp src/CIsing.h src/visual/Mesh.cpp src/visual/Mesh.h src/visual/Display.cpp src/visual/Display.h src/visual/Visual.cpp src/visual/Visual.h)

# executable and linking
add_executable(MCMC ${SOURCE_FILES})
target_link_libraries(MCMC libopenblas.so ${OpenGlLinkers})
