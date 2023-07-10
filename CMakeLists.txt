cmake_minimum_required(VERSION 3.17)
project(virus_research)

set(CMAKE_CXX_STANDARD 14)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

find_package (Eigen3 3.3 REQUIRED)

include_directories(eigen-3.4.0)

add_executable(virus_research virus.cpp Matrix6fFileReader.cpp IcosahedralGroup.cpp TetrahedralGroup.cpp TetrahedralGroup.hpp Matrix6fGroup.cpp Matrix6fGroup.hpp PermutationGroup.cpp PermutationGroup.hpp std_vector_functions.cpp std_vector_functions.hpp MatrixFunctions.cpp MatrixFunctions.hpp DihedralGroupOnTen.cpp DihedralGroupOnTen.hpp DihedralGroupOnSix.cpp DihedralGroupOnSix.hpp GeneratingVectorsForViruses.cpp GeneratingVectorsForViruses.hpp outputResults.cpp outputResults.hpp)
add_executable(print_output_results printOutputResults.cpp outputResults.cpp Matrix6fFileReader.cpp MatrixFunctions.cpp outputResults.hpp)
add_executable(entry_sampling sampling.cpp GeneratingVectorsForViruses.cpp IcosahedralGroup.cpp Matrix6fFileReader.cpp Matrix6fGroup.cpp MatrixFunctions.cpp std_vector_functions.cpp)
target_link_libraries (virus_research Eigen3::Eigen)
