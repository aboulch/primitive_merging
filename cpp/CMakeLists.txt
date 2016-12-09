cmake_minimum_required(VERSION 2.8)
project(ComparePrim)


#use openMP
find_package( OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fPIC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#Boost
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lboost_program_options")

add_executable(TestPrim main.cpp 
surface_comparator_MW.h 
surface_comparator_KS.h)
target_link_libraries(TestPrim ${RECONS_LIBS} ${Boost_LIBRARIES})
