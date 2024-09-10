cmake_minimum_required(VERSION 3.16)

project(math)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR}/include)

link_directories(${PROJECT_SOURCE_DIR}/lib)

# add_executable(${PROJECT_NAME} src/main.cpp src/parser.cpp)
add_executable(${PROJECT_NAME} src/ast.cpp src/parser.cpp src/main.cpp)

# Add compiler flags
target_compile_options(${PROJECT_NAME} PRIVATE -Wall)
if (WIN32)
	target_link_libraries(${PROJECT_NAME} PRIVATE -lraylib -lgdi32 -lwinmm)
endif (WIN32)
