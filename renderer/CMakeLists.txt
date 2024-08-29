# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name
project(renderer)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable target with the source files
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/model.cpp
    src/tgaimage.cpp
)

# Include directories
target_include_directories(${PROJECT_NAME} PUBLIC
    ${PROJECT_SOURCE_DIR}/src   # Include the src directory
)
