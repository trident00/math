@echo off
REM Create the build directory if it does not exist
if not exist "build" (
    mkdir build
)

REM Change to the build directory
cd build

REM Run CMake to configure the project
cmake -G "MinGW Makefiles" -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

REM Build the project using MinGW Makefiles
mingw32-make
