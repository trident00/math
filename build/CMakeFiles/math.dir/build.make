# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\math

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\math\build

# Include any dependencies generated for this target.
include CMakeFiles/math.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/math.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/math.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/math.dir/flags.make

CMakeFiles/math.dir/src/ast.cpp.obj: CMakeFiles/math.dir/flags.make
CMakeFiles/math.dir/src/ast.cpp.obj: CMakeFiles/math.dir/includes_CXX.rsp
CMakeFiles/math.dir/src/ast.cpp.obj: C:/math/src/ast.cpp
CMakeFiles/math.dir/src/ast.cpp.obj: CMakeFiles/math.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\math\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/math.dir/src/ast.cpp.obj"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/math.dir/src/ast.cpp.obj -MF CMakeFiles\math.dir\src\ast.cpp.obj.d -o CMakeFiles\math.dir\src\ast.cpp.obj -c C:\math\src\ast.cpp

CMakeFiles/math.dir/src/ast.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/math.dir/src/ast.cpp.i"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\math\src\ast.cpp > CMakeFiles\math.dir\src\ast.cpp.i

CMakeFiles/math.dir/src/ast.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/math.dir/src/ast.cpp.s"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\math\src\ast.cpp -o CMakeFiles\math.dir\src\ast.cpp.s

CMakeFiles/math.dir/src/parser.cpp.obj: CMakeFiles/math.dir/flags.make
CMakeFiles/math.dir/src/parser.cpp.obj: CMakeFiles/math.dir/includes_CXX.rsp
CMakeFiles/math.dir/src/parser.cpp.obj: C:/math/src/parser.cpp
CMakeFiles/math.dir/src/parser.cpp.obj: CMakeFiles/math.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\math\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/math.dir/src/parser.cpp.obj"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/math.dir/src/parser.cpp.obj -MF CMakeFiles\math.dir\src\parser.cpp.obj.d -o CMakeFiles\math.dir\src\parser.cpp.obj -c C:\math\src\parser.cpp

CMakeFiles/math.dir/src/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/math.dir/src/parser.cpp.i"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\math\src\parser.cpp > CMakeFiles\math.dir\src\parser.cpp.i

CMakeFiles/math.dir/src/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/math.dir/src/parser.cpp.s"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\math\src\parser.cpp -o CMakeFiles\math.dir\src\parser.cpp.s

CMakeFiles/math.dir/src/main2.cpp.obj: CMakeFiles/math.dir/flags.make
CMakeFiles/math.dir/src/main2.cpp.obj: CMakeFiles/math.dir/includes_CXX.rsp
CMakeFiles/math.dir/src/main2.cpp.obj: C:/math/src/main2.cpp
CMakeFiles/math.dir/src/main2.cpp.obj: CMakeFiles/math.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:\math\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/math.dir/src/main2.cpp.obj"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/math.dir/src/main2.cpp.obj -MF CMakeFiles\math.dir\src\main2.cpp.obj.d -o CMakeFiles\math.dir\src\main2.cpp.obj -c C:\math\src\main2.cpp

CMakeFiles/math.dir/src/main2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/math.dir/src/main2.cpp.i"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\math\src\main2.cpp > CMakeFiles\math.dir\src\main2.cpp.i

CMakeFiles/math.dir/src/main2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/math.dir/src/main2.cpp.s"
	C:\ProgramData\mingw64\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\math\src\main2.cpp -o CMakeFiles\math.dir\src\main2.cpp.s

# Object files for target math
math_OBJECTS = \
"CMakeFiles/math.dir/src/ast.cpp.obj" \
"CMakeFiles/math.dir/src/parser.cpp.obj" \
"CMakeFiles/math.dir/src/main2.cpp.obj"

# External object files for target math
math_EXTERNAL_OBJECTS =

math.exe: CMakeFiles/math.dir/src/ast.cpp.obj
math.exe: CMakeFiles/math.dir/src/parser.cpp.obj
math.exe: CMakeFiles/math.dir/src/main2.cpp.obj
math.exe: CMakeFiles/math.dir/build.make
math.exe: CMakeFiles/math.dir/linkLibs.rsp
math.exe: CMakeFiles/math.dir/objects1.rsp
math.exe: CMakeFiles/math.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=C:\math\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable math.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\math.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/math.dir/build: math.exe
.PHONY : CMakeFiles/math.dir/build

CMakeFiles/math.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\math.dir\cmake_clean.cmake
.PHONY : CMakeFiles/math.dir/clean

CMakeFiles/math.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\math C:\math C:\math\build C:\math\build C:\math\build\CMakeFiles\math.dir\DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/math.dir/depend

