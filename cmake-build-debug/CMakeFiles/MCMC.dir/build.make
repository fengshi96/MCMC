# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/shifeng/IDE/clion-2020.2.5/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/shifeng/IDE/clion-2020.2.5/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Barn/Lab/MCMC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Barn/Lab/MCMC/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MCMC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MCMC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MCMC.dir/flags.make

CMakeFiles/MCMC.dir/main.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MCMC.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/main.cpp.o -c /Barn/Lab/MCMC/main.cpp

CMakeFiles/MCMC.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/main.cpp > CMakeFiles/MCMC.dir/main.cpp.i

CMakeFiles/MCMC.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/main.cpp -o CMakeFiles/MCMC.dir/main.cpp.s

CMakeFiles/MCMC.dir/src/Matrix.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/Matrix.cpp.o: ../src/Matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MCMC.dir/src/Matrix.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/Matrix.cpp.o -c /Barn/Lab/MCMC/src/Matrix.cpp

CMakeFiles/MCMC.dir/src/Matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/Matrix.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/Matrix.cpp > CMakeFiles/MCMC.dir/src/Matrix.cpp.i

CMakeFiles/MCMC.dir/src/Matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/Matrix.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/Matrix.cpp -o CMakeFiles/MCMC.dir/src/Matrix.cpp.s

CMakeFiles/MCMC.dir/src/Helper.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/Helper.cpp.o: ../src/Helper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MCMC.dir/src/Helper.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/Helper.cpp.o -c /Barn/Lab/MCMC/src/Helper.cpp

CMakeFiles/MCMC.dir/src/Helper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/Helper.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/Helper.cpp > CMakeFiles/MCMC.dir/src/Helper.cpp.i

CMakeFiles/MCMC.dir/src/Helper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/Helper.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/Helper.cpp -o CMakeFiles/MCMC.dir/src/Helper.cpp.s

CMakeFiles/MCMC.dir/src/Vector.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/Vector.cpp.o: ../src/Vector.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/MCMC.dir/src/Vector.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/Vector.cpp.o -c /Barn/Lab/MCMC/src/Vector.cpp

CMakeFiles/MCMC.dir/src/Vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/Vector.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/Vector.cpp > CMakeFiles/MCMC.dir/src/Vector.cpp.i

CMakeFiles/MCMC.dir/src/Vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/Vector.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/Vector.cpp -o CMakeFiles/MCMC.dir/src/Vector.cpp.s

CMakeFiles/MCMC.dir/src/CIsing.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/CIsing.cpp.o: ../src/CIsing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/MCMC.dir/src/CIsing.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/CIsing.cpp.o -c /Barn/Lab/MCMC/src/CIsing.cpp

CMakeFiles/MCMC.dir/src/CIsing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/CIsing.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/CIsing.cpp > CMakeFiles/MCMC.dir/src/CIsing.cpp.i

CMakeFiles/MCMC.dir/src/CIsing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/CIsing.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/CIsing.cpp -o CMakeFiles/MCMC.dir/src/CIsing.cpp.s

CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o: ../src/visual/Mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o -c /Barn/Lab/MCMC/src/visual/Mesh.cpp

CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/visual/Mesh.cpp > CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.i

CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/visual/Mesh.cpp -o CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.s

CMakeFiles/MCMC.dir/src/visual/Display.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/visual/Display.cpp.o: ../src/visual/Display.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/MCMC.dir/src/visual/Display.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/visual/Display.cpp.o -c /Barn/Lab/MCMC/src/visual/Display.cpp

CMakeFiles/MCMC.dir/src/visual/Display.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/visual/Display.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/visual/Display.cpp > CMakeFiles/MCMC.dir/src/visual/Display.cpp.i

CMakeFiles/MCMC.dir/src/visual/Display.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/visual/Display.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/visual/Display.cpp -o CMakeFiles/MCMC.dir/src/visual/Display.cpp.s

CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o: CMakeFiles/MCMC.dir/flags.make
CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o: ../src/visual/Visual.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o -c /Barn/Lab/MCMC/src/visual/Visual.cpp

CMakeFiles/MCMC.dir/src/visual/Visual.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MCMC.dir/src/visual/Visual.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Barn/Lab/MCMC/src/visual/Visual.cpp > CMakeFiles/MCMC.dir/src/visual/Visual.cpp.i

CMakeFiles/MCMC.dir/src/visual/Visual.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MCMC.dir/src/visual/Visual.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Barn/Lab/MCMC/src/visual/Visual.cpp -o CMakeFiles/MCMC.dir/src/visual/Visual.cpp.s

# Object files for target MCMC
MCMC_OBJECTS = \
"CMakeFiles/MCMC.dir/main.cpp.o" \
"CMakeFiles/MCMC.dir/src/Matrix.cpp.o" \
"CMakeFiles/MCMC.dir/src/Helper.cpp.o" \
"CMakeFiles/MCMC.dir/src/Vector.cpp.o" \
"CMakeFiles/MCMC.dir/src/CIsing.cpp.o" \
"CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o" \
"CMakeFiles/MCMC.dir/src/visual/Display.cpp.o" \
"CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o"

# External object files for target MCMC
MCMC_EXTERNAL_OBJECTS =

MCMC: CMakeFiles/MCMC.dir/main.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/Matrix.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/Helper.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/Vector.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/CIsing.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/visual/Mesh.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/visual/Display.cpp.o
MCMC: CMakeFiles/MCMC.dir/src/visual/Visual.cpp.o
MCMC: CMakeFiles/MCMC.dir/build.make
MCMC: CMakeFiles/MCMC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Barn/Lab/MCMC/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable MCMC"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MCMC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MCMC.dir/build: MCMC

.PHONY : CMakeFiles/MCMC.dir/build

CMakeFiles/MCMC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MCMC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MCMC.dir/clean

CMakeFiles/MCMC.dir/depend:
	cd /Barn/Lab/MCMC/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Barn/Lab/MCMC /Barn/Lab/MCMC /Barn/Lab/MCMC/cmake-build-debug /Barn/Lab/MCMC/cmake-build-debug /Barn/Lab/MCMC/cmake-build-debug/CMakeFiles/MCMC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MCMC.dir/depend

