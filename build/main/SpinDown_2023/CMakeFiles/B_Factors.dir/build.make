# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 4.0

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build"

# Include any dependencies generated for this target.
include main/SpinDown_2023/CMakeFiles/B_Factors.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include main/SpinDown_2023/CMakeFiles/B_Factors.dir/compiler_depend.make

# Include the progress variables for this target.
include main/SpinDown_2023/CMakeFiles/B_Factors.dir/progress.make

# Include the compile flags for this target's objects.
include main/SpinDown_2023/CMakeFiles/B_Factors.dir/flags.make

main/SpinDown_2023/CMakeFiles/B_Factors.dir/codegen:
.PHONY : main/SpinDown_2023/CMakeFiles/B_Factors.dir/codegen

main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o: main/SpinDown_2023/CMakeFiles/B_Factors.dir/flags.make
main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/B_Factors.cpp
main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o: main/SpinDown_2023/CMakeFiles/B_Factors.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o -MF CMakeFiles/B_Factors.dir/B_Factors.cpp.o.d -o CMakeFiles/B_Factors.dir/B_Factors.cpp.o -c "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/B_Factors.cpp"

main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/B_Factors.dir/B_Factors.cpp.i"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/B_Factors.cpp" > CMakeFiles/B_Factors.dir/B_Factors.cpp.i

main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/B_Factors.dir/B_Factors.cpp.s"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/B_Factors.cpp" -o CMakeFiles/B_Factors.dir/B_Factors.cpp.s

# Object files for target B_Factors
B_Factors_OBJECTS = \
"CMakeFiles/B_Factors.dir/B_Factors.cpp.o"

# External object files for target B_Factors
B_Factors_EXTERNAL_OBJECTS =

main/SpinDown_2023/B_Factors: main/SpinDown_2023/CMakeFiles/B_Factors.dir/B_Factors.cpp.o
main/SpinDown_2023/B_Factors: main/SpinDown_2023/CMakeFiles/B_Factors.dir/build.make
main/SpinDown_2023/B_Factors: libCompactStar.a
main/SpinDown_2023/B_Factors: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/dependencies/lib/Zaki/Darwin/arm64/libZaki.a
main/SpinDown_2023/B_Factors: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/dependencies/lib/Confind/Darwin/arm64/libConfind.a
main/SpinDown_2023/B_Factors: /opt/local/lib/libgsl.dylib
main/SpinDown_2023/B_Factors: /opt/local/lib/libgslcblas.dylib
main/SpinDown_2023/B_Factors: /Library/Frameworks/Python.framework/Versions/3.12/lib/libpython3.12.dylib
main/SpinDown_2023/B_Factors: main/SpinDown_2023/CMakeFiles/B_Factors.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable B_Factors"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/B_Factors.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/SpinDown_2023/CMakeFiles/B_Factors.dir/build: main/SpinDown_2023/B_Factors
.PHONY : main/SpinDown_2023/CMakeFiles/B_Factors.dir/build

main/SpinDown_2023/CMakeFiles/B_Factors.dir/clean:
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && $(CMAKE_COMMAND) -P CMakeFiles/B_Factors.dir/cmake_clean.cmake
.PHONY : main/SpinDown_2023/CMakeFiles/B_Factors.dir/clean

main/SpinDown_2023/CMakeFiles/B_Factors.dir/depend:
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023/CMakeFiles/B_Factors.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : main/SpinDown_2023/CMakeFiles/B_Factors.dir/depend

