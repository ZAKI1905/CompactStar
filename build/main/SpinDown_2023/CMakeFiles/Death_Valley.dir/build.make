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
include main/SpinDown_2023/CMakeFiles/Death_Valley.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include main/SpinDown_2023/CMakeFiles/Death_Valley.dir/compiler_depend.make

# Include the progress variables for this target.
include main/SpinDown_2023/CMakeFiles/Death_Valley.dir/progress.make

# Include the compile flags for this target's objects.
include main/SpinDown_2023/CMakeFiles/Death_Valley.dir/flags.make

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/codegen:
.PHONY : main/SpinDown_2023/CMakeFiles/Death_Valley.dir/codegen

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o: main/SpinDown_2023/CMakeFiles/Death_Valley.dir/flags.make
main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/Death_Valley.cpp
main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o: main/SpinDown_2023/CMakeFiles/Death_Valley.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o -MF CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o.d -o CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o -c "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/Death_Valley.cpp"

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/Death_Valley.dir/Death_Valley.cpp.i"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/Death_Valley.cpp" > CMakeFiles/Death_Valley.dir/Death_Valley.cpp.i

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/Death_Valley.dir/Death_Valley.cpp.s"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && /opt/homebrew/opt/llvm/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023/Death_Valley.cpp" -o CMakeFiles/Death_Valley.dir/Death_Valley.cpp.s

# Object files for target Death_Valley
Death_Valley_OBJECTS = \
"CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o"

# External object files for target Death_Valley
Death_Valley_EXTERNAL_OBJECTS =

main/SpinDown_2023/Death_Valley: main/SpinDown_2023/CMakeFiles/Death_Valley.dir/Death_Valley.cpp.o
main/SpinDown_2023/Death_Valley: main/SpinDown_2023/CMakeFiles/Death_Valley.dir/build.make
main/SpinDown_2023/Death_Valley: libCompactStar.a
main/SpinDown_2023/Death_Valley: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/dependencies/lib/Zaki/Darwin/arm64/libZaki.a
main/SpinDown_2023/Death_Valley: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My\ Drive/Research/Tools/Coding/CompactStar/dependencies/lib/Confind/Darwin/arm64/libConfind.a
main/SpinDown_2023/Death_Valley: /opt/local/lib/libgsl.dylib
main/SpinDown_2023/Death_Valley: /opt/local/lib/libgslcblas.dylib
main/SpinDown_2023/Death_Valley: /Library/Frameworks/Python.framework/Versions/3.12/lib/libpython3.12.dylib
main/SpinDown_2023/Death_Valley: main/SpinDown_2023/CMakeFiles/Death_Valley.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Death_Valley"
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Death_Valley.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/SpinDown_2023/CMakeFiles/Death_Valley.dir/build: main/SpinDown_2023/Death_Valley
.PHONY : main/SpinDown_2023/CMakeFiles/Death_Valley.dir/build

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/clean:
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" && $(CMAKE_COMMAND) -P CMakeFiles/Death_Valley.dir/cmake_clean.cmake
.PHONY : main/SpinDown_2023/CMakeFiles/Death_Valley.dir/clean

main/SpinDown_2023/CMakeFiles/Death_Valley.dir/depend:
	cd "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/SpinDown_2023" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023" "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/SpinDown_2023/CMakeFiles/Death_Valley.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : main/SpinDown_2023/CMakeFiles/Death_Valley.dir/depend

