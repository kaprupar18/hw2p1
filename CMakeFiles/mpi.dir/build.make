# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jariy/hw2p1t/hw2p1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jariy/hw2p1t/hw2p1

# Include any dependencies generated for this target.
include CMakeFiles/mpi.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mpi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mpi.dir/flags.make

CMakeFiles/mpi.dir/mpi.cpp.o: CMakeFiles/mpi.dir/flags.make
CMakeFiles/mpi.dir/mpi.cpp.o: mpi.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jariy/hw2p1t/hw2p1/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mpi.dir/mpi.cpp.o"
	/usr/lib64/ccache/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mpi.dir/mpi.cpp.o -c /home/jariy/hw2p1t/hw2p1/mpi.cpp

CMakeFiles/mpi.dir/mpi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mpi.dir/mpi.cpp.i"
	/usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/jariy/hw2p1t/hw2p1/mpi.cpp > CMakeFiles/mpi.dir/mpi.cpp.i

CMakeFiles/mpi.dir/mpi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mpi.dir/mpi.cpp.s"
	/usr/lib64/ccache/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/jariy/hw2p1t/hw2p1/mpi.cpp -o CMakeFiles/mpi.dir/mpi.cpp.s

CMakeFiles/mpi.dir/mpi.cpp.o.requires:
.PHONY : CMakeFiles/mpi.dir/mpi.cpp.o.requires

CMakeFiles/mpi.dir/mpi.cpp.o.provides: CMakeFiles/mpi.dir/mpi.cpp.o.requires
	$(MAKE) -f CMakeFiles/mpi.dir/build.make CMakeFiles/mpi.dir/mpi.cpp.o.provides.build
.PHONY : CMakeFiles/mpi.dir/mpi.cpp.o.provides

CMakeFiles/mpi.dir/mpi.cpp.o.provides.build: CMakeFiles/mpi.dir/mpi.cpp.o

# Object files for target mpi
mpi_OBJECTS = \
"CMakeFiles/mpi.dir/mpi.cpp.o"

# External object files for target mpi
mpi_EXTERNAL_OBJECTS =

mpi: CMakeFiles/mpi.dir/mpi.cpp.o
mpi: CMakeFiles/mpi.dir/build.make
mpi: /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib/libmpifort.so
mpi: /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib/release_mt/libmpi.so
mpi: /opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib/libmpigi.a
mpi: /usr/lib64/libdl.so
mpi: /usr/lib64/librt.so
mpi: /usr/lib64/libpthread.so
mpi: CMakeFiles/mpi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable mpi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mpi.dir/build: mpi
.PHONY : CMakeFiles/mpi.dir/build

CMakeFiles/mpi.dir/requires: CMakeFiles/mpi.dir/mpi.cpp.o.requires
.PHONY : CMakeFiles/mpi.dir/requires

CMakeFiles/mpi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mpi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mpi.dir/clean

CMakeFiles/mpi.dir/depend:
	cd /home/jariy/hw2p1t/hw2p1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jariy/hw2p1t/hw2p1 /home/jariy/hw2p1t/hw2p1 /home/jariy/hw2p1t/hw2p1 /home/jariy/hw2p1t/hw2p1 /home/jariy/hw2p1t/hw2p1/CMakeFiles/mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mpi.dir/depend

