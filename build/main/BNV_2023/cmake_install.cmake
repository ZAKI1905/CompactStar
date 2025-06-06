# Install script for directory: /Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/opt/homebrew/opt/llvm/bin/llvm-objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/B_Chi_Transition_2023.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/B_Chi_Transition_2023")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Transition_2023" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Transition_2023")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Transition_2023")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/B_Chi_Transition_2023.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/B_Chi_Photon_Decay.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/B_Chi_Photon_Decay")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Photon_Decay" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Photon_Decay")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Photon_Decay")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/B_Chi_Photon_Decay.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/B_Chi_Combo.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/B_Chi_Combo")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Combo" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Combo")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/B_Chi_Combo")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/B_Chi_Combo.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/compose_micro_fig.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/compose_micro_fig")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compose_micro_fig" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compose_micro_fig")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compose_micro_fig")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/compose_micro_fig.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/compare_JMB_comb.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/compare_JMB_comb")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_comb" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_comb")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_comb")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/compare_JMB_comb.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/compare_JMB_indiv.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/compare_JMB_indiv")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_indiv" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_indiv")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/compare_JMB_indiv")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/compare_JMB_indiv.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023" TYPE FILE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/BNV_2023/b_factors.cpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin" TYPE EXECUTABLE FILES "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/b_factors")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/b_factors" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/b_factors")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/opt/local/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/main/BNV_2023/bin/b_factors")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/CMakeFiles/b_factors.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/build/main/BNV_2023/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
