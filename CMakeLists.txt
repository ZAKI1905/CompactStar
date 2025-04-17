cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(CompactStar VERSION 2.0 DESCRIPTION "Compact Star")

string(TIMESTAMP CompactStar_RELEASE_DATE "%m, %d, %Y")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_compile_options(-Wall -Wextra -pedantic -Werror -g)
# add_compile_options(-pthread -Wall -Wextra -pedantic -ggdb3)
add_compile_options(-pthread)

# For using gnuplot in the future
# find_package(Boost REQUIRED COMPONENTS
#     iostreams system filesystem
# )

# set(root_cmake_path "/opt/local/libexec/root6/share/root/cmake")
# set(CMAKE_PREFIX_PATH ${root_cmake_path})
# find_package(ROOT REQUIRED)

# set(ROOT_LIBS
#     ROOT::Core
#     ROOT::Gpad
#     ROOT::Graf
#     ROOT::Hist
# )

find_package(GSL REQUIRED)

set(DEP_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib)
set(DEP_INC_DIR dependencies/include)

set(ZAKI_LIB ${DEP_LIB_DIR}/Zaki/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libZaki.a)
set(CONFIND_LIB ${DEP_LIB_DIR}/Confind/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}/libConfind.a)

# list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# for matplotlib:
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(CompactStarConfig.h.in ${CMAKE_CURRENT_SOURCE_DIR}/CompactStar/Core/CompactStarConfig.h)

find_package(OpenMP REQUIRED)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

add_subdirectory(CompactStar)

add_library(CompactStar STATIC ${CompactStar_SRC_Files})
set_target_properties(CompactStar PROPERTIES VERSION ${PROJECT_VERSION})

target_include_directories(CompactStar 
                            PUBLIC
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${GSL_INCLUDE_DIRS}
                            ${DEP_INC_DIR}
                            # ${ROOT_INCLUDE_DIRS}
                            OpenMP::OpenMP_CXX
                            )

target_link_libraries(CompactStar PUBLIC 
                        ${ZAKI_LIB} 
                        # ${ROOT_LIBS}
                        ${CONFIND_LIB}
                        GSL::gsl
                        GSL::gslcblas
                        Python3::Python # for matplotlib
                        Python3::Module # for matplotlib
                        Python3::NumPy
                        # Boost::iostreams
                        # Boost::system
                        # Boost::filesystem
                        ) 
                            
install(TARGETS CompactStar DESTINATION lib)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include"
        DESTINATION dependencies)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib" 
        DESTINATION dependencies)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/License" DESTINATION Docs)

# Option for generating the doxygen docs
option(GEN_DOCS "Generate documents using doxygen." OFF)

if(GEN_DOCS)
    execute_process (
        COMMAND
        doxygen ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/html" DESTINATION Docs)

add_subdirectory(main)

# foreach(src ${CompactStar_main_SRC_Files})
#     add_executable(${src} "main/${src}.cpp")  
#     target_link_libraries(${src} PUBLIC CompactStar)
#     target_include_directories(${src} PUBLIC  
#                                 ${CMAKE_CURRENT_SOURCE_DIR}/include
#                                 ${GSL_INCLUDE_DIRS}
#                                 ${DEP_INC_DIR}
#                                 )
#     install(FILES "main/${src}.cpp" DESTINATION main)
#     install(TARGETS ${src} DESTINATION main/${src})
# endforeach()

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${CompactStar_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${CompactStar_VERSION_MINOR}")
include(CPack)