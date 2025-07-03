cmake_minimum_required(VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(roptlite LANGUAGES CXX)

# Build all binaries in a separate directory
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#add extra search paths for libraries and includes
SET(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "The directory the headers are installed in")
SET(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE STRING "Directory where lib will install")
SET(DOC_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/share/doc/${PROJECT_NAME}" CACHE PATH "Path to the documentation")


# add_compile_definitions(__cplusplus=201703L)
# message(STATUS "C++ Standard: ${__cplusplus}")
set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
add_compile_options(-march=native) # Enable faster instruction sets

# include(GNUInstallDirs)

# include(GNUInstallDirs)
# set(ROPTLITE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/roptlite)

################################################################################
# Export package for use from the build tree
# EXPORT( PACKAGE ${PROJECT_NAME} )

# # Directory for built libraries
# set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# # Directory for built executables
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")





# include(version)

add_subdirectory(src)
# add_subdirectory(tests)

# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

install(
    TARGETS ${PROJECT_NAME}
    EXPORT roptliteTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/roptlite)

export(
    EXPORT roptliteTargets
    FILE roptliteTargets.cmake
    NAMESPACE roptlite::)
include(CMakePackageConfigHelpers)


# configure_package_config_file(
#     cmake/roptliteConfig.cmake "${CMAKE_CURRENT_BINARY_DIR}/roptliteConfig.cmake"
#     INSTALL_DESTINATION "${LIB_INSTALL_DIR}/cmake/")


#   install(
#   FILES ${CMAKE_CURRENT_BINARY_DIR}/roptliteConfig.cmake
#         ${CMAKE_CURRENT_BINARY_DIR}/roptliteTargets.cmake
#   DESTINATION share/roptlite/cmake)

# install header files
# install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
# 		DESTINATION ${ROPTLITE_INSTALL_INCLUDEDIR}
# 		FILES_MATCHING PATTERN "*.h"
# 		)
