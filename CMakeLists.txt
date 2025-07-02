cmake_minimum_required(VERSION 3.15)
project(roptlite VERSION 0.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set module path for custom CMake scripts
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Compiler flags
add_compile_options(-march=native)

# Find required libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(MATH_LIB m)

# Include directories
include(GNUInstallDirs)
set(ROPTLITE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/roptlite)

# # Directory for built libraries
# set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# # Directory for built executables
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")


# Versioning
include(version)

# Add subdirectories
add_subdirectory(src)
# add_subdirectory(tests) # Uncomment if tests are available

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------
include(CMakePackageConfigHelpers)

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/roptlite)

configure_package_config_file(
    cmake/roptliteConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/roptliteConfig.cmake"
    INSTALL_DESTINATION "${INSTALL_CONFIGDIR}"
)

# # Export targets
# export(
#     TARGETS ${PROJECT_NAME}
#     EXPORT roptliteTargets
#     FILE "${CMAKE_CURRENT_BINARY_DIR}/roptliteTargets.cmake"
#     NAMESPACE roptlite::
# )

# Install targets
install(
    TARGETS ${PROJECT_NAME}
    EXPORT roptliteTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install header files
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${ROPTLITE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

# Install CMake config files
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/roptliteConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/roptliteTargets.cmake"
    DESTINATION ${INSTALL_CONFIGDIR}
)

# Optional: Enable testing
include(CTest)
# if(BUILD_TESTING)
#     add_subdirectory(tests)
# endif()

# Print summary
message(STATUS "roptlite version: ${PROJECT_VERSION}")
message(STATUS "Install include dir: ${ROPTLITE_INSTALL_INCLUDEDIR}")
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
