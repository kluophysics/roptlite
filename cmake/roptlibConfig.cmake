@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

include("${CMAKE_CURRENT_LIST_DIR}/roptliteTargets.cmake")

find_dependency(M REQUIRED)
find_dependency(BLAS REQUIRED)
find_dependency(LAPACK REQUIRED)

check_required_components(roptlite)