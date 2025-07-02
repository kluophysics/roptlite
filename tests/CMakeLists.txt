set(SRC_LIST

test/TestCFRankQ2FBlindDecon2D.cpp
test/TestCSFRQPhaseRetrieval.cpp
test/TestCStieBrockett.cpp
test/TestElement.cpp
test/TestEucQuadratic.cpp
test/TestFRankE3FMatCompletion.cpp
test/TestFRankESparseApprox.cpp
test/TestFRankETextureInpainting.cpp
test/TestFRankEWeightApprox.cpp
test/TestFRankQ2FMatCompletion.cpp
test/TestGrassMatCompletion.cpp
test/TestGrassPCA.cpp
test/TestGrassRQ.cpp
test/TestGrassSVPCA.cpp
test/TestObliqueSPCA.cpp
test/TestPoincareEmbeddings.cpp
test/TestProdStieSumBrockett.cpp
test/TestSFRQLyapunov.cpp
test/TestSPDKarcherMean.cpp
test/TestSphereSparsestVector.cpp
test/TestStieBrockett.cpp
test/TestStieSPCA.cpp
test/TestStieSoftICA.cpp
)

target_include_directories(roptlite
    PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/Manifolds
      ${CMAKE_CURRENT_SOURCE_DIR}/Others
      ${CMAKE_CURRENT_SOURCE_DIR}/Others/SparseBLAS
      ${CMAKE_CURRENT_SOURCE_DIR}/Others/fftw
      ${CMAKE_CURRENT_SOURCE_DIR}/Others/wavelet
      ${CMAKE_CURRENT_SOURCE_DIR}/Problems
      ${CMAKE_CURRENT_SOURCE_DIR}/Solvers
      ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper
      ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/blas
      ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/lapack
      # ${CMAKE_CURRENT_SOURCE_DIR}/test
)



# if(ENABLE_COVERAGE)
#   add_coverage(surchem)
# endif()

# IF (BUILD_TESTING)
#   if(ENABLE_MPI)
#     add_subdirectory(test)
#   endif()
# endif()
