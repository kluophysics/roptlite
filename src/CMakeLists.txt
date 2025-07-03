# roptlite requires lapack, blas, and m libraries
find_library(LAPACK lapack DOC "LAPACK library to link against" REQUIRED)
find_library(BLAS blas DOC "BLAS library to link against" REQUIRED)
find_library(M m "M library to link against" REQUIRED)


set(SRC_LIST

Manifolds/CFixedRankQ2F.cpp
Manifolds/CStiefel.cpp
Manifolds/CSymFixedRankQ.cpp
Manifolds/Element.cpp
Manifolds/Euclidean.cpp
Manifolds/FixedRankE.cpp
Manifolds/FixedRankE3F.cpp
Manifolds/FixedRankQ2F.cpp
Manifolds/Grassmann.cpp
Manifolds/Manifold.cpp
Manifolds/MultiManifolds.cpp
Manifolds/PoincareBall.cpp
Manifolds/SPDManifold.cpp
Manifolds/SmartSpace.cpp
Manifolds/Sphere.cpp
Manifolds/SphereTx.cpp
Manifolds/Stiefel.cpp
Manifolds/SymFixedRankQ.cpp

Others/BlasLapackCppWrapper.cpp
Others/ForDebug.cpp
Others/MinPNormConHull.cpp
Others/SparseMatrix.cpp
Others/Spline.cpp
Others/Timer.cpp
Others/randgen.cpp

Others/SparseBLAS/nist_spblas.cpp

Others/wavelet/wavelet.cpp

Problems/CFRankQ2FBlindDecon2D.cpp
Problems/CSFRQPhaseRetrieval.cpp
Problems/CStieBrockett.cpp
Problems/EucQuadratic.cpp
Problems/FRankE3FMatCompletion.cpp
Problems/FRankESparseApprox.cpp
Problems/FRankETextureInpainting.cpp
Problems/FRankEWeightApprox.cpp
Problems/FRankQ2FMatCompletion.cpp
Problems/GrassMatCompletion.cpp
Problems/GrassPCA.cpp
Problems/GrassRQ.cpp
Problems/GrassSVPCA.cpp
Problems/ObliqueSPCA.cpp
Problems/PoincareEmbeddings.cpp
Problems/Problem.cpp
Problems/ProdStieSumBrockett.cpp
Problems/SFRQLyapunov.cpp
Problems/SPDKarcherMean.cpp
Problems/SphereConvexHull.cpp
Problems/SphereSparsestVector.cpp
Problems/SphereTxRQ.cpp
Problems/StieBrockett.cpp
Problems/StieSPCA.cpp
Problems/StieSoftICA.cpp
Problems/juliaProblem.cpp
Problems/mexProblem.cpp

Solvers/IARPG.cpp
Solvers/IRPG.cpp
Solvers/LRBFGS.cpp
Solvers/LRBFGSSub.cpp
Solvers/LRBroydenFamily.cpp
Solvers/LRTRSR1.cpp
Solvers/LRTRSR1woR.cpp
Solvers/RADAM.cpp
Solvers/RADAMSP.cpp
Solvers/RAMSGRAD.cpp
Solvers/RAMSGRADSP.cpp
Solvers/RBFGS.cpp
Solvers/RBFGSSub.cpp
Solvers/RBroydenFamily.cpp
Solvers/RCG.cpp
Solvers/RGS.cpp
Solvers/RNewton.cpp
Solvers/RSD.cpp
Solvers/RSGD.cpp
Solvers/RSVRG.cpp
Solvers/RTRNewton.cpp
Solvers/RTRSD.cpp
Solvers/RTRSR1.cpp
Solvers/RWRBFGS.cpp
Solvers/SVRLRBFGS.cpp
Solvers/SVRLRBroydenFamily.cpp
Solvers/Solvers.cpp
Solvers/SolversNSM.cpp
Solvers/SolversNSMPGLS.cpp
Solvers/SolversNSMSub.cpp
Solvers/SolversNSMSubLS.cpp
Solvers/SolversSM.cpp
Solvers/SolversSMLS.cpp
Solvers/SolversSMSVRG.cpp
Solvers/SolversSMSto.cpp
Solvers/SolversSMTR.cpp
)


add_library(roptlite  ${SRC_LIST})
# set_target_properties(roptlite PROPERTIES OUTPUT_NAME "roptlite")

target_include_directories(roptlite
    PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Manifolds>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Others>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Others/SparseBLAS>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Others/fftw>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Others/wavelet>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Problems>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cwrapper>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/blas>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/lapack>

    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Manifolds>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Others>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Others/SparseBLAS>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Others/fftw>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Others/wavelet>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Problems>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/Solvers>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/cwrapper>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/cwrapper/blas>
    $<INSTALL_INTERFACE:${ROPTLITE_INSTALL_INCLUDEDIR}/cwrapper/lapack>
  
    )

target_link_libraries(
    roptlite
    PUBLIC ${M_LIBRARIES}
    PUBLIC ${BLAS_LIBRARIES}
    PUBLIC ${LAPACK_LIBRARIES}
    )
