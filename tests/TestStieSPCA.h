/*
This is the test file for the Sparse PCA problem defined in StieSPCA.h and StieSPCA.cpp.

---- WH
*/

#ifndef TESTSTIESPCA_H
#define TESTSTIESPCA_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "Others/randgen.h"
//#include "Others/MyMatrix.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test/DriverMexProb.h"

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/StieSPCA.h"
//#include "Problems/SphereTxRQ/SphereTxRQ.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
//#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel.h"
//#include "Manifolds/SphereTx/SphereTx.h"

/*Linesearch based solvers*/
#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/IRPG.h"
#include "Solvers/IARPG.h"

/*Trust-region based solvers*/
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

/*The global head file*/
#include "Others/def.h"

using namespace ROPTLITE;

/*The main test function*/
void testStieSPCA(void);
Vector &invHV(const Vector &B, realdp alpha, const Vector &V, Vector *result);
Vector SOCforSPCA(Vector A, realdp mu, Vector initX);

#endif // end of TESTSTIESPCA_H
