/*
This file defines the class of the Riemannian trust-region Newton method in [ABG2007]
	[ABG2007]: P.-A. Absil, C. G. Baker, and K. A. Gallivan. Trust-region methods on Riemannian manifolds. 
		Foundations of Computational Mathematics, 7(3):303?30, 2007.

Solvers --> SolversSM --> SolversSMTR --> RTRNewton

---- WH
*/

#ifndef RTRNEWTON_H
#define RTRNEWTON_H

#include "Solvers/SolversSMTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLITE{

	class RTRNewton : public SolversSMTR{
	public:
		/*The contructor of RTRNewton method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		RTRNewton(const Problem *prob, const Variable *initialx);

	protected:
        /*Call Solvers::SetProbX function; initialize temporary vectors; and indicate RTRNewton needs action of Hessian.
        INPUT:    prob is the problem which defines the cost function, gradient and possible the action of Hessian
        and specifies the manifold of domain.
        initialx is the initial iterate.*/
        virtual void SetProbX(const Problem *prob, const Variable *initialx);

        /*Setting parameters (member variables) to be default values */
        virtual void SetDefaultParams();
        
		/*Compute result = H[Eta], where H is the Hessian*/
		virtual Vector &HessianEta(const Vector &Eta, Vector *result);
	};
}; /*end of ROPTLITE namespace*/

#endif /* end of RTRNEWTON_H */
