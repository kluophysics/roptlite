
#include "Problems/FRankQ2FMatCompletion.h"

/*Define the namespace*/
namespace ROPTLITE{

	FRankQ2FMatCompletion::FRankQ2FMatCompletion(unsigned long *inir, unsigned long *injc, unsigned long *injcc, realdp *invals, integer innz, integer inm, integer inn, integer inr)
	{
        ir = inir;
        jc = injc;
        jcc = injcc;
        vals = invals;
        nz = innz;
		m = inm;
		n = inn;
		r = inr;
        
        NumGradHess = false;
	};

	FRankQ2FMatCompletion::~FRankQ2FMatCompletion(void)
	{
	};

    void FRankQ2FMatCompletion::ProjecOmegaGHT(const realdp *G, const realdp *H, integer inm, integer inn, integer inr, unsigned long *inir, unsigned long *injc, unsigned long *injcc, integer nz, realdp *result) const
    {
        for(integer i = 0; i < nz; i++)
        {
            result[i] = 0;
            for(integer j = 0; j < inr; j++)
            {
                /*row:inir[i], col:injc[i]*/
                integer row = inir[i], col = injc[i];
                result[i] += G[row + j * inm] * H[col + j * inn];
            }
        }
    };

	realdp FRankQ2FMatCompletion::f(const Variable &x) const
	{
        Element PGHTmA(nz);
        realdp *PGHTmAptr = PGHTmA.ObtainWriteEntireData();
        Vector G = x.GetElement(0), H = x.GetElement(1);
        ProjecOmegaGHT(G.ObtainReadData(), H.ObtainReadData(), m, n, r, ir, jc, jcc, nz, PGHTmAptr);

        /*P_Omage(G H^T - A)*/
        axpy_(const_cast<integer *> (&nz), &GLOBAL::DNONE, vals, &GLOBAL::IONE, PGHTmAptr, &GLOBAL::IONE);
        
        realdp result = 0.5 * dot_(&nz, PGHTmAptr, &GLOBAL::IONE, PGHTmAptr, &GLOBAL::IONE);
        
        x.AddToFields("PGHTmA", PGHTmA);
        
        return result;
	};
	
	Vector &FRankQ2FMatCompletion::EucGrad(const Variable &x, Vector *result) const
	{
        Vector PGHTmA = x.Field("PGHTmA");
        const realdp *PGHTmAptr = PGHTmA.ObtainReadData();
        SparseMatrix SM(m, n, ir, jc, jcc, const_cast<realdp *> (PGHTmAptr), nz);
        
        Vector G = x.GetElement(0), H = x.GetElement(1);
        result->NewMemoryOnWrite(); /*For multimanifolds, it is important to new memory at the beginning!*/
        result->GetElement(0) = SM * H;
        result->GetElement(1) = (G.GetTranspose() * SM).GetTranspose();
        
        SM.Setnullptr();
        return *result;
	};
	
	Vector &FRankQ2FMatCompletion::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
	{
        Vector PGHTmA = x.Field("PGHTmA");
        const realdp *PGHTmAptr = PGHTmA.ObtainReadData();
        SparseMatrix SM1(m, n, ir, jc, jcc, const_cast<realdp *> (PGHTmAptr), nz);
        
        Vector G = x.GetElement(0), H = x.GetElement(1);
        Vector dG = etax.GetElement(0), dH = etax.GetElement(1);
        const realdp *Gptr = G.ObtainReadData();
        const realdp *Hptr = H.ObtainReadData();
        const realdp *dGptr = dG.ObtainReadData();
        const realdp *dHptr = dH.ObtainReadData();
        realdp *tmp1 = new realdp [2 * nz];
        realdp *tmp2 = tmp1 + nz;
        ProjecOmegaGHT(dGptr, Hptr, m, n, r, ir, jc, jcc, nz, tmp1);
        ProjecOmegaGHT(Gptr, dHptr, m, n, r, ir, jc, jcc, nz, tmp2);
        axpy_(&nz, &GLOBAL::DONE, tmp1, &GLOBAL::IONE, tmp2, &GLOBAL::IONE);
        SparseMatrix SM2(m, n, ir, jc, jcc, tmp2, nz);
        
        result->NewMemoryOnWrite(); /*For multimanifolds, it is important to new memory at the beginning!*/
        result->GetElement(0).AlphaABaddBetaThis(1, SM2, GLOBAL::N, H, GLOBAL::N, 0);
        result->GetElement(0).AlphaABaddBetaThis(1, SM1, GLOBAL::N, dH, GLOBAL::N, 1);
        result->GetElement(1).AlphaABaddBetaThis(1, SM2, GLOBAL::T, G, GLOBAL::N, 0);
        result->GetElement(1).AlphaABaddBetaThis(1, SM1, GLOBAL::T, dG, GLOBAL::N, 1);
//        result->GetElement(0) = SM2 * H + SM1 * dH;
//        result->GetElement(1) = (G.GetTranspose() * SM2 + dG.GetTranspose() * SM1).GetTranspose();
        
        delete [] tmp1;
        
        SM1.Setnullptr();
        SM2.Setnullptr();
        return *result;
	};
}; /*end of ROPTLITE namespace*/
