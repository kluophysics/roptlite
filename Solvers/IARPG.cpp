
#include "Solvers/IARPG.h"

/*Define the namespace*/
namespace ROPTLITE{

	IARPG::IARPG(const Problem *prob, const Variable *initialx)
	{
        Initialization(prob, initialx);
	};

    IARPG::~IARPG()
    {
    };

	void IARPG::SetProbX(const Problem *prob, const Variable *initialx)
	{
        SolversNSMPGLS::SetProbX(prob, initialx);
        z1 = *initialx;
        z2 = *initialx;
        gfz1 = Prob->GetDomain()->GetEMPTY();
        gfz2 = Prob->GetDomain()->GetEMPTY();
        SolversNSMPGLS::SetProbX(prob, initialx);
        y1 = *initialx;
        y2 = *initialx;
        gfy1 = Prob->GetDomain()->GetEMPTY();
        gfy2 = Prob->GetDomain()->GetEMPTY();
	};

	void IARPG::SetDefaultParams()
	{
        SolversNSMPGLS::SetDefaultParams();
        Variant = LSPG_REGULAR;
        SGIterGap = 5;
		SolverName.assign("IARPG");
	};

    void IARPG::CheckParams(void)
    {
        SolversNSMPGLS::CheckParams();

        char YES[] = "YES";
        char NO[] = "NO";
        char *status;

        printf("IARPG METHODS PARAMETERS:\n");
        status = (SGIterGap > 0) ? YES : NO;
        printf("SGIterGap     :%15d[%s],\n", SGIterGap, status);
    };

    void IARPG::PrintInfo(void)
    {
        printf("i:%d,f:%.10e,df/f:%.3e,", iter, f2, ((f1 - f2) / std::fabs(f2)));

        printf("|nd|:%.3e,si:%.2e,t0:%.2e,t:%.2e,s0:%.2e,s:%.2e,time:%.2e,", ndir1, s1, initiallength, stepsize, initialslope, newslope,  static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

        printf("nf:%d,ng:%d,", nf, ng);
        
        if (nH != 0)
            printf("nH:%d,", nH);
        
        printf("nR:%d,", nR);
        
        if (nV != 0)
            printf("nV(nVp):%d(%d),", nV, nVp);
        
        printf("LSstatus:%s,", LSstatusSetnames[LSstatus].c_str());
        
        printf("SMlambda:%.1e,SMtol:%.1e,SMiter:%d,SMCGiter:%d,PMiter:%d,", SMlambda, SMtol, SMiter, SMCGiter, PMiter);
        printf("\n");
    };

    void IARPG::PrintFinalInfo(void)
    {
        if (Verbose >= FINALRESULT)
        {
            printf("Iter:%d,f:%.10e,ndir:%.3e,|nd|/|nd0|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", iter, f1, ndir1, ndir1/ndir0, ComTime, nf, ng, nR);
            if (nH != 0)
            {
                printf("nH:%d,", nH);
            }
            if (nV != 0)
            {
                printf("nV(nVp):%d(%d),", nV, nVp);
            }
            printf("SMlambda:%.1e,SMtol:%.1e,totalSMiter:%d,totalSMCGiter:%d,totalPMiter:%d,numrestart:%d,", SMlambda, SMtol, totalSMiter, totalSMCGiter, totalPMiter, numrestart);
            printf("\n");
        }
    };

    void IARPG::Run(void)
    {
        Variable xTemp;
        Vector gfTemp;
        
        SolversNSM::Run();
        
        numrestart = 0;
        LSstatus = LSPG_SUCCESS;
        f1 = Prob->f(x1) + Prob->g(x1); nf++;
        f2 = f1;
        Prob->Grad(x1, &gf1); ng++;
        ndir0 = 0; ndir1 = 1;
        newslope = 0;
        y1 = x1; gfy1 = gf1; fy1 = f1;
        z1 = x1; gfz1 = gf1; fz1 = f1;
        
        s1 = 1;
        iter = 0;
        
        Vector initialstepsize(1 + Max_Iteration), acceptedstepsize(1 + Max_Iteration);
        realdp *timeSeriesptr = timeSeries.ObtainWritePartialData();
        realdp *funSeriesptr = funSeries.ObtainWritePartialData();
        realdp *dirSeriesptr = dirSeries.ObtainWritePartialData();
        realdp *initialstepsizeptr = initialstepsize.ObtainWritePartialData();
        realdp *acceptedstepsizeptr = acceptedstepsize.ObtainWritePartialData();
        
        if (Verbose >= FINALRESULT)
            printf("i:%d,f:%.3e,\n", iter, f1);
            
        if (Verbose >= ITERRESULT)
        {
            timeSeriesptr[iter] = static_cast<realdp > (getTickCount() - starttime) / CLK_PS;
            funSeriesptr[iter] = f1;
            dirSeriesptr[iter] = 0;
            initialstepsizeptr[iter] = 0;
            acceptedstepsizeptr[iter] = 0;
        }
        
        bool isstop = IsStopped();
        realdp Wadavalue = 0;
        Vector Weight;

        /*Start the loop*/
        while ((((! isstop) && iter < Max_Iteration) || iter < Min_Iteration) && LSstatus == LSPG_SUCCESS)
        {
            if(iter % SGIterGap == 0)
            {
//                 std::cout << "h1" << std::endl;//--
//                printf("time7:%f\n", static_cast<realdp>(getTickCount() - starttime) / CLK_PS);

                Prob->PreConditioner(x1, Weight, &Weight);
                Weight.ScalarTimesThis(adavalue); /*Weight = adavalue * Weight*/
                Wadavalue = Weight.Onenorm() / Weight.Getlength() * adavalue;
                Mani->TangentSpaceProximalMap(x1, gf1, SMtol, SMlambda, ProxMapType, iter, Weight, Prob, &initD, &SMiter, &SMCGiter, &PMiter, &eta1);
//                if(SMiter >= 20)
//                {
//                    printf("SSN early termination\n");
//                    break;
//                }
//                printf("time8:%f\n", static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
//                break;//---
                totalSMiter += SMiter;
                totalSMCGiter += SMCGiter;
                totalPMiter += PMiter;
                ndir1 = sqrt(Mani->Metric(x2, eta1, eta1)) * Wadavalue;
                if(ndir0 == 0)
                { /*initialize ndir0*/
                    ndir0 = ndir1;
                }
                
                initialslope = Mani->Metric(x1, gf1, eta1);
                /*Initial step size is 1*/
                stepsize = 1;
                initiallength = stepsize;

                /*Start a line search algorithm.*/
                LinesearchArmijo();
                if(LSstatus == LSPG_MINSTEPSIZE)
                { /*If line search reaches the minimum step size, then solving the proximal mapping more accurate. */
                    if(SMlambda == 1e-10 && SMtol == 1e-24)
                    {
                        break;
                    }
                    SMtol = ((SMtol * 0.1 > 1e-10) ? SMtol * 0.1 : 1e-10);
                    SMlambda = ((SMlambda * 0.1 > 1e-10) ? SMlambda * 0.1 : 1e-10);
                    LSstatus = LSPG_SUCCESS;
//                 std::cout << "h2" << std::endl;//--
                    continue;
                }
                
                if(Variant == LSPG_ADALIPSCHITZ)
                {
                    if(stepsize == initiallength)
                    {
                        adavalue /= 1.01;
                    }
                    else
                    {
                        adavalue = ((adavalue * 1.01 < 1.0) ? adavalue * 1.01 : 1.0);
                    }
                }
                /*Output debug information if necessary.*/
                if (LSstatus < LSPG_SUCCESS && Verbose > FINALRESULT)
                {
                    printf("Linesearch fails! LSstatus:%s\n", LSstatusSetnames[LSstatus].c_str());
                }

                if (std::isnan(f2) || std::isinf(f2)) /*Stop when got a nan or inf*/
                {
                    printf("New function value is either nan or inf. Stop!\n");
                    break;
                }
                
                /*if f2 is smaller than the function value at new iterate z1, then safeguard takes effect and restart.*/
                if(iter != 0 && f2 < fz2)
                {
                    y1 = x2;
                    gfy1 = gf2;
                    fy1 = f2;
                    z1 = x2;
                    s1 = 1;
                    if (Verbose >= ITERRESULT)
                    {
                        numrestart++;
//                        std::cout << "restart takes effect!" << std::endl;
                    }
                }
                /*update the safeguard iterate*/
                x1 = z1;
                f1 = Prob->f(x1) + Prob->g(x1); nf++;
                Prob->Grad(x1, &gf1); ng++;
            }

//            printf("time8:%f\n", static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
//                 std::cout << "h3" << std::endl;//--
            Prob->PreConditioner(y1, Weight, &Weight);
            Weight.ScalarTimesThis(adavalue); /*Weight = adavalue * Weight*/
//            Wadavalue = Weight.Onenorm() / Weight.Getlength() * adavalue;
            Mani->TangentSpaceProximalMap(y1, gfy1, SMtol, SMlambda, ProxMapType, iter, Weight, Prob, &initDy, &SMiter, &SMCGiter, &PMiter, &eta1);
//            if(SMiter >= 20)
//            {
//                printf("SSN early termination\n");
//                break;
//            }
//            printf("time9:%f\n", static_cast<realdp>(getTickCount() - starttime) / CLK_PS);
            totalSMiter += SMiter;
            totalSMCGiter += SMCGiter;
            totalPMiter += PMiter;
            
            Mani->Retraction(y1, eta1, &z2);
            fz2 = Prob->f(z2) + Prob->g(z2); nf++;
            s2 = (1.0 + std::sqrt(1.0 + 4.0 * s1 * s1)) / 2.0;
            Vector zeta(gf1); Mani->InvRetraction(z2, z1, &zeta);
            Mani->ScalarTimesVector(z2, -(s1 - 1) / s2, zeta, &zeta);
            
            Mani->Retraction(z2, zeta, &y2);
            
            fy2 = Prob->f(y2) + Prob->g(y2); nf++;
            Prob->Grad(y2, &gfy2); ng++;
            
//            if(iter % 100 == 0)
//                std::cout << "iter:" << iter << ", norm(z1-z2):" << (z1 - z2).Fnorm() << std::endl;//-----
            /*Call the function to check whether the stopping criterion is satisfied or not.
            The default function is written in Solvers.h and Solvers.cpp*/
            isstop = IsStopped();
            
            iter++;

            if (Verbose >= ITERRESULT)
            {
                /*Output information*/
                if (iter % OutputGap == 0)
                {
                    PrintInfo();
                }
                /*Store debug information in the arrays*/
                timeSeriesptr[iter] = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
				funSeriesptr[iter] = (fz2 < f2) ? fz2 : f2; dirSeriesptr[iter] = ndir1;
                acceptedstepsizeptr[iter] = stepsize;
                initialstepsizeptr[iter] = initiallength;
            }
            
            xTemp = z1; z1 = z2; z2 = xTemp;
            gfTemp = gfz1; gfz1 = gfz2; gfz2 = gfTemp;
            xTemp = y1; y1 = y2; y2 = xTemp;
            gfTemp = gfy1; gfy1 = gfy2; gfy2 = gfTemp;
            s1 = s2;
        }
        
        ComTime = static_cast<realdp>(getTickCount() - starttime) / CLK_PS;
        if (Verbose >= ITERRESULT)
        {
            lengthSeries = iter + 1;
            SolverInfo.AddToFields("initialstepsize", initialstepsize);
            SolverInfo.AddToFields("acceptedstepsize", acceptedstepsize);
        }
        PrintFinalInfo();
    };

    void IARPG::SetParams(PARAMSMAP params)
    {
        SolversNSMPGLS::SetParams(params);
        PARAMSMAP::iterator iter;
        for (iter = params.begin(); iter != params.end(); iter++)
        {
            if (iter->first == static_cast<std::string> ("SGIterGap"))
            {
                SGIterGap = static_cast<integer> (static_cast<integer> (iter->second));
            }
        }
    };
}; /*end of ROPTLITE namespace*/
