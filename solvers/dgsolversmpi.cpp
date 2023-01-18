#ifndef __DGSOLVERSMPI
#define __DGSOLVERSMPI

#include "pDLAK.cpp"
#include "DIRKcoeff.cpp"
#include "solutionupdate.cpp"
//#include "newtoninit.cpp"
//#include "../application/FM/getAVfield.cpp"
//#include "../application/FM/nsNew/computeDynamicSmagConstant.cpp"
//#include "../utilities/computeForces.cpp"
//#include "../utilities/computeAvgSolution.cpp"

// Written by: C. Nguyen & P. Fernandez

void computeOrdering(sysstruct &sys, Int preconditioner);

void computePreconditioner(sysstruct &sys, Int preconditioner);

void gmresmpi(sysstruct &sys, Int *convFlag);

void recomputeJacobian(sysstruct &sys, appstruct &app)
{
    if (sys.my_rank == 0 && app.reuseJacobian == 1)
        printf("\n\n*** The Jacobian matrix will be updated in the next iteration.***\n\n");
    
    app.reuseJacobian = 0;
    sys.linearSolvesClocks = 0;
}

void decideIfReuseJacobian(sysstruct &sys, appstruct &app)
{
    // Make decision:
    Int reuseJacobian;
    if (sys.my_rank == 0) {
        if (sys.linearSolvesClocks < sys.lastAssemblyAndPrecClocks)
            reuseJacobian = 1;
        else
            reuseJacobian = 0;
    }

    // Broadcast decision:
    Int root = 0, count = 1;
    MPI_Bcast(&reuseJacobian, count, MPI_INT, root, MPI_COMM_WORLD);

    // Update structures:
    if ((int) reuseJacobian == 1)
        app.reuseJacobian = 1;
    else
        recomputeJacobian(sys, app);
}

void decideIfAdaptGMREStol(sysstruct &sys, appstruct app, Int iter)
{
    if (sys.adaptGMREStol == 1 && app.reuseJacobian == 0 && iter > 1)
        sys.adaptGMREStol = 0;
    else
        sys.adaptGMREStol = 1;
}

void solveLinearProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{
    clock_t t;
    double assemblyTime, precondTime, orderingTime, GMREStime;

    // Assembly linear system / RHS
//     if (app.quasiNewton == 1)
//         MPI_Barrier( MPI_COMM_WORLD );
    t = clock();
    if (app.reuseJacobian == 0) {
        app.reusePreconditioner = 0;
        if (sys.my_rank == 0)
            printf("Assembling linear system...\n");
        assembleLinearSystemMPI(sys, elems, mesh, master, sol, app, temps);
        if (sys.print >= 1 || app.quasiNewton == 1)
            MPI_Barrier( MPI_COMM_WORLD );
        assemblyTime = (((clock() - t)*1.0e3)/CLOCKS_PER_SEC);
        if (app.quasiNewton == 1)
            sys.lastAssemblyAndPrecClocks = (long) (clock() - t);
    }
    else {
        app.reusePreconditioner = 1;
        if (sys.my_rank == 0)
            printf("Assembling RHS vector...\n");
        assembleRHSMPI(sys, elems, mesh, master, sol, app, temps);
//         if (sys.print >= 1)
//             MPI_Barrier( MPI_COMM_WORLD );
        assemblyTime = (((clock() - t)*1.0e3)/CLOCKS_PER_SEC);
    }
    
    // Compute ordering (if necessary)
    t = clock();
    if (app.reuseOrdering == 0) {
        if (sys.my_rank == 0)
            printf("Computing ordering...\n");
        computeOrdering(sys, sys.preconditioner);
        app.reuseOrdering = 1;
    }
    if (sys.print >= 1 || app.quasiNewton == 1)
        MPI_Barrier( MPI_COMM_WORLD );
    orderingTime = (((clock() - t)*1.0e3)/CLOCKS_PER_SEC);

    // Compute preconditioner (if necessary)
    t = clock();
    if (app.reusePreconditioner == 0) {
        if (sys.my_rank == 0)
            printf("Computing preconditioner...\n");
        computePreconditioner(sys, sys.preconditioner);
    }
    if (sys.print >= 1 || app.quasiNewton == 1)
        MPI_Barrier( MPI_COMM_WORLD );
    precondTime = (((clock() - t)*1.0e3)/CLOCKS_PER_SEC);
    if (app.quasiNewton == 1 && app.reusePreconditioner == 0)
        sys.lastAssemblyAndPrecClocks += (long) (clock() - t);

    // Solve linear system:
    t = clock();
    if (sys.my_rank == 0)
        printf("Solving linear system...\n");
    gmresmpi(sys, convFlag);
//     if (sys.print >= 1 || app.quasiNewton == 1)
//         MPI_Barrier( MPI_COMM_WORLD );
    GMREStime = (((clock() - t)*1.0e3)/CLOCKS_PER_SEC);
    if (app.quasiNewton == 1)
        sys.linearSolvesClocks += (long) (clock() - t);

    if (sys.print >= 1 && sys.my_rank == 0) {
        printf("\n\nBREAKDOWN OF CPU TIMES:\n");
        printf("1. Matrix assembly: %g ms\n", assemblyTime);
        printf("2. Compute preconditioner: %g ms\n", precondTime);
        printf("3. Compute ordering: %g ms\n", orderingTime);
        printf("4. GMRES: %g ms\n\n", GMREStime);
    }
}

void solveNonlinearProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{
    Int i, j, n, ri, inc = 1;
    Int convLinearFlag[1];
    double alphaMin = 1.0e-4;
    double rNorm = sys.NewtonTol + 1.0, alpha, oldNorm;
    
    /* Newton iteration */
    sys.robustMode = 0;
    Int iter = 0, trueNewtonIter = 0;
    while (rNorm > sys.NewtonTol && iter < sys.NewtonMaxiter && trueNewtonIter < sys.trueNewtonMaxiter) {
        iter += 1;
        if (sys.my_rank == 0)
            printf("\nNewton iteration:  %d\n", iter);
        if (app.reuseJacobian == 0)
            trueNewtonIter += 1;
        
        if (sys.linearProblem == 0 && sys.adaptiveGMREStol == 1)
            decideIfAdaptGMREStol(sys, app, iter);
                
        /* Solve linearized problem */
        sys.rNorm = rNorm;
        solveLinearProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convLinearFlag);
        MPI_Allreduce(&sol.rNorm, &oldNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        oldNorm = sqrt(oldNorm);
        if ((isnan(oldNorm) || isinf(oldNorm)) && sys.my_rank == 0) {
            error("\nResidual norm in current Newton's iterate is NaN or Inf.\n");
        }
        
//         if (sys.my_rank==0) {
//             print2darray(&sys.x[0],app.nch,3*mesh.npfmax);
//             print2darray(&sol.UH[0],app.nch,3*mesh.npfmax);
//             print2darray(&sol.UDG[0],mesh.npemax,app.nch);
//         }
        
        if (*convLinearFlag == -1)      // Failure (GMRES breakdown or orthogonalization did not converge)
            break;
        
        // Decide whether or not to reuse the Jacobian (only applies for quasi-Newton):
        if (app.quasiNewton == 1)
            decideIfReuseJacobian(sys, app);
        
        // Check if initial guess is below Newton tolerance:
        if (iter == 1 && oldNorm < sys.NewtonTol) {
            rNorm = oldNorm;
            sys.rNorm = rNorm;
            if (sys.my_rank == 0)
                printf("Old residual: %g,   New residual: %g    \n", oldNorm, rNorm);
            break;
        }
        
//         printf("Old residual: %g   \n", oldNorm);
//         computeResidualNormMPI(sys, elems, mesh, master, sol, app, temps, ndims);                
//         MPI_Allreduce(&sol.rNorm, &rNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                  
//         rNorm = sqrt(rNorm);
//         printf("New residual: %g    \n", rNorm);    
//         error("here");
        
        /* Perform line search */
        sol.alpha = 1.0;
        while (1) {
            /* update solution */
            if (sol.alpha < 0.999999) 
                updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha);
            else
                updateSolMPI(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha);   
            
//             if (sys.my_rank==0) {
//                 print2darray(&sys.x[0],app.nch,3*mesh.npfmax);
//                 print2darray(&sol.UH[0],app.nch,3*mesh.npfmax);
//                 print2darray(&sol.UDG[0],mesh.npemax,app.nch);
//             }
//             error("here");
            
            /* Compute the residual norm */            
            computeResidualNormMPI(sys, elems, mesh, master, sol, app, temps, ndims);                
            MPI_Allreduce(&sol.rNorm, &rNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                  
            rNorm = sqrt(rNorm);
            if ((isnan(rNorm) || isinf(rNorm)) && sys.my_rank == 0) {
                printf("\n\n***********************\n");
                printf("Residual norm in line search of Newton's method is NaN or Inf.\n");
                printf("\n***********************\n\n");
            }
            
            if (isnan(rNorm) || isinf(rNorm) || (rNorm > oldNorm && rNorm > sys.NewtonTol)) {
                alpha = fabs(sol.alpha/2);
                sol.alpha = -alpha;
                if (sys.my_rank == 0)
                    printf("Damped Newton coefficient: alpha = %g || Residual = %g || Old residual = %g\n", alpha, rNorm, oldNorm);
                if (app.quasiNewton == 1)
                    recomputeJacobian(sys, app);
                if (alpha < alphaMin) {
                    if (sys.my_rank == 1)
                        printf("Damped Newton does not converge for alpha > %g\n.", alphaMin);
                    if (sys.robustMode == 0 || app.quasiNewton == 1) {
                        // Go back to previous Newton iteration and try with a more robust version of the code:
                        sys.robustMode = 1;     // The code will remain in robust mode for the entire nonlinear solve
                        sol.alpha = - 2.0 * fabs(sol.alpha);
                        updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha);
//                         if (app.AVflag == 6 || app.AVflag == 8 || app.AVflag == 9 || app.AVflag == 10) {
//                             restartNewton = 1;
//                             app.AVflag = 0;
//                             Int avgType = 2;
// //                            getAVfield(&sol.avField_p1CG[0], &sol.UDGi[0], &mesh.dgnodes[0], sys, mesh, master, app, temp, ndims, avgType, app.AVflag);
//                         }
                    }
                    break;
                }
            }
            else
                break;
        }
        if (sys.my_rank == 0)
            printf("Old residual: %g,   New residual: %g    \n", oldNorm, rNorm);
    }

    if (rNorm < sys.NewtonTol)
        *convFlag = 1;
    else {
        *convFlag = 0;
        if (app.quasiNewton == 1)
            recomputeJacobian(sys, app);
        if (sys.my_rank == 0)
            printf("Newton's method did not converge to tolerance %g in %d iterations (%d true Newton itererations).\n", sys.NewtonTol, iter, trueNewtonIter);
    }
}

void solveSteadyProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{
    if (app.linearProblem == 0)
        solveNonlinearProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
    else if (app.linearProblem == 1)
        solveLinearProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
        sol.alpha=1;
        updateSolMPI(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha);
        /* TODO: Need to update solution (updateSol) and decide if the Jacobian is reused (the latter only applies for quasi-Newton) */
    //else if (app.linearProblem == 2)
    //    acceleratedResidualDescent(sys, elems, mesh, master, sol, app, temps, ndims);
}

void LPcoeff(double * an, double * tn, double t, Int k) 
{
    Int i, j;
    
    for (i=0; i<k; i++) {
        an[i] = 1.0;
        for (j=0; j<k; j++)
            if (i != j)
                an[i] = an[i]*((t-tn[j])/(tn[i]-tn[j]));
    }
}

void BDFsolveUnsteadyProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
//     error("BDF solver seems to work but has not been validated properly yet.\n");
    
    Int i, j, k, m, ie, jstar, numSolvesAttemped;
    Int ne = ndims[5];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int npv = ndims[9];
    Int szu = sol.UDG.size();
    Int szh = sol.UH.size();
    Int BDFsteps_tmp;
    Int BDFsteps = app.BDFsteps;
    Int nein = sys.elempartpts[0] + sys.elempartpts[1];
    Int nfin = sys.entpartpts[0] + sys.entpartpts[1];
    Int bsz = sys.blkSize;
    Int convFlag[1];
    clock_t t1, tNewton;
    
    Int sza = 3;                                // Number of snapshots for Minimal Residual algorithm
    Int saveSolFreq = 1000;                       // Every how many time-steps the solution is saved to a file
    Int recomputeOrderingFreq = 50;             // Every how many time-steps the ordering is recomputed
    Int writeQflag = 1;                         // 0: Only u is written to a binary file. 1: Both u and q are written to a binary file
    Int writeAvgSolution = 0;                   // 0: The time average solution is NOT saved to a file. 1: The time average solution is saved to a file
    Int timeStepToStartAvg = 1000;              // Time step to start averaging the solution
    string timeStepsFileName = app.filein + "_timeSteps.out";
    
    /* Initialize structures */
    sol.UDGi.resize(szu);
    sol.UHi.resize(szh);
    if (app.wave==1)
        error("BDF schemes not implemented for wave problems.\n");
    vector<double> UDGprev(BDFsteps * szu);
    vector<double> UHprev(BDFsteps * szh);
    if (app.wave == 1)
        error("BDF schemes not implemented for wave problems.\n");
    vector<double> UDGpre(szu*sza);
    vector<double> UHpre(szh*sza);
    vector<double> a(sza);
    sys.Ru_MR.resize(nein*npv*ncu);
    sys.dRuda_MR.resize(sza*nein*npv*ncu);
    sys.a2_MR.resize(sza);
    sys.a2s_MR.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.a2s_MR[i].resize(sys.a2_MR.size());
    sys.da_MR.resize(sza);
    sys.C_MR.resize(sza*sza);
    sys.Clocal_MR.resize(sza*sza);
    
    double time = 0.0, dt;
    for (i=0; i < app.dt.size(); i++) {
        if (sys.my_rank==0)
            printf("Timestep :  %d,  Time : %g,  dt:  %g\n", i+1, time, app.dt[i]);
        
        t1 = clock();
        dt = app.dt[i];
        app.time = time + dt;
        BDFsteps_tmp = min(i+1 , BDFsteps);     // BDF scheme to be used in current time step.
        
        /* get solution from the previous timestep */
        for (j = 0; j< sol.UDG.size(); j++)
            sol.UDGi[j] = sol.UDG[j];
        for (j = 0; j< sol.UH.size(); j++)
            sol.UHi[j] = sol.UH[j];
        if (app.wave == 1)
            error("BDF schemes not implemented for wave problems.\n");
        
        for (k = 0; k < BDFsteps_tmp-1; k++) {
            for (j = 0; j< szu; j++)
                UDGprev[(BDFsteps_tmp-1-k)*szu + j] = UDGprev[(BDFsteps_tmp-1-k-1)*szu + j];
            for (j = 0; j< szh; j++)
                UHprev[(BDFsteps_tmp-1-k)*szh + j] = UHprev[(BDFsteps_tmp-1-k-1)*szh + j];
            if (app.wave == 1)
                error("BDF schemes not implemented for wave problems.\n");
        }
        for (j = 0; j< szu; j++)
            UDGprev[j] = sol.UDG[j];
        for (j = 0; j< szh; j++)
            UHprev[j] = sol.UH[j];
        if (app.wave == 1)
            error("BDF schemes not implemented for wave problems.\n");
        
//         error("Code properly some line in BDFsolveUnsteadyProblemMPI.\n");
//         if (app.AVflag == 6 || app.AVflag == 8 || app.AVflag == 9 || app.AVflag == 10) {
//             Int avgType = 2;
//             getAVfield(&sol.avField_p1CG[0], &sol.UDGi[0], &mesh.dgnodes[0], sys, mesh, master, app, temps[0], ndims, avgType, app.AVflag);
//         }
//         if (app.SGSmodel == 4) {
//             Int minimizationTarget = 0;
//             Int projPorder = app.porder-1;
//             computeDynamicSmagConstant(&sol.Cs[0], sys, mesh, master, app, &sol.UDGi[0], &mesh.hAvg[0], projPorder, &ndims[0], minimizationTarget);
//         }
        
        if (i == 0 || (i % recomputeOrderingFreq) == 0) {
            app.reuseOrdering = 0;
        }
        
        /* Compute the source term */
        switch (BDFsteps_tmp) {
            case 1:
                app.fc_u = 1.0 / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = UDGprev[0*szu+j] / dt;
                break;
            case 2:
                app.fc_u = (3.0/2.0) / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = ( (4.0/2.0) * UDGprev[0*szu+j] - (1.0/2.0) * UDGprev[1*szu+j] ) / dt;
                break;
            case 3:
                app.fc_u = (11.0/6.0) / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = ( (18.0/6.0) * UDGprev[0*szu+j] - (9.0/6.0) * UDGprev[0*szu+j] 
                                + (2.0/6.0) * UDGprev[0*szu+j] ) / dt;
                break;
            case 4:
                app.fc_u = (25.0/12.0) / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = ( (48.0/12.0) * UDGprev[0*szu+j] - (36.0/12.0) * UDGprev[1*szu+j] 
                                + (16.0/12.0) * UDGprev[2*szu+j] - (3.0/12.0)  * UDGprev[3*szu+j] ) / dt;
                break;
            case 5:
                app.fc_u = (137.0/60.0) / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = ( (300.0/60.0) * UDGprev[0*szu+j] - (300.0/60.0) * UDGprev[1*szu+j] 
                                + (200.0/60.0) * UDGprev[2*szu+j] - (75.0/60.0)  * UDGprev[3*szu+j] 
                                + (12.0/60.0) * UDGprev[4*szu+j] ) / dt;
                break;
            case 6:
                app.fc_u = (147.0/60.0) / dt;
                for (j=0; j<szu; j++)
                    sol.SH[j] = ( (360.0/60.0) * UDGprev[0*szu+j] - (450.0/60.0) * UDGprev[1*szu+j]
                                + (400.0/60.0) * UDGprev[2*szu+j] - (225.0/60.0) * UDGprev[3*szu+j]
                                + (72.0/60.0) * UDGprev[4*szu+j]  - (10.0/60.0)  * UDGprev[5*szu+j] ) / dt;
                break;
            default:
                error("BDF scheme not implemented yet.");
        }
        
        if (app.wave == 1)
            error("BDF schemes not implemented for wave problems.\n");
        else
            app.fc_q = 1;
        
        /* Set q source term to zero for non-wave problems */
        if (app.wave == 0) {
            for (m = 0; m < ne; m++)
                for (j = ncu; j < nc; j++)
                    for (k = 0; k < npv; k++)
                        sol.SH[k + j*npv + m*npv*nc] = 0.0;
        }
        
        jstar = (i % sza);
        for (j=0; j<szu; j++)
            UDGpre[jstar*szu+j] = sol.UDG[j];
        for (j=0; j<szh; j++)
            UHpre[jstar*szh+j] = sol.UH[j];
        
        /* Compute the initial guess */
        if (i+1 >= sza) {
            for (j=0; j<sza; j++)
                a[j] = 0.0;
            a[jstar] = 1.0;
            tNewton = clock();
            // NewtonInitMPI(sys, elems, mesh, master, sol, app, temps, ndims, &UDGpre[0], &UHpre[0], &a[0], sza);
            if (sys.my_rank == 0)
                printf("Time to compute initial guess with minimal residual algorithm: %g ms\n", ((clock()-tNewton)*1.0e3)/CLOCKS_PER_SEC);
        }
        
        /* Solve the steady problem */
        *convFlag = 0;
        numSolvesAttemped = 0;
        while (*convFlag != 1) {
            solveSteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
            if (*convFlag != 1) {
                numSolvesAttemped ++;
                if (i+1 < sza || numSolvesAttemped >= 2) {
                    app.BDFnotConverged[i] = 1;

                    if (sys.my_rank == 0)
                        printf("\n\nWARNING: Unsteady solve at time step %d did not converge.\nThe time-step size will be reduced by half.\n",i+1);
//                         error("\n");

                    /////////// HACK START ///////////
                    for (j=0; j<szu; j++)
                        sol.UDG[j] = UDGpre[jstar*szu+j];
                    for (j=0; j<szh; j++)
                        sol.UH[j] = UHpre[jstar*szh+j];

                    app.dt[i+1] = app.dt[i+1] / 2.0;
                    break;
                    /////////// HACK END ///////////
                }
                // Try without minimal residual algorithm for initial guess:
                if (sys.my_rank == 0)
                    printf("\n\nAttemping nonlinear solve without Minimal Resigual algorithm for initial guess.\n");
                for (j=0; j<szu; j++)
                    sol.UDG[j] = UDGpre[jstar*szu+j];
                for (j=0; j<szh; j++)
                    sol.UH[j] = UHpre[jstar*szh+j];
            }
        }
        
        if (app.wave == 1)
            error("BDF schemes not implemented for wave problems.\n");

        time = time + dt;
        
        // Compute average solution:
        //if ((i+1 >= timeStepToStartAvg) && writeAvgSolution == 1)
            //computeAvgSolution(sol, (i+1)-(timeStepToStartAvg-1));
        
        // Write solution to file:
        if (((i+1) % saveSolFreq) == 0) {
            string filename = app.fileout + "_t" + NumberToString(i+1) + "_np" + NumberToString(sys.my_rank) + ".bin";
            sol.writeSol2File(filename, nein, bsz*nfin, ndims, writeQflag);
            
            //if ((i+1 >= timeStepToStartAvg) && (writeAvgSolution == 1)) {
            //    string filename_avg = app.fileout + "_avg_np" + NumberToString(sys.my_rank) + ".bin";
            //    sol.writeAvgSol2File(filename_avg, nein, bsz*nfin, ndims, writeQflag);
            //}
        }
        
        if (sys.my_rank == 0) {
            printf("\n\nTOTAL TIME TO SOLVE TIME-STEP NO. %d: %g ms\n\n", i+1, ((clock() - t1)*1.0e3)/CLOCKS_PER_SEC);
            writeTimeStepSize2File(timeStepsFileName, i+1, i+1, app.time, dt);
        }
    }

    // Report time-steps in which nonlinear solver did not converge:
    if (sys.my_rank == 0) {
        printf("\n\nLIST OF TIME-STEP IN WHICH NONLINEAR SOLVER DID NOT CONVERGE:\n");
        for (i = 0; i < app.BDFnotConverged.size(); i++) {
            if (app.BDFnotConverged[i] == 1)
                printf("Time-step No. %d.\n",i+1);
        }
    }
}

void DIRKsolveUnsteadyProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
    Int i, j, k, m, istage, ie, jstar, numSolvesAttemped;
    Int ne = ndims[5];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int npv = ndims[9];
    Int dirkStage = app.dirkStage;
    Int dirkOrder = app.dirkOrder;
    Int nein = sys.elempartpts[0]+sys.elempartpts[1];
    Int nfin = sys.entpartpts[0]+sys.entpartpts[1];
    Int bsz = sys.blkSize;
    Int convFlag[1];

    clock_t t1, t2, tNewton;
 
    /* Get DIRK coefficients */    
    double *dirkc = new double[dirkStage];
    double *dirkd = new double[dirkStage*dirkStage];
    double *dirkt = new double[dirkStage];    
    DIRKcoeff(dirkc, dirkd, dirkt, dirkStage, dirkOrder);
    
    /* Initialize the previous timestep solution */
    sol.UDGi.resize(sol.UDG.size());
    sol.UHi.resize(sol.UH.size());
    if (app.wave==1) {
        sol.PDGi.resize(sol.PDG.size());
        sol.SP.resize(sol.PDG.size());
    }

    /* Initial guess for nonlinear problem */
    Int sza = 3;
    Int szu = sol.UDG.size();
    Int szh = sol.UH.size();
    vector<double> UDGpre(szu*sza);
    vector<double> UHpre(szh*sza);
    vector<double> a(sza);
    sys.Ru_MR.resize(nein*npv*ncu);
    sys.dRuda_MR.resize(sza*nein*npv*ncu);
    sys.a2_MR.resize(sza);
    sys.a2s_MR.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.a2s_MR[i].resize(sys.a2_MR.size());
    sys.da_MR.resize(sza);
    sys.C_MR.resize(sza*sza);
    sys.Clocal_MR.resize(sza*sza);
    
    Int saveSolFreq = 200;                       // Every how many time-steps the solution is saved to a file
    //Int saveSolFreq = app.flag[28];
    Int recomputeOrderingFreq = 50;             // Every how many time-steps the ordering is recomputed
    Int writeQflag = 1;                         // 0: Only u is written to a binary file. 1: Both u and q are written to a binary file
    Int writeAvgSolution = 1;                   // 0: The time average solution is NOT saved to a file. 1: The time average solution is saved to a file
    Int timeStepToStartAvg = 1000;              // Time step to start averaging the solution
    double time = 0.0, dt;
    int ii = 0;
    string timeStepsFileName = app.filein + "_timeSteps.out";
    for (i=0; i<app.dt.size(); i++) {
        t1 = clock();
        
        if (sys.my_rank==0)
            printf("Timestep :  %d,  Time : %g,  dt:  %g\n",i+1,time,app.dt[i]);
        
        dt = app.dt[i];
        
        /* get solution from the previous timestep */
        for (j = 0; j< sol.UDG.size(); j++)
            sol.UDGi[j] = sol.UDG[j];
        for (j = 0; j< sol.UH.size(); j++)
            sol.UHi[j] = sol.UH[j];
        if (app.wave == 1)
            for (j = 0; j< sol.PDG.size(); j++)
                sol.PDGi[j] = sol.PDG[j];
//         error("Code properly line 459 in dgsolversmpi.cpp\n");
//         if (app.AVflag == 6 || app.AVflag == 8 || app.AVflag == 9 || app.AVflag == 10) {
//             Int avgType = 2;
//             getAVfield(&sol.avField_p1CG[0], &sol.UDGi[0], &mesh.dgnodes[0], sys, mesh, master, app, temps[0], ndims, avgType, app.AVflag);
//         }
//         if (app.SGSmodel == 4) {
//             Int minimizationTarget = 0;
//             Int projPorder = app.porder-1;
//             computeDynamicSmagConstant(&sol.Cs[0], sys, mesh, master, app, &sol.UDGi[0], &mesh.hAvg[0], projPorder, &ndims[0], minimizationTarget);
//         }
        
        if (i == 0 || (i % recomputeOrderingFreq) == 0) {
            app.reuseOrdering = 0;
        }
        
        /* Loop over each DIRK stage */
        for (istage = 0; istage < app.dirkStage; istage++) {
            t2 = clock();

            if (sys.my_rank==0)
                printf("DIRK stage :  %d\n", istage+1);
            
            /* Update app */
            app.time = time+dt*dirkt[istage];
            app.fc_u = dirkd[istage*dirkStage+istage]/dt;
            if (app.wave==1) 
                app.fc_q = app.fc_u;
            else
                app.fc_q = 1;            
            
            /* Compute the source term */
            switch (istage) {
                case 0:
                    for (j=0; j<sol.UDG.size(); j++)
                        sol.SH[j] = (dirkd[0]/dt)*sol.UDGi[j];
                    break;
                case 1:
                    for (j=0; j<sol.UDG.size(); j++)
                        sol.SH[j] = ((dirkd[1]+dirkd[dirkStage+1])/dt)*sol.UDGi[j] - sol.UDG[j]*dirkd[1]/dt;
                    break;
                case 2:
                    for (j=0; j<sol.UDG.size(); j++)
                        sol.SH[j] = ((dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2])/dt)*sol.UDGi[j] 
                                     - (dirkd[2]/dirkd[1])*((dirkd[1]+dirkd[dirkStage+1])*sol.UDGi[j]/dt-sol.SH[j]) 
                                     - sol.UDG[j]*dirkd[dirkStage+2]/dt;
                    break;
                default:
                    error("DIRK scheme not implemented yet.");
            }

            /* Set q sourcce term to zero for non-wave problems */
            if (app.wave == 0) {
                for (m = 0; m < ne; m++)
                    for (j = ncu; j < nc; j++)
                        for (k = 0; k < npv; k++)
                            sol.SH[k + j*npv + m*npv*nc] = 0.0;
            }
            
            jstar = (ii % sza);
            for (j=0; j<szu; j++)
                UDGpre[jstar*szu+j] = sol.UDG[j];
            for (j=0; j<szh; j++)
                UHpre[jstar*szh+j] = sol.UH[j];
            ii += 1;
            
            /* Compute the initial guess */
            if (ii>=sza) {
                tNewton = clock();
                for (j=0; j<sza; j++)
                    a[j] = 0.0;
                a[jstar] = 1.0;
                //NewtonInitMPI(sys, elems, mesh, master, sol, app, temps, ndims, &UDGpre[0], &UHpre[0], &a[0], sza);
                if (sys.my_rank == 0)
                    printf("Time to compute initial guess with minimal residual algorithm: %g ms\n", ((clock()-tNewton)*1.0e3)/CLOCKS_PER_SEC);
            }
            
            /* Solve the steady problem corresponding to the istage stage */
            *convFlag = 0;
            numSolvesAttemped = 0;
            while (*convFlag != 1) {
                solveSteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
                if (*convFlag != 1) {
                    numSolvesAttemped ++;
                    if (ii < sza || numSolvesAttemped >= 2) {
                        app.DIRKnotConverged[ii-1] = 1;
                        
                        if (sys.my_rank == 0)
                            printf("\n\nWARNING: Unsteady solve at time step %d, DIRK stage %d did not converge.\nThe time-step size will be reduced by half.\n",i+1,istage+1);
//                         error("\n");
                        
                        /////////// HACK START ///////////
                        for (j=0; j<szu; j++)
                            sol.UDG[j] = UDGpre[jstar*szu+j];
                        for (j=0; j<szh; j++)
                            sol.UH[j] = UHpre[jstar*szh+j];
                        
                        if (istage == app.dirkStage - 1)
                            app.dt[i+1] = app.dt[i+1] / 2;
                        else {
                            app.dt[i] = app.dt[i] / 2;
                            dt = app.dt[i];
                        }
                        break;
                        /////////// HACK END ///////////
                    }
                    // Try without minimal residual algorithm for initial guess:
                    if (sys.my_rank == 0)
                        printf("\n\nAttemping nonlinear solve without Minimal Resigual algorithm for initial guess.\n");
                    for (j=0; j<szu; j++)
                        sol.UDG[j] = UDGpre[jstar*szu+j];
                    for (j=0; j<szh; j++)
                        sol.UH[j] = UHpre[jstar*szh+j];
                }
            }
            
            if (app.wave) {
                switch (istage) {
                    case 0:
                        for (j=0; j<sol.PDG.size(); j++)
                            sol.SP[j] = (dirkd[0]/dt)*sol.PDGi[j];
                        break;
                    case 1:
                        for (j=0; j<sol.PDG.size(); j++)
                            sol.SP[j] = ((dirkd[1]+dirkd[dirkStage+1])/dt)*sol.PDGi[j] - sol.PDG[j]*dirkd[1]/dt;
                        break;
                    case 2:
                        for (j=0; j<sol.PDG.size(); j++)
                            sol.SP[j] = ((dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2])/dt)*sol.PDGi[j] 
                                        - (dirkd[2]/dirkd[1])*((dirkd[1]+dirkd[dirkStage+1])*sol.PDGi[j]/dt-sol.SP[j]) 
                                        - sol.PDG[j]*dirkd[dirkStage+2]/dt;
                        break;
                    default:
                        error("DIRK scheme not implemented yet.\n");
                }
                for (j=0; j<ne; j++)
                    for (k=0; k<ncu; k++)
                        sol.PDG[j*ncu+k] = (dt/app.fc_u)*(sol.UDG[j*nc+k]+sol.SP[j*ncu+k]);                 
            }

            if (sys.my_rank == 0) {
                printf("\n\nTOTAL TIME TO SOLVE TIME STAGE NO. %d IN TIME-STEP NO. %d: %g ms\n\n", istage+1, i+1, ((clock() - t2)*1.0e3)/CLOCKS_PER_SEC);
                writeTimeStepSize2File(timeStepsFileName, i+1, istage+1, app.time, dt);
            }
        }
        time = time + dt;
        
        // Compute average solution:
        //if ((i+1 >= timeStepToStartAvg) && writeAvgSolution == 1)
            //computeAvgSolution(sol, (i+1)-(timeStepToStartAvg-1));
        
        // Write solution to file:
        if (((i+1) % saveSolFreq) == 0) {
            string filename = app.fileout + "_t" + NumberToString(i+1) + "_np" + NumberToString(sys.my_rank) + ".bin";
            sol.writeSol2File(filename, nein, bsz*nfin, ndims, writeQflag);
            
            //if ((i+1 >= timeStepToStartAvg) && (writeAvgSolution == 1)) {
            //    string filename_avg = app.fileout + "_avg_np" + NumberToString(sys.my_rank) + ".bin";
            //    sol.writeAvgSol2File(filename_avg, nein, bsz*nfin, ndims, writeQflag);
            //}
        }
        
        if (sys.my_rank == 0)
            printf("\n\nTOTAL TIME TO SOLVE TIME-STEP NO. %d: %g ms\n\n", i+1, ((clock() - t1)*1.0e3)/CLOCKS_PER_SEC);
    }

    // Report time-steps in which nonlinear solver did not converge:
    if (sys.my_rank == 0) {
        printf("\n\nLIST OF TIME-STEP IN WHICH NONLINEAR SOLVER DID NOT CONVERGE:\n");
        for (i = 0; i < app.DIRKnotConverged.size(); i++) {
            if (app.DIRKnotConverged[i] == 1)
                printf("Time-step No. %d, DIRK stage No. %d.\n",(i/app.dirkStage)+1,(i%app.dirkStage)+1);
        }
    }
    
    delete[] dirkd;
    delete[] dirkt;
    delete[] dirkc;
}

void solveUnsteadyProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
    DIRKsolveUnsteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims);
//     BDFsolveUnsteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims);
}

void solveProblemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
    // TODO: Use a more sophisticated strategy for reuseOrdering
    
    clock_t t = clock();
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int convFlag[1];
    
    app.reuseOrdering = 0;
    app.reuseJacobian = 0;
    app.reuseResidual = 0;
    sys.linearSolvesClocks = 0;         // Time spent in linear solves since the Jacobian matrix was computed the last time (only applies for quasi-Newton)

    /* Initialize x */
    for (Int i=0; i<sys.x.size(); i++)
        sys.x[i] = 0.0;
    
    /* Compute Q to make Rq = 0 */
    if (app.flag_q == 1) {
// // //         #pragma omp parallel num_threads(sys.noThreads)
// // //         {
// // //             int this_thread = omp_get_thread_num();
            int this_thread = 0;
// // //             #pragma omp for
            for (Int ie=0; ie<ne; ie++)
                getQ(elems[this_thread], mesh, master, app, sol, temps[this_thread], 
                     &sol.UDG[ie*npv*nc], &sol.UDG[ie*npv*nc+npv*ncu], &sol.UH[0], &sol.SH[ie*npv*nc+npv*ncu], ie);
// // //         }
    }
    
    if (app.tdep == 1)    
        solveUnsteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims);    
    else {
        solveSteadyProblemMPI(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
        if (*convFlag != 1 && sys.my_rank == 0)
            printf("\n\nWARNING: Steady solve did not converge.\n\n");
        
        Int writeQflag = 1; 
        Int nein = sys.elempartpts[0]+sys.elempartpts[1];
        Int nfin = sys.entpartpts[0]+sys.entpartpts[1];
        Int bsz = sys.blkSize;
        string filename = app.fileout + "_np" + NumberToString(sys.my_rank) + ".bin";
        sol.writeSol2File(filename, nein, bsz*nfin, ndims, writeQflag);
    }
    
    if (sys.my_rank==0) {
        printf("\n\n\n*************************************************\n");
        printf("******** TOTAL SIMULATION TIME: %g s ********", (double) ((clock()-t)/CLOCKS_PER_SEC));
        printf("\n*************************************************\n");
    }
}

#endif
