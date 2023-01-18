#ifndef __DGSOLVERS
#define __DGSOLVERS

//#include "pDLAK.cpp"
#include "DIRKcoeff.cpp"
#include "solutionupdate.cpp"
//#include "../application/FM/getAVfield.cpp"
//#include "../application/FM/nsNew/computeDynamicSmagConstant.cpp"
//#include "../utilities/computeForces.cpp"
//#include "../utilities/computeAvgSolution.cpp"

// Written by: C. Nguyen & P. Fernandez

// TODO: Implement reduced-basis approximation for the non-linear initial guess in the serial version of the code.

void computeOrdering(sysstruct &sys, Int preconditioner);

void computePreconditioner(sysstruct &sys, Int preconditioner);

void gmres(sysstruct &sys, Int *convFlag);

void recomputeJacobian(sysstruct &sys, appstruct &app)
{
    printf("\n\n*** The Jacobian matrix will be updated in the next iteration.***\n\n");
    
    app.reuseJacobian = 0;
    sys.linearSolvesClocks = 0;
}

void decideIfAdaptGMREStol(sysstruct &sys, appstruct app, Int iter)
{
    if (sys.adaptGMREStol == 1 && app.reuseJacobian == 0 && iter > 1) {
        sys.adaptGMREStol = 0;
    }
    else {
        sys.adaptGMREStol = 1;
    }
}

void solveLinearProblem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{
    clock_t t;
    double assemblyTime, precondTime, orderingTime, GMREStime;

    // Assembly linear system (or RHS only)
    t = clock();
    if (app.reuseJacobian == 0) {
        /* compute Jacobian matrix and RHS */
        app.reusePreconditioner = 0;
        printf("Assembling linear system...\n");
        assembleLinearSystem(sys, elems, mesh, master, sol, app, temps);
        assemblyTime = ((clock() - t)/CLOCKS_PER_SEC)*1.0e3;
        if (app.quasiNewton == 1) {
            sys.lastAssemblyAndPrecClocks = (long) (clock() - t);
        }
    }
    else {
        /* compute RHS */
        app.reusePreconditioner = 1;
        if (app.reuseResidual == 0) { /* TODO: Why this line? I think what we need is a flag to determine if the true residual is available, compute it if not, and then assemble the RHS */
        }
        printf("Assembling RHS vector...\n");
        assembleRHS(sys, elems, mesh, master, sol, app, temps);
        assemblyTime = ((clock() - t)/CLOCKS_PER_SEC)*1.0e3;
    }

    // Compute ordering (if necessary)
    t = clock();
    if (app.reuseOrdering == 0) {
        printf("Computing ordering...\n");
        computeOrdering(sys, sys.preconditioner);
        app.reuseOrdering = 1;
    }
    orderingTime = ((clock() - t)/CLOCKS_PER_SEC)*1.0e3;

    // Compute preconditioner (if necessary)
    t = clock();
    if (app.reusePreconditioner == 0) {
        printf("Computing preconditioner...\n");
        computePreconditioner(sys, sys.preconditioner);
    }
    precondTime = ((clock() - t)/CLOCKS_PER_SEC)*1.0e3;
    if (app.quasiNewton == 1 && app.reusePreconditioner == 0) {
        sys.lastAssemblyAndPrecClocks += (long) (clock() - t);
    }

    // Solve linear system:
    printf("Solving linear system...\n");
    t = clock();
    gmres(sys, convFlag);
    GMREStime = ((clock() - t)/CLOCKS_PER_SEC)*1.0e3;
    if (app.quasiNewton == 1) {
        sys.linearSolvesClocks += (long) (clock() - t);
    }

    if (sys.print >= 1) {
        printf("\n\nBREAKDOWN OF CPU TIMES:\n");
        printf("1. Matrix assembly: %g ms\n", assemblyTime);
        printf("2. Compute preconditioner: %g ms\n", precondTime);
        printf("3. Compute ordering: %g ms\n", orderingTime);
        printf("4. GMRES: %g ms\n\n", GMREStime);
    }
}

void solveNonlinearProblem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{
    Int j;
    Int convLinearFlag[1];
    double alphaMin = 1.0e-4, rNorm = sys.NewtonTol + 1.0, alpha, oldNorm;
    
    /* Newton iteration */
    sys.robustMode = 0;
    Int iter = 0, trueNewtonIter = 0;
    while (rNorm > sys.NewtonTol && iter < sys.NewtonMaxiter && trueNewtonIter < sys.trueNewtonMaxiter) {
        iter += 1;
        printf("\nNewton iteration:  %d\n", iter);
        if (app.reuseJacobian == 0)
            trueNewtonIter += 1;

        if (sys.linearProblem == 0 && sys.adaptiveGMREStol == 1)
            decideIfAdaptGMREStol(sys, app, iter);

        /* Solve the linearized problem */
        sys.rNorm = rNorm;
        solveLinearProblem(sys, elems, mesh, master, sol, app, temps, ndims, convLinearFlag);
        oldNorm = sol.rNorm;
        if (isnan(oldNorm) || isinf(oldNorm)) {
            error("\nResidual norm in current Newton's iterate is NaN or Inf.\n");
        }
        
        if (*convLinearFlag == -1)      // Failure (GMRES breakdown or orthogonalization did not converge)
            break;
        
        // Decide whether or not to reuse the Jacobian (only applies for quasi-Newton):
        if (app.quasiNewton == 1) {
            if (sys.linearSolvesClocks < 3*sys.lastAssemblyAndPrecClocks)
                app.reuseJacobian = 1;
            else
                recomputeJacobian(sys, app);
        }
        
        // Check if initial guess is below Newton tolerance:
        if (iter == 1 && oldNorm < sys.NewtonTol) {
            rNorm = oldNorm;
            sys.rNorm = rNorm;
            printf("Old residual: %g,   New residual: %g    \n", oldNorm, rNorm);
            break;
        }
        
        /* Perform line search */
        sol.alpha = 1.0;
        while (1) {                        
            /* update solution */
            updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha);      /* This approach to update the solution saves storage but increases flop count */
            
            /* Compute the residual norm */
            computeResidualNorm(sys, elems, mesh, master, sol, app, temps);
            rNorm = sol.rNorm;
            if (isnan(rNorm) || isinf(rNorm)) {
                printf("\n\n***********************\n");
                printf("Residual norm in line search of Newton's method is NaN or Inf.\n");
                printf("\n***********************\n\n");
            }
            
            if (isnan(rNorm) || isinf(rNorm) || (rNorm > oldNorm && rNorm > sys.NewtonTol)) {
                alpha = fabs(sol.alpha/2);
                sol.alpha = -alpha;
                printf("Damped Newton coefficient: alpha = %g || Residual = %g || Old residual = %g\n", alpha, rNorm, oldNorm);
                if (app.quasiNewton == 1)
                    recomputeJacobian(sys, app);
                if (alpha < alphaMin)  {
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
        
        printf("Old residual: %g,   New residual: %g    \n", oldNorm, rNorm);                   
    }
    
    if (rNorm < sys.NewtonTol)
        *convFlag = 1;
    else {
        *convFlag = 0;
        if (app.quasiNewton == 1)
            recomputeJacobian(sys, app);
        printf("Newton's method did not converge to tolerance %g in %d iterations (%d true Newton iter).\n", sys.NewtonTol, iter, trueNewtonIter);
    }
}

void solveSteadyProblem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, Int* convFlag)
{       
    if (app.linearProblem == 0)
        solveNonlinearProblem(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);    
    else if (app.linearProblem == 1) {
        printf("Solve linear problem \n");
        solveLinearProblem(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
        /* TODO: Need to update solution (updateSol) and decide if the Jacobian is reused (the latter only applies for quasi-Newton) */
        sol.alpha = 1.0;
        updateSol(sys, elems, mesh, master, sol, app, temps, ndims, &sol.UDG[0], &sol.UH[0], &sys.x[0], sol.alpha); 
    }    
}

void solveUnsteadyProblem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
    /* TODO: We may want to save the solution at the end of some time-steps */
    Int i, j, k, m, istage, ie, numSolvesAttemped;
    Int ne = ndims[5];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int npv = ndims[9];
    Int numEntities = sys.numEntities;
    Int bsz = sys.blkSize;
    Int dirkStage = app.dirkStage;
    Int dirkOrder = app.dirkOrder;
    Int convFlag[1];
    
    clock_t t1, t2, tNewton;
    
    /* Get DIRK coefficients */
    double *dirkc = new double[dirkStage];
    double *dirkd = new double[dirkStage*dirkStage];
    double *dirkt = new double[dirkStage];
    DIRKcoeff(dirkc, dirkd, dirkt, app.dirkStage, dirkOrder);
    
    /* Initialize the previous timestep solution */
    sol.UDGi.resize(sol.UDG.size());
    sol.UHi.resize(sol.UH.size());
    if (app.wave==1) {
        sol.PDGi.resize(sol.PDG.size());
        sol.SP.resize(sol.PDG.size());
    }
    
    /* Initial guess for nonlinear problem. TODO: Implement the algorithm in the serial version of the code */
    Int sza = 3;
    Int szu = sol.UDG.size();
    Int szh = sol.UH.size();
    vector<double> UDGpre(szu*sza);
    vector<double> UHpre(szh*sza);
    vector<double> a(sza);
    sys.Ru_MR.resize(ne*npv*ncu);
    sys.dRuda_MR.resize(sza*ne*npv*ncu);
    sys.a2_MR.resize(sza);
    sys.a2s_MR.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.a2s_MR[i].resize(sys.a2_MR.size());
    sys.da_MR.resize(sza);
    sys.C_MR.resize(sza*sza);
    sys.Clocal_MR.resize(sza*sza);
    
    Int saveSolFreq = 700;                       // Every how many time-steps the solution is saved to a file
    Int recomputeOrderingFreq = 20;        // Every how many time-steps the ordering is recomputed
    Int writeQflag = 1;                         // 0: Only u is written to a binary file. 1: Both u and q are written to a binary file
    Int writeAvgSolution = 1;                   // 0: The time average solution is NOT saved to a file. 1: The time average solution is saved to a file
    Int timeStepToStartAvg = 1000;              // Time step to start averaging the solution
    double time = 0.0, dt;
    string timeStepsFileName = app.filein + "_timeSteps.out";
    for (i=0; i<app.dt.size(); i++) {
      t1 = clock();
        
      printf("Timestep :  %d ||  Time :   %g ||  dt :   %g \n", i+1, time, app.dt[i]);
        
      /* get timestep size */
      dt = app.dt[i];
        
      /* get solution from the previous timestep */
      for (j = 0; j< sol.UDG.size(); j++)
	sol.UDGi[j] = sol.UDG[j];
      for (j = 0; j< sol.UH.size(); j++)
	sol.UHi[j] = sol.UH[j];        
      if (app.wave == 1)
	for (j = 0; j< sol.PDG.size(); j++)
                sol.PDGi[j] = sol.PDG[j];
        //error("Fix line 400 in dgsolvers.cpp\n");
        if (app.AVflag == 6 || app.AVflag == 8 || app.AVflag == 9 || app.AVflag == 10) {
            Int avgType = 2;
            //getAVfield(&sol.avField_p1CG[0], &sol.UDGi[0], &mesh.dgnodes[0], sys, mesh, master, app, temp, ndims, avgType, app.AVflag);
            /////////DG_2_p1CG(&sol.avField_p1CG[0], avField_DG, mesh, master, ndims);
        }
        if (app.SGSmodel == 4) {
            Int minimizationTarget = 0;
            Int projPorder = app.porder[0]-1;
            //computeDynamicSmagConstant(&sol.Cs[0], sys, mesh, master, app, &sol.UDGi[0], &mesh.hAvg[0], projPorder, &ndims[0], minimizationTarget);
        }

        if (i == 0 || (i % recomputeOrderingFreq) == 0) {
            app.reuseOrdering = 0;
        }
        
        /* Loop over each DIRK stage */
        for (istage = 0; istage < app.dirkStage; istage++) {
            t2 = clock();
            
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
            
            // TODO: Code serial minimal residual algorithm

            /* Solve the steady problem corresponding to the istage stage */
            *convFlag = 0;
            numSolvesAttemped = 0;
            while (*convFlag != 1) {
                solveSteadyProblem(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
                if (*convFlag != 1) {
//                     numSolvesAttemped ++;
//                     if (ii < sza || numSolvesAttemped >= 2) {
//                         app.DIRKnotConverged[ii-1] = 1;
                        
                        printf("\n\nWARNING: Unsteady solve at time step %d, DIRK stage %d did not converge.\n",i+1,istage+1);
                        error("\n");
                        
// //                         /////////// HACK START ///////////
// //                         for (j=0; j<szu; j++)
// //                             sol.UDG[j] = UDGpre[jstar*szu+j];
// //                         for (j=0; j<szh; j++)
// //                             sol.UH[j] = UHpre[jstar*szh+j];
// //                         
// //                         if (istage == app.dirkStage - 1)
// //                             app.dt[i+1] = app.dt[i+1] / 2;
// //                         else {
// //                             app.dt[i] = app.dt[i] / 2;
// //                             dt = app.dt[i];
// //                         }
// //                         break;
// //                         /////////// HACK END ///////////
//                     }
//                     // Try without minimal residual algorithm for initial guess:
//                     for (j=0; j<szu; j++)
//                         sol.UDG[j] = UDGpre[jstar*szu+j];
//                     for (j=0; j<szh; j++)
//                         sol.UH[j] = UHpre[jstar*szh+j];
                }
            }
                        
            if (app.wave)   {
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
                        error("DIRK scheme not implemented yet.");
                }
                for (j=0; j<ne; j++)
                    for (k=0; k<ncu; k++)
                        sol.PDG[j*ncu+k] = (dt/app.fc_u)*(sol.UDG[j*nc+k]+sol.SP[j*ncu+k]);                 
            }
            printf("\n\nTOTAL TIME TO SOLVE TIME STAGE NO. %d IN TIME-STEP NO. %d: %g ms\n\n", istage+1, i+1, ((clock() - t2)*1.0e3)/CLOCKS_PER_SEC);
            writeTimeStepSize2File(timeStepsFileName, i+1, istage+1, app.time, dt);
        }        
        time = time + dt;
        
        // Compute average solution:
        //if ((i+1 >= timeStepToStartAvg) && writeAvgSolution == 1)
        //    computeAvgSolution(sol, (i+1)-(timeStepToStartAvg-1));
        
        // Write solution to file:
        if (((i+1) % saveSolFreq) == 0) {
            string filename = app.fileout + "_t" + NumberToString(i+1) + ".bin";
            sol.writeSol2File(filename, ne, bsz*numEntities, ndims, writeQflag);
            
            //if ((i+1 >= timeStepToStartAvg) && (writeAvgSolution == 1)) {
            //    string filename_avg = app.fileout + "_avg.bin";
            //    sol.writeAvgSol2File(filename_avg, ne, bsz*numEntities, ndims, writeQflag);
            //}
        }
        
        printf("\n\nTOTAL TIME TO SOLVE TIME-STEP NO. %d: %g ms\n\n", i+1, ((clock() - t1)*1.0e3)/CLOCKS_PER_SEC);
    }
    
//     // Report time-steps in which nonlinear solver did not converge:
//     printf("\n\nLIST OF TIME-STEP IN WHICH NONLINEAR SOLVER DID NOT CONVERGE:\n");
//     for (i = 0; i < app.DIRKnotConverged.size(); i++) {
//         if (app.DIRKnotConverged[i] == 1)
//             printf("Time-step No. %d, DIRK stage No. %d.\n",(i/app.dirkStage)+1,(i%app.dirkStage)+1);
//     }
    
    delete[] dirkd;
    delete[] dirkt;
    delete[] dirkc;
}

void solveProblem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
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
        solveUnsteadyProblem(sys, elems, mesh, master, sol, app, temps, ndims);    
    else {
        solveSteadyProblem(sys, elems, mesh, master, sol, app, temps, ndims, convFlag);
        if (*convFlag != 1)
            printf("\n\nWARNING: Steady solve did not converge.\n\n");
        
        Int writeQflag = 1;         
        string filename = app.fileout + "sol.bin";
        sol.writeSol2File(filename, mesh.ne, app.nch*mesh.ndofuh, ndims, writeQflag);
    }

    printf("\n\n\n*************************************************\n");
    printf("******** TOTAL SIMULATION TIME: %g s ********", ((clock()-t)/CLOCKS_PER_SEC)*1.0e0);
    printf("\n*************************************************\n");
}

#endif
