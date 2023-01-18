#ifndef __NEWTONINIT
#define __NEWTONINIT

// Written by: C. Nguyen & P. Fernandez

void computeUResidualMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, 
        double* Ru,double* UDGpre, double* UHpre, double* a, Int sza)
{
    Int inc = 1, i, j, npv, ncu, nch, ne, nfe, npf, ncf;
    ne = ndims[5];
    npv = ndims[9];
    ncu = ndims[20];
    nch = ndims[23];
    nfe = ndims[2];
    npf = ndims[10];
    ncf = ndims[75];
    
    double one = 1.0, zero = 0.0;
    char chn = 'N';

    Int n1 = npv*ncu;
    Int n5 = nch*npf*nfe;
    Int nein = sys.elempartpts[0] + sys.elempartpts[1];
    
    Int szu = sol.UDG.size();
    Int szh = sol.UH.size();    
    Int szr = n1*nein;
    
    for (j=0; j<szu; j++) 
        sol.UDG[j] = a[0]*UDGpre[j];        
    for (i=1; i<sza; i++)
        for (j=0; j<szu; j++) 
            sol.UDG[j] += a[i]*UDGpre[i*szu+j];
//     DGEMV(&chn, &szu, &sza, &one, &UDGpre[0], &szu, &a[0],
//                   &inc, &zero, &sol.UDG[0], &inc);
    
    for (j=0; j<szh; j++) 
        sol.UH[j] = a[0]*UHpre[j];        
    for (i=1; i<sza; i++)
        for (j=0; j<szh; j++) 
            sol.UH[j] += a[i]*UHpre[i*szh+j];
//     DGEMV(&chn, &szh, &sza, &one, &UHpre[0], &szh, &a[0],
//               &inc, &zero, &sol.UH[0], &inc);
    
    
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for
        for (Int ie=0; ie<nein; ie++) {
            /* compute Ru only */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);        

            /* store *Ru  */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&Ru[ie*n1],&inc);
        }
// // //     }
}

void NewtonInitMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, 
        double* UDGpre, double* UHpre, double* a, Int sza)
{
    Int inc = 1, i, j, ii, jj, npv, ncu, nch, ne, nfe, npf, ncf, info;
    ne = ndims[5];
    npv = ndims[9];
    ncu = ndims[20];
    nch = ndims[23];
    nfe = ndims[2];
    npf = ndims[10];
    ncf = ndims[75];

    Int n1 = npv*ncu;
    Int n5 = nch*npf*nfe;    
    Int nein = sys.elempartpts[0] + sys.elempartpts[1];    
    
    Int szu = sol.UDG.size();
    Int szh = sol.UH.size();
    Int szr = n1*nein;
    Int nb = 1;
    Int *ipiv = &temps[0].ipiv[0];
    Int noThreads = sys.noThreads;
    
    char chn = 'N', cht = 'T';
    
    /* Ru, dRuda, a2, da, C, Clocal */
    double *Ru = &sys.Ru_MR[0];
    double *dRuda = &sys.dRuda_MR[0];
    double *a2 = &sys.a2_MR[0];
    vector<vector<double> >* a2s = &sys.a2s_MR;
    double *da = &sys.da_MR[0];
    double *C = &sys.C_MR[0];
    double *Clocal = &sys.Clocal_MR[0];

    Int jstar = 0;
    for (i=0; i<sza; i++) {
        if (a[i] == 1.0) {
            jstar = i;
            break;
        }
    }

    // Copy of initial guess for MR algorithm (i.e. the solution at the previous DIRK stage)
    for (j=0; j<szu; j++) 
        sol.UDG_initMR[j] = sol.UDG[j];        
    
    for (j=0; j<szh; j++) 
        sol.UH_initMR[j] = sol.UH[j];
    
    double one = 1.0, zero = 0.0, minusone = -1.0;
    double oldNorm, rNorm, alphaMin = 1e-2;
    double da_n;
    Int iter = 0, itermax = 1;
    Int restrict2hyperplane = 0;
    if (app.AVflag > 0)
        restrict2hyperplane = 1;
    while (iter < itermax) {
        /* Compute residual vector */
        computeUResidualMPI(sys, elems, mesh, master, sol, app, temps, ndims, 
             Ru, UDGpre, UHpre, a, sza);
        oldNorm = sqrt(DDOTMPI(szr, &Ru[0], inc, &Ru[0], inc, noThreads));
        if ((isnan(oldNorm) || isinf(oldNorm)) && sys.my_rank == 0) {
            printf("\n\n***********************\n");
            printf("Residual norm in current iterate of Minimal Residual algorithm is NaN or Inf.\n");
            printf("The solution at the previous DIRK stage will be used as the initial guess of the nonlinear problem");
            printf("\n***********************\n\n");
        }
        
        /* Compute jacobian matrix using finite differences */
        for (i=0; i<sza; i++) {
            for (j=0; j<sza; j++)
                a2[j] = a[j];
            a2[i] = a2[i] + 1.0e-6;
            
            computeUResidualMPI(sys, elems, mesh, master, sol, app, temps, ndims, 
             &dRuda[i*szr], UDGpre, UHpre, a2, sza);

            for (j=0; j<szr; j++)
               dRuda[i*szr+j] = (dRuda[i*szr+j]-Ru[j])/1e-6;        
        }
        
        /* RHS */
        DGEMVMPI(cht, szr, sza, one, dRuda, szr, Ru, inc, zero, da, inc, a2, a2s, noThreads);
        
        /* Matrix */
        DGEMMMPI(cht, chn, sza, sza, szr, one, dRuda, szr, dRuda, szr, zero, C, sza, Clocal, noThreads);
                
        /* Solve the linear system */
        DGESV(&sza, &nb, C, &sza, &ipiv[0], da, &sza, &info);
        
        /* Project update direction onto feasible hyperplane (if such flag is activated) */
        if (restrict2hyperplane == 1) {
            da_n = 0.0;
            for (j = 0; j < sza; j++)
                da_n += da[j]/sqrt(sza);
            for (j = 0; j < sza; j++)
                da[j] -= da_n/sqrt(sza);
        }

        for (j=0; j<sza; j++)
            a2[j] = a[j];        
        
        /* line search */
        sol.alpha = 1;
        while (fabs(sol.alpha) > alphaMin) {
            for (j=0; j<sza; j++)                
                a[j] = a[j] - sol.alpha*da[j];
            
            /* Compute the residual norm */            
            computeUResidualMPI(sys, elems, mesh, master, sol, app, temps, ndims, 
              Ru, UDGpre, UHpre, a, sza);            
            rNorm = sqrt(DDOTMPI(szr, &Ru[0], inc, &Ru[0], inc, noThreads));
            
            if ((isnan(rNorm) || isinf(rNorm)) && sys.my_rank == 0) {
                printf("\n\n***********************\n");
                printf("Residual norm in line search of Minimal Residual algorithm is NaN or Inf.\n");
                printf("The solution at the previous DIRK stage will be used as the initial guess of the nonlinear problem");
                printf("\n***********************\n\n");
            }

            if (rNorm > oldNorm) {
                sol.alpha = -fabs(sol.alpha/2.0);
            }
            else {
                break;
            }            
        }

        if (rNorm <= oldNorm) {
            if (sys.my_rank == 0 && sys.print >= 1)
                printf("Newton initialization iter No. %d.   Old residual: %g,   New residual: %g,   alpha: %g \n", iter+1, oldNorm, rNorm, fabs(sol.alpha));
        }
        else {
            // Make initial guess equal to (UDG, UH) in previous iteration
            for (jj=0; jj<szu; jj++)
                sol.UDG[jj] = sol.UDG_initMR[jj];

            for (jj=0; jj<szh; jj++)
                sol.UH[jj] = sol.UH_initMR[jj];

            break;
        }
        
        if (abs(rNorm) < sys.NewtonTol)
            break;
        
        iter += 1; 
    }
}

#endif
