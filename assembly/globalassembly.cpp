#ifndef __GLOBALASSEMBLY
#define __GLOBALASSEMBLY

// Written by: C. Nguyen & P. Fernandez

// TODO: elemtimes and schurtimes may not be accurate in the threaded version of the code due to data race issues.

void mapElement2LinearSystem(sysstruct &sys, elemstruct &elem, meshstruct &mesh, 
        appstruct &app, tempstruct &temp, double *Hg, double *Rg, double *Rr, Int ie)
{
    double *He = &elem.H[0];
    double *Re = &elem.Rh[0];
    double *Ro = &elem.Rhonly[0];
    
    Int *rp = &sys.ent2entStart[0];
    Int *cj = &sys.ent2ent[0];
    Int *ind = &temp.ind[0];
    Int *jnd = &temp.jnd[0];
    
    Int e = mesh.elementtype[ie];
    Int *cg2dg = &mesh.cg2dg[e][0];
    Int *dg2cg = &mesh.dg2cg[e][0];

    Int inc = 1, a, b, i, j, k, m1, n1, m2, n2, t1, t2, m, n;
//     nfe = ndims[2];
//     ne = ndims[5];
//     npf = ndims[10];
//     nch = ndims[23];
//     ncf = ndims[75];    
    Int ne  = mesh.ne;
    Int nfe = mesh.nfes[e];
    Int npf = mesh.npfs[e][0];
    Int ndf = mesh.ndf[e];    
    Int ncf = mesh.ncf[e];       
    Int nch = app.nch;    

    Int isEDGelement = mesh.isEDGelement[ie];

    Int hybrid = app.problem[0];
    if (hybrid == 0) {          // HDG

        Int bsz = nch*npf;
        Int sz0 = npf*nfe;
        Int sz1 = nch*npf*nfe;
        Int sz2 = nch*bsz;
        Int sz3 = bsz*bsz;
        Int sz4 = nch*sz1;
        Int sz5 = npf*nch*sz1;
        
        for (a = 0; a<nfe; a++) {
            //i = mesh.t2f[a*ne+ie];
            i = mesh.t2f[ie*nfe+a];
            
            for (m1=0; m1<npf; m1++)
                ind[m1] = mesh.elcon[ie*sz0+a*npf+m1] - i*npf;
            
            for (m1=0; m1<npf; m1++)
                for (n1=0; n1<nch; n1++) {
// // //                     #pragma omp atomic
                    Rg[i*bsz+m1*nch+n1] += Re[a*bsz+ind[m1]*nch+n1];
// // //                     #pragma omp atomic
                    Rr[i*bsz+m1*nch+n1] += Ro[a*bsz+ind[m1]*nch+n1];
                }

            for (b=0; b<nfe; b++) {
                //j = mesh.t2f[b*ne+ie];
                j = mesh.t2f[ie*nfe+b];

                for (k = rp[i]; k<rp[i+1]; k++)
                    if (cj[k] == j) break;

                if (k == rp[i+1]) {
                    //cout<<app.my_rank<<" . "<<ie<<" . "<<rp[i]
                    error("Error No. 1 in mapElement2LinearSystem.\n");
                }

                for (m2=0; m2<npf; m2++)
                    jnd[m2] = mesh.elcon[ie*sz0+b*npf+m2] - j*npf;

                for (m1=0; m1<npf; m1++)
                    for (n1=0; n1<nch; n1++)
                        for (m2=0; m2<npf; m2++)
                            for (n2=0; n2<nch; n2++) {
// // //                                 #pragma omp atomic
                                Hg[k*sz3 + m1*sz2 + n1*bsz + m2*nch + n2] += He[b*sz5+jnd[m1]*sz4+n1*sz1+a*bsz+ind[m2]*nch+n2];
                            }
            }
        }
    }
    else if (hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {         // EDG
        Int nn  = ncf;
        Int nnn = npf*nfe;
        Int bsz = nch;
        Int sz0 = bsz*bsz;
        Int sz1 = bsz*nn;
        Int sz2 = bsz*nn*bsz;

        Int *elc = &mesh.elcon[ie*npf*nfe];
        for (a=0; a<nn; a++) {
            i = elc[cg2dg[a]];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rg[i*bsz+m] += Re[a*bsz+m];
// // //                 #pragma omp atomic
                Rr[i*bsz+m] += Ro[a*bsz+m];
            }
            
            for (b=0; b<nn; b++) {
                j = elc[cg2dg[b]];
                for (k = rp[i]; k<rp[i+1]; k++)
                    if (cj[k] == j) break;

                if (k == rp[i+1]) {
                    error("Error No. 2 in mapElement2LinearSystem.\n");
                }

                for (m=0; m<bsz; m++)
                    for (n = 0; n < bsz; n++) {
// // //                         #pragma omp atomic
                        Hg[k*sz0 + m*bsz + n] += He[b * sz2 + m * sz1 + a * bsz + n];
                    }
            }
        }
    }
    else if (app.hybrid == 2 && isEDGelement == 0) {
        Int nn  = npf*nfe;
        Int bsz = nch;
        Int sz0 = bsz*bsz;
        Int sz1 = bsz*nn;
        Int sz2 = bsz*nn*bsz;

        Int *elc = &mesh.elcon[ie*nn];
        for (a=0; a<nn; a++) {
            i = elc[a];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rg[i*bsz + m] += Re[a*bsz+m];
// // //                 #pragma omp atomic
                Rr[i*bsz + m] += Ro[a*bsz+m];
            }

            for (b=0; b<nn; b++) {
                j = elc[b];
                for (k = rp[i]; k<rp[i+1]; k++)
                    if (cj[k] == j) break;

                if (k == rp[i+1]) {
                    error("Error No. 3 in mapElement2LinearSystem.\n");
                }

                for (m=0; m<bsz; m++)
                    for (n = 0; n < bsz; n++) {
// // //                         #pragma omp atomic
                        Hg[k*sz0 + m*bsz + n] += He[b * sz2 + m * sz1 + a * bsz + n];
                    }
            }
        }
    }
}

void mapElement2RHS(sysstruct &sys, elemstruct &elem, meshstruct &mesh, 
        appstruct &app, tempstruct &temp, double* Rg, double* Rr, Int ie)
{
    double *Re = &elem.Rh[0];
    double *Ro = &elem.Rhonly[0];
    
    Int *rp = &sys.ent2entStart[0];
    Int *cj = &sys.ent2ent[0];
    Int *ind = &temp.ind[0];
    Int *jnd = &temp.jnd[0];
    
    Int e = mesh.elementtype[ie];
    Int *cg2dg = &mesh.cg2dg[e][0];    

    Int a, b, i, j, k, m1, n1, m2, n2, t1, t2, m, n;
//     nfe = ndims[2];
//     ne = ndims[5];
//     npf = ndims[10];
//     nch = ndims[23];
//     ncf = ndims[75];
    Int ne  = mesh.ne;
    Int nfe = mesh.nfes[e];
    Int npf = mesh.npfs[e][0];
    Int ndf = mesh.ndf[e];    
    Int ncf = mesh.ncf[e];       
    Int nch = app.nch;        

    Int isEDGelement = mesh.isEDGelement[ie];

    Int hybrid = app.problem[0];
    if (hybrid == 0) {            // HDG
        Int bsz = nch*npf;
        Int sz0 = npf*nfe;

        for (a = 0; a<nfe; a++) {
            //i = mesh.t2f[a*ne+ie];
            i = mesh.t2f[ie*nfe+a];

            for (m1=0; m1<npf; m1++)
                ind[m1] = mesh.elcon[ie*sz0+a*npf+m1] - i * npf;

            for (m1=0; m1<npf; m1++) {
                for (n1 = 0; n1 < nch; n1++) {
// // //                     #pragma omp atomic
                    Rg[i*bsz + m1*nch + n1] += Re[a * bsz + ind[m1] * nch + n1];
// // //                     #pragma omp atomic
                    Rr[i*bsz + m1*nch + n1] += Ro[a * bsz + ind[m1] * nch + n1];
                }
            }
        }
    }
    else if (hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // EDG or IEDG/HEDG with EDG element
        Int nn  = ncf;
        Int bsz = nch;
        Int *elc = &mesh.elcon[ie*npf*nfe];

        for (a=0; a<nn; a++) {
            i = elc[cg2dg[a]];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rg[i*bsz + m] += Re[a*bsz+m];
// // //                 #pragma omp atomic
                Rr[i*bsz + m] += Ro[a*bsz+m];
            }
        }
    }
    else if (app.hybrid == 2 && isEDGelement == 0) {        // IEDG/HDG with non-EDG element
        Int nn  = npf*nfe;
        Int bsz = nch;

        Int *elc = &mesh.elcon[ie*nn];
        for (a=0; a<nn; a++) {
            i = elc[a];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rg[i*bsz + m] += Re[a*bsz+m];
// // //                 #pragma omp atomic
                Rr[i*bsz + m] += Ro[a*bsz+m];
            }
        }
    }
}

void mapElement2Rh(sysstruct &sys, elemstruct &elem, meshstruct &mesh, appstruct &app, 
        tempstruct &temp, double* Rh, Int ie)
{    
    double *Ro = &elem.Rhonly[0];
    
    Int *rp = &sys.ent2entStart[0];
    Int *cj = &sys.ent2ent[0];
    Int *ind = &temp.ind[0];
    Int *jnd = &temp.jnd[0];
    
    Int e = mesh.elementtype[ie];
    Int *cg2dg = &mesh.cg2dg[e][0];    

    Int a, b, i, j, k, m1, n1, m2, n2, t1, t2, m, n;
//     nfe = ndims[2];
//     ne = ndims[5];
//     npf = ndims[10];
//     nch = ndims[23];
//     ncf = ndims[75];
    Int ne  = mesh.ne;
    Int nfe = mesh.nfes[e];
    Int npf = mesh.npfs[e][0];
    Int ndf = mesh.ndf[e];    
    Int ncf = mesh.ncf[e];       
    Int nch = app.nch;    


    Int isEDGelement = mesh.isEDGelement[ie];

    Int hybrid = app.problem[0];
    if (hybrid == 0) {            // HDG
        Int bsz = nch*npf;
        Int sz0 = npf*nfe;

        for (a = 0; a<nfe; a++) {
            //i = mesh.t2f[a*ne+ie];
            i = mesh.t2f[ie*nfe+a];

            for (m1=0; m1<npf; m1++)
                ind[m1] = mesh.elcon[ie*sz0+a*npf+m1] - i * npf;

            for (m1=0; m1<npf; m1++) {
                for (n1 = 0; n1 < nch; n1++) {
// // //                     #pragma omp atomic
                    Rh[i*bsz + m1*nch + n1] += Ro[a * bsz + ind[m1] * nch + n1];
                }
            }
        }
    }
    else if (hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // EDG or IEDG/HEDG with EDG element
        Int nn  = ncf;
        Int bsz = nch;
        Int *elc = &mesh.elcon[ie*npf*nfe];

        for (a=0; a<nn; a++) {
            i = elc[cg2dg[a]];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rh[i*bsz + m] += Ro[a*bsz+m];
            }
        }
    }
    else if (app.hybrid == 2 && isEDGelement == 0) {        // IEDG/HDG with non-EDG element
        Int nn  = npf*nfe;
        Int bsz = nch;

        Int *elc = &mesh.elcon[ie*nn];
        for (a=0; a<nn; a++) {
            i = elc[a];

            for (m=0; m<bsz; m++) {
// // //                 #pragma omp atomic
                Rh[i*bsz + m] += Ro[a*bsz+m];
            }
        }
    }
}

void updateQuasiNewtonStructures(sysstruct &sys, elemstruct &elem, meshstruct &mesh, solstruct &sol, appstruct &app, tempstruct &temp, Int ie)
{
    Int inc = 1, i;
//     nfe = ndims[2];
//     npv = ndims[9];
//     npf = ndims[10];
//     ncu = ndims[20];
//     nch = ndims[23];
//     ncf = ndims[75];        
    Int e = mesh.elementtype[ie];
    Int nfe = mesh.nfes[e];
    Int npv = mesh.npes[e];
    Int npf = mesh.npfs[e][0];
    Int ndf = mesh.ndf[e];    
    Int ncf = mesh.ncf[e];       
    Int nch = app.nch;    
    Int ncu = app.ncu;    

    Int n1 = npv*ncu;
    Int n2 = n1*n1;
    Int n4;
    if (app.hybrid == 0 || app.hybrid == 2)
        n4 = nch*npf*nfe*npv*ncu;
    else if (app.hybrid == 1)
        n4 = nch*ncf*npv*ncu;

    if (sys.quasiNewtonAccuracy == 0) {         // Dinv and K are stored in single precision
        if (sys.schurImplementation == 0) {
            /* Store inv(D) into sol structure */
            for (i=0; i<n2; i++)
                sol.DinvFloat[ie*n2+i] = (float) elem.BD[i];

            /* Store pivots for LU factors of D into sol structure */
            for (i=0; i<n1; i++)
                sol.ipivD[ie*n1+i] = temp.ipiv[i];
        }
        else if (sys.schurImplementation == 1) {
            /* Store block LU factors of D into sol structure */
            for (i=0; i<n2; i++)
                sol.DinvFloat[ie*n2+i] = (float) elem.D_LU[i];
        }
        /* Store K into sol structure */
        for (i=0; i<n4; i++)
            sol.Kfloat[ie*n4+i] = (float) elem.GK[i];
    }
    else if (sys.quasiNewtonAccuracy == 1) {    // Dinv and K are stored in double precision
        if (sys.schurImplementation == 0) {
            /* Store inv(D) into sol structure */
            DCOPY(&n2,&elem.BD[0],&inc,&sol.Dinv[ie*n2],&inc);

            /* Store pivots for LU factors of D into sol structure */
            for (i=0; i<n1; i++)
                sol.ipivD[ie*n1+i] = temp.ipiv[i];
        }
        else if (sys.schurImplementation == 1) {
            /* Store block LU factors of D into sol structure */
            DCOPY(&n2,&elem.D_LU[0],&inc,&sol.Dinv[ie*n2],&inc);
        }

        /* Store K into sol structure */
        DCOPY(&n4,&elem.GK[0],&inc,&sol.K[ie*n4],&inc);
    }
}

void collapseRowsVector(double * vector, double * vector_tmp, meshstruct &mesh, appstruct &app, Int ie)
{
    Int inc = 1, i;
//     Int nfe = ndims[2];
//     Int npf = ndims[10];
//     Int nch = ndims[23];
//     Int ncf = ndims[75];
//     Int ndf = npf*nfe;
    Int e = mesh.elementtype[ie];
    Int nfe = mesh.nfes[e];
    Int npf = mesh.npfs[e][0];
    Int ndf = mesh.ndf[e];    
    Int ncf = mesh.ncf[e];       
    Int nch = app.nch;                
    Int *dg2cg = &mesh.dg2cg[e][0];        

    double one = 1.0;

    Int n1 = ncf*nch;
    for (i=0; i<n1; i++)
        vector_tmp[i] = 0.0;
//    memset(&vector_tmp[0], 0.0, sizeof(double) * n1);
    n1 = nch;
    Int i_vector = 0;
    for (i = 0; i < ndf; i++) {
        DAXPY(&n1, &one, &vector[i_vector], &inc, &vector_tmp[dg2cg[i]*nch], &inc);
        i_vector += n1;
    }
    n1 = ncf*nch;
    DCOPY(&n1, &vector_tmp[0], &inc, &vector[0], &inc);
}

void assembleLinearSystem(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps)
{
    Int inc = 1, i, e, npv, ndf, ncf, ncu, nch, n1, n2, n3, n5;
//     nfe = ndims[2];
//     ne = ndims[5];
//     npv = ndims[9];
//     npf = ndims[10];
//     ncu = ndims[20];
//     nch = ndims[23];
//     ncf = ndims[75];
//     Int e = mesh.elementtype[ie];
//     Int ne  = mesh.ne;
//     Int nfe = mesh.nfes[e];
//     Int npv = mesh.npes[e];
//     Int npf = mesh.npfs[e][0];
//     Int ndf = mesh.ndf[e];    
//     Int ncf = mesh.ncf[e];       
//     Int nch = app.nch;    
//     Int ncu = app.ncu;        
    double normRuElem, normRhonly;

    double schurtimes[6];
    double elemtimes[32];
    for (i = 0; i < 6; i++) {
        schurtimes[i] = 0.0;
    }
    for (i = 0; i < 32; i++) {
        elemtimes[i] = 0.0;
    }

    clock_t t;
    double elemtime = 0.0, schurtime = 0.0, maptime = 0.0;

//     Int n1 = npv*ncu;
//     Int n2 = n1*n1;
//     Int n3, n5 = nch*npf*nfe;
//     if (app.hybrid == 0 || app.hybrid == 2) {          // HDG, IEDG and HEDG
//         n3 = npv*ncu*nch*npf*nfe;
//     }
//     else if (app.hybrid == 1) {     // EDG
//         n3 = npv*ncu*nch*ncf;
//     }

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );
    std::fill( sys.Hg.begin(), sys.Hg.end(), 0.0 );
    
    long long noFlops = 0;
    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(t,normRuElem) reduction(+:rNorm,elemtime,schurtime,maptime)
        for (Int ie=0; ie<mesh.ne; ie++) {
            e = mesh.elementtype[ie];            
            //nfe = mesh.nfes[e];
            npv = mesh.npes[e];            
            ndf = mesh.ndf[e];    
            ncf = mesh.ncf[e];       
            nch = app.nch;    
            ncu = app.ncu;      
            n1 = npv*ncu;
            n2 = n1*n1;
            n5 = nch*ndf;
            if (app.hybrid == 0 || app.hybrid == 2) {          // HDG, IEDG and HEDG
                n3 = npv*ncu*nch*ndf;
            }
            else if (app.hybrid == 1) {     // EDG
                n3 = npv*ncu*nch*ncf;
            }
            
            /* Compute elemental matrices and vectors */
            t = clock();
            assembleElementMatrixVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie, &elemtimes[0]);
            elemtime += clock() - t;
            
            normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
            rNorm += normRuElem*normRuElem;

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);            
            
            /* Peform Schur complement */
            t = clock();
            schurElementMatrixVector(elems[this_thread], mesh, app, temps[this_thread], sys.schurImplementation, ie, &schurtimes[0], &noFlops);
            schurtime += clock() - t;            
                        
            /* Store inv(D)*Ru into sol structure */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&sol.DinvRu[ie*n1],&inc);

            /* Store inv(D)*F into sol structure */
            DCOPY(&n3,&elems[this_thread].F[0],&inc,&sol.DinvF[ie*n3],&inc);

            if (app.quasiNewton == 1) {
                updateQuasiNewtonStructures(sys, elems[this_thread], mesh, sol, app, temps[this_thread], ie);
            }
            
            /* Map elemental matrices and vectors to form the global linear system */
            t = clock();
            mapElement2LinearSystem(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Hg[0], &sys.Rg[0], &sys.r[0], ie);
            maptime += clock() - t;
        }
// // //     }
    elemtime /= sys.noThreads;
    schurtime /= sys.noThreads;
    maptime /= sys.noThreads;

    // Copy Hg to single precision, if necessary:
    if (sys.matvecPrecision == 0) {
        for (i = 0; i < sys.Hg.size(); i++)
            sys.Hg_sp[i] = (float) sys.Hg[i];
    }

    /* Compute the residual norm */
    n1 = sys.Rg.size();
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly * normRhonly;
    sol.rNorm = sqrt(sol.rNorm);
    
    // Report CPU time breakdown for matrix assembly and Schur complement:
    if (sys.print >= 2) {
        printf("\nBREAKDOWN OF ASSEMBLY TIMES:\n");
        printf("\n1. Assemble full system: %g ms\n", (elemtime/CLOCKS_PER_SEC)*1.0e3);
        for (i = 0; i < 3; i++) {
            printf("1.%d: %g ms\n", i, (elemtimes[i]/CLOCKS_PER_SEC)*1.0e3);
        }
        printf("\n2. Perform Schur complement: %g ms\n", (schurtime/CLOCKS_PER_SEC)*1.0e3);
        for (i = 0; i < 6; i++) {
            printf("2.%d: %g ms\n", i, (schurtimes[i]/CLOCKS_PER_SEC)*1.0e3);
        }
        printf("\n3. Local-to-global mapping: %g ms\n\n", (maptime/CLOCKS_PER_SEC)*1.0e3);
    }
}

/* This function is used when QuasiNewton = 1 to assemble the RHS when the matrix is reused */
void assembleRHS(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps)
{
    Int inc = 1, i, e, npv, ndf, ncu, nch, n1, n5;
//     ne = ndims[5];
//     npv = ndims[9];
//     ncu = ndims[20];
//     nch = ndims[23];
//     nfe = ndims[2];
//     npf = ndims[10];
//     ncf = ndims[75];
//     Int n1 = npv*ncu;
//     Int n5 = nch*npf*nfe;
    nch = app.nch;    
    ncu = app.ncu;      
    double normRuElem, normRhonly;

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );

    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(normRuElem) reduction(+:rNorm)
        for (Int ie=0; ie<mesh.ne; ie++) {
            e = mesh.elementtype[ie];                        
            npv = mesh.npes[e];            
            ndf = mesh.ndf[e];    
            n1 = npv*ncu;
            n5 = nch*ndf;
            
            /* compute Ru and Rh only */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);
            
            normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
            rNorm += normRuElem*normRuElem;

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            /* perform Schur complement to compute elemental vectors for UH */
            schurElementVector(sys, elems[this_thread], mesh, sol, app, temps[this_thread], ie);

            /* store inv(D)*Ru into sol structure */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&sol.DinvRu[ie*n1],&inc);

            /* map elemental vectors to form the global RHS */
            mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie);
        }
// // //     }

    /* Compute the residual norm */
    n1 = sys.Rg.size();
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly * normRhonly;
    sol.rNorm = sqrt(sol.rNorm);
}

void computeResidualNorm(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps)
{
    Int inc = 1, i, e, npv, ncu, nch, ne, ndf, n1, n5;
//     ne = ndims[5];
//     npv = ndims[9];
//     ncu = ndims[20];
//     nch = ndims[23];
//     nfe = ndims[2];
//     npf = ndims[10];
//     ncf = ndims[75];
//     Int n1 = npv*ncu;
//     Int n5 = nch*npf*nfe;    
    nch = app.nch;    
    ncu = app.ncu;          
    Int isEDGelement;

    double normRuElem, normRhonly;

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );

    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(normRuElem) reduction(+:rNorm)
        for (Int ie=0; ie<mesh.ne; ie++) {
            e = mesh.elementtype[ie];                        
            npv = mesh.npes[e];            
            ndf = mesh.ndf[e];    
            n1 = npv*ncu;
            n5 = nch*ndf;
            
            /* compute Ru and Rh only */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);

            normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
            rNorm += normRuElem*normRuElem;

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            isEDGelement = mesh.isEDGelement[ie];

            // Collapse rows in Rhonly (EDG, IEDG and HEDG only)
            if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {
                collapseRowsVector(&elems[this_thread].Rhonly[0], &temps[this_thread].Rh_tmp[0], mesh, app, ie);
            }

            /* map elemental vectors to form the global RHS */
            mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie);
        }
// // //     }
    
    /* Compute the residual norm */
    n1 = sys.Rg.size();
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly * normRhonly;
    sol.rNorm = sqrt(sol.rNorm);
}


void computeResidualSerial(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, double* U, double* R)
{
    Int inc = 1, i, e, npv, ncu, nc, nch, ne, nh, nfe, ndf;
//     ne = ndims[5];
//     nh = ndims[8];
//     npv = ndims[9];
//     nc  = ndims[19];
//     ncu = ndims[20];
//     nch = ndims[23];
//     nfe = ndims[2];
//     npf = ndims[10];
//     ncf = ndims[75];
    ne = mesh.ne;
    nh = mesh.ndofuh;
    nch = app.nch;    
    ncu = app.ncu;            
    nc  = app.nc;
    e = mesh.elementtype[0];                        
    npv = mesh.npes[e];                        
    Int bsz = sys.blkSize;    
    Int n1 = npv*ncu;
    Int n2 = npv*nc;
    Int n4 = nh*nch;
    Int n5 = nch*ndf;
    Int n6 = sys.numEntities*bsz;
    Int *ind = &temps[0].ind[0];
    
    double normRuElem, normRhonly;
    
    double *Ru = &R[0];
    double *Rh = &R[n1*ne];
    double *Uu = &U[0];
    double *Uh = &U[n1*ne];
    double *UDG = &sol.UDG[0];
    double *UH = &sol.UH[0];
    double *SH = &sol.SH[0];
    
    // Copy U solution to sol structure and compute Q
    DCOPY(&n4, Uh, &inc, UH, &inc);
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for
        for (Int ie = 0; ie < ne; ie++) {                        
            DCOPY(&n1,&Uu[ie*n1],&inc,&UDG[ie*n2],&inc);        

            /* Compute Q */
            if (app.flag_q==1) {
                getQ(elems[this_thread], mesh, master, app, sol, temps[this_thread], 
                     &UDG[ie*n2], &UDG[ie*n2+n1], &UH[0], &SH[ie*n2+n1], ie);
            }
        }
// // //     }

    // Initialize Rh residual
    for (i = 0; i < n4; i++)
        Rh[i] = 0.0;
    
    // Compute residual vector:
    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        this_thread = 0;
// // //         #pragma omp for private(normRuElem) reduction(+:rNorm)
        for (Int ie=0; ie<ne; ie++) {
            /* Compute Ru and Rh */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);
            
            /* Store Ru */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&Ru[ie*n1],&inc);        
            
            // Add contribution of Ru in element to full residual norm
            normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
            rNorm += normRuElem*normRuElem;        
            
            /* Copy Rh to Rhonly */
            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);
            
            /* Map elemental vectors to form the residual */
            mapElement2Rh(sys, elems[this_thread], mesh, app, temps[this_thread], &Rh[0], ie);        
        }
// // //     }
    
    // Add contribution of Rh to full residual norm
    normRhonly = DNRM2(&n6, &Rh[0], &inc);
    sol.rNorm = rNorm + normRhonly*normRhonly;
}

#ifdef  HAVE_MPI

void computeEnt2entWeight(sysstruct &sys)
{
    Int i, inc = 1;
    Int ri, rj, nri, pij = -1, pijW = -1, pijHg = -1, pijKg = -1;
    Int bsz = sys.blkSize;
    Int bsz2 = bsz*bsz;
    double * ent2entWeight = &sys.ent2entWeight[0];
    
    if (sys.preconditioner == 0) {
        // RAS preconditioner
        for (ri = 0; ri < sys.BJ_nrows; ri++) {
            nri = sys.ent2entStart[ri+1] - sys.ent2entStart[ri];
            for (i = 0; i < nri; i++) {
                pijW += 1;
                pijHg += 1;
                ent2entWeight[pijW] = DNRM2(&bsz2, &sys.Hg[bsz2*pijHg], &inc);
            }
        }
    }
    else if (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2) {
        // Subdomain-wise or entity-wise BJ preconditioner
        for (ri = 0; ri < sys.BJ_nrows; ri++) {
            nri = sys.ent2entStart[ri+1] - sys.ent2entStart[ri];
            for (i = 0; i < nri; i++) {
                rj = sys.ent2ent[pij];
                if (rj < sys.BJ_nrows) {
                    pijW += 1;
                    pijHg += 1;
                    ent2entWeight[pijW] = DNRM2(&bsz2, &sys.Hg[bsz2*pijHg], &inc);
                }
                else {
                    pijW += 1;
                    pijKg += 1;
                    ent2entWeight[pijW] = DNRM2(&bsz2, &sys.Kg[bsz2*pijKg], &inc);
                }
            }
        }
    }
}

void saveEnt2entWeight(string filename, sysstruct &sys)
{
    if (sys.my_rank == 0) {
        printf("Saving interaction between different entities in the mesh.\n");
    }

    ofstream out(filename.c_str(), ios::out | ios::binary);
    if (!out) {
        cout <<"Unable to open file" << filename << endl;
    }
    if (out) {
        out.write( reinterpret_cast<char*>( &sys.ent2entWeight[0] ), sizeof(double) * sys.ent2entWeight.size() );
    }
    out.close();

    // Terminate execution
    if (sys.my_rank == 0) {
        printf("Execution will be terminated.\n");
    }
    error("\n");
}

void assembleLinearSystemMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps)
{
    Int inc = 1, i, e, n, ri, pi, neighbor, npv, ncu, nch, ndf, ncf, ne, n1, n2, n3, n5;
    Int bsz2 = sys.blkSize*sys.blkSize;
//     nfe = ndims[2];
//     ne = ndims[5];
//     npv = ndims[9];
//     npf = ndims[10];
//     ncu = ndims[20];
//     nch = ndims[23];
//     ncf = ndims[75];

    double normRuElem, normRhonly;
    double schurtimes[6];
    double elemtimes[32];
    for (i = 0; i < 6; i++)
        schurtimes[i] = 0.0;
    for (i = 0; i < 32; i++)
        elemtimes[i] = 0.0;

    clock_t t;
    double elemtime = 0.0, schurtime = 0.0, maptime = 0.0, waittime = 0.0;
    
//     Int n1 = npv*ncu;
//     Int n2 = n1*n1;
//     Int n3, n5 = nch*npf*nfe;
//     if (app.hybrid == 0 || app.hybrid == 2) {          // HDG, IEDG and HEDG
//         n3 = npv*ncu*nch*npf*nfe;
//     }
//     else if (app.hybrid == 1) {     // EDG
//         n3 = npv*ncu*nch*ncf;
//     }

    Int nein = sys.elempartpts[0] + sys.elempartpts[1];
    Int ne0 = sys.elempartpts[0];    

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );
    std::fill( sys.Hg.begin(), sys.Hg.end(), 0.0 );
    std::fill( sys.Kg.begin(), sys.Kg.end(), 0.0 );

    double rNorm = 0.0;
    long long noFlops = 0;
    // interface elements + other neighboring elements
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(t,normRuElem) reduction(+:rNorm,elemtime,schurtime,maptime)
        for (Int ie=ne0; ie<mesh.ne; ie++) {
            /* compute elemental matrices and vectors */
            t = clock();
                        
            e = mesh.elementtype[ie];            
            //nfe = mesh.nfes[e];
            npv = mesh.npes[e];            
            ndf = mesh.ndf[e];    
            ncf = mesh.ncf[e];       
            nch = app.nch;    
            ncu = app.ncu;      
            n1 = npv*ncu;
            n2 = n1*n1;
            n5 = nch*ndf;
            if (app.hybrid == 0 || app.hybrid == 2) {          // HDG, IEDG and HEDG
                n3 = npv*ncu*nch*ndf;
            }
            else if (app.hybrid == 1) {     // EDG
                n3 = npv*ncu*nch*ncf;
            }
            
            assembleElementMatrixVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie, &elemtimes[0]);
            elemtime += clock() - t;
            
            if (ie < nein) {
                normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
                rNorm += normRuElem*normRuElem;
            }

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            /* peform Schur complement */
            t = clock();
            schurElementMatrixVector(elems[this_thread], mesh, app, temps[this_thread], sys.schurImplementation, ie, &schurtimes[0], &noFlops);
            schurtime += clock() - t;

            /* store inv(D)*Ru into sol structure */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&sol.DinvRu[ie*n1],&inc);

            /* Store inv(D)*F into sol structure */
            DCOPY(&n3,&elems[this_thread].F[0],&inc,&sol.DinvF[ie*n3],&inc);

            if (app.quasiNewton == 1)
                updateQuasiNewtonStructures(sys, elems[this_thread], mesh, sol, app, temps[this_thread], ie);

            /* Map elemental matrices and vectors to form the global linear system */
            t = clock();
            mapElement2LinearSystem(sys, elems[this_thread], mesh, app, temps[this_thread], 
                         &sys.Hg[0], &sys.Rg[0], &sys.r[0], ie);
//             if (sys.preconditioner == 0)
//                 mapElement2LinearSystem(sys, elems[this_thread], mesh, app, temps[this_thread], 
//                         &ndims[0], &sys.Hg[0], &sys.Rg[0], &sys.r[0], ie);
//             else if (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2)
//                 mapElement2LinearSystemMPI(sys, elems[this_thread], mesh, app, temps[this_thread], 
//                         &ndims[0], &sys.Hg[0], &sys.Kg[0], &sys.Rg[0], &sys.r[0], ie);
            maptime += clock() - t;
        }
// // //     }
    
    // Send / receive blocks of global system from other processors (for RAS preconditioner only)
     Int request_counter = 0;
     if (sys.preconditioner == 0) {
         /* Copy some portion of Hg to buffsendmat */
         for (i=0; i<sys.nmatsend; i++) {
             pi = sys.matsend[i];
             DCOPY(&bsz2, &sys.Hg[bsz2*pi], &inc, &sys.buffsendmat[bsz2*i], &inc);
         }

         /* non-blocking send */
         Int nsend, psend = 0;
         for (n=0; n<sys.nnbsd; n++) {
             neighbor = sys.nbsd[n];
             nsend = sys.matsendpts[n]*bsz2;
             if (nsend>0) {
                MPI_Isend(&sys.buffsendmat[psend], nsend, MPI_DOUBLE, neighbor, 0,
                       MPI_COMM_WORLD, &sys.requests[request_counter]);
                 psend += nsend;
                 request_counter += 1;
             }
         }
        
         /* non-blocking receive */
         Int nrecv, precv = 0;
         for (n=0; n<sys.nnbsd; n++) {
             neighbor = sys.nbsd[n];
             nrecv = sys.matrecvpts[n]*bsz2;
             if (nrecv>0) {
                MPI_Irecv(&sys.buffrecvmat[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                       MPI_COMM_WORLD, &sys.requests[request_counter]);
                 precv += nrecv;
                 request_counter += 1;
             }
         }
      }
    
    // interior elements
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        this_thread = 0;
// // //         #pragma omp for private(t,normRuElem) reduction(+:rNorm,elemtime,schurtime,maptime)
        for (Int ie=0; ie<ne0; ie++) {
            /* compute elemental matrices and vectors */
            t = clock();
            assembleElementMatrixVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie, &elemtimes[0]);
            elemtime += clock() - t;

            normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
            rNorm += normRuElem*normRuElem;

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            /* Peform Schur complement */
            t = clock();
            schurElementMatrixVector(elems[this_thread], mesh, app, temps[this_thread], sys.schurImplementation, ie, &schurtimes[0], &noFlops);
            schurtime += clock() - t;

            /* Store inv(D)*Ru into sol structure */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&sol.DinvRu[ie*n1],&inc);

            /* Store inv(D)*F into sol structure */
            DCOPY(&n3,&elems[this_thread].F[0],&inc,&sol.DinvF[ie*n3],&inc);

            if (app.quasiNewton == 1) {
                updateQuasiNewtonStructures(sys, elems[this_thread], mesh, sol, app, temps[this_thread], ie);
            }

            /* Map elemental matrices and vectors to form the global linear system */
            t = clock();
            mapElement2LinearSystem(sys, elems[this_thread], mesh, app, temps[this_thread], 
                    &sys.Hg[0], &sys.Rg[0], &sys.r[0], ie);            
//             if (sys.preconditioner == 0)
//                 mapElement2LinearSystem(sys, elems[this_thread], mesh, app, temps[this_thread], 
//                         &ndims[0], &sys.Hg[0], &sys.Rg[0], &sys.r[0], ie);
//             else if (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2)
//                 mapElement2LinearSystemMPI(sys, elems[this_thread], mesh, app, temps[this_thread], 
//                         &ndims[0], &sys.Hg[0], &sys.Kg[0], &sys.Rg[0], &sys.r[0], ie);
            maptime += clock() - t;
        }
// // //     }
    elemtime /= sys.noThreads;
    schurtime /= sys.noThreads;
    maptime /= sys.noThreads;

    if (sys.my_rank == 0 && sys.print >= 2) {
        printf("No. flops for last operation in Schur complement: %lli (GFLOPs: %g)\n", noFlops, ((double) noFlops) / (1.0e9*((double) schurtimes[3]/CLOCKS_PER_SEC)));
    }
    
    /* Compute the residual norm */
    n1 = sys.BJ_nrows*sys.blkSize;
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly*normRhonly;
    
    // Copy Kg to single precision, if necessary:
    if (sys.matvecPrecision == 0 && (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2)) {
        for (i = 0; i < sys.Kg.size(); i++)
            sys.Kg_sp[i] = (float) sys.Kg[i];
    }

    // Wait until blocks of global system from other processors are received (for RAS preconditioner only)
    t = clock();
     if (sys.preconditioner == 0) {
         /* Wait until all send and receive operations are completely done */
         MPI_Waitall(request_counter, sys.requests, sys.statuses);

         /* Copy buffrecvmat to Hg */
         for (i=0; i<sys.nmatrecv; i++) {
             pi = sys.matrecv[i];
             DCOPY(&bsz2, &sys.buffrecvmat[bsz2*i], &inc, &sys.Hg[bsz2*pi], &inc);
         }
         waittime = clock() - t;
     }
     else {
         waittime = 0;
     }
    
    // Copy Hg to single precision, if necessary:
    if (sys.matvecPrecision == 0) {
        for (i = 0; i < sys.Hg.size(); i++)
            sys.Hg_sp[i] = (float) sys.Hg[i];
    }

    // Report CPU time breakdown for matrix assembly:
    if (sys.print >= 2) {
        if (sys.my_rank == 0) {
            printf("\nBREAKDOWN OF ASSEMBLY TIMES:\n");
            printf("\n1. Assemble full system: %g ms\n", (elemtime/CLOCKS_PER_SEC)*1.0e3);
            for (i = 0; i < 3; i++) {
                printf("1.%d: %g ms\n", i, (elemtimes[i]/CLOCKS_PER_SEC)*1.0e3);
            }
            printf("\n2. Perform Schur complement: %g ms\n", (schurtime/CLOCKS_PER_SEC)*1.0e3);
            for (i = 0; i < 6; i++) {
                printf("2.%d: %g ms\n", i, (schurtimes[i]/CLOCKS_PER_SEC)*1.0e3);
            }
            printf("\n3. Local-to-global mapping: %g ms\n", (maptime/CLOCKS_PER_SEC)*1.0e3);
            printf("\n4. Wait for matrix communication: %g ms\n\n", (waittime/CLOCKS_PER_SEC)*1.0e3);
        }
    }

    if (sys.computeGlobalEnt2entWeight == 1) {  // Compute global entity-to-entity weights and write them in a file
        // Compute local entity-to-entity weights (local = in current processor)
        computeEnt2entWeight(sys);

        // Compute global entity-to-entity weights and write them in a file (global = all processors)
//        string fileout = app.fileout + "_globalEnt2entWeight.bin";
//        saveGlobalEnt2entWeight(fileout, sys);
        string fileout = app.filein + "_ent2entWeight_np" + NumberToString(sys.my_rank) + ".bin";
        saveEnt2entWeight(fileout, sys);
    }
}

void assembleRHSMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps)
{
    Int inc = 1, i, e, npv, ncu, nch, ndf, n1, n5;
//     ne = ndims[5];
//     npv = ndims[9];
//     ncu = ndims[20];
//     nch = ndims[23];
//     nfe = ndims[2];
//     npf = ndims[10];
//     ncf = ndims[75];
    nch = app.nch;    
    ncu = app.ncu;              
//     Int n1 = npv*ncu;
//     Int n5 = nch*ndf;

    Int nein = sys.elempartpts[0] + sys.elempartpts[1];
    double normRuElem, normRhonly;

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );

    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(normRuElem) reduction(+:rNorm)
        for (Int ie=0; ie<mesh.ne; ie++) {
            e = mesh.elementtype[ie];                        
            npv = mesh.npes[e];            
            ndf = mesh.ndf[e];    
            n1 = npv*ncu;
            n5 = nch*ndf;
                        
            /* Compute Ru and Rh only */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);

            if (ie < nein) {
                normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
                rNorm += normRuElem*normRuElem;
            }

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            /* Perform Schur complement to compute elemental vectors for UH */
            schurElementVector(sys, elems[this_thread], mesh, sol, app, temps[this_thread], ie);

            /* Store inv(D)*Ru into sol structure */
            DCOPY(&n1,&elems[this_thread].Ru[0],&inc,&sol.DinvRu[ie*n1],&inc);
            
            /* Map elemental vectors to form the global RHS */
            mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie); 
//             if (sys.preconditioner == 0)
//                 mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &ndims[0], &sys.Rg[0], &sys.r[0], ie); 
//             else if (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2)            
//                 mapElement2RHSMPI(sys, elems[this_thread], mesh, app, temps[this_thread], &ndims[0], &sys.Rg[0], &sys.r[0], ie);       
        }
// // //     }

    /* Compute the residual norm */
    n1 = sys.BJ_nrows*sys.blkSize;
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly*normRhonly;
}

void computeResidualNormMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims)
{
    Int inc = 1, i, npv, ncu, nch, ne, nfe, npf, ncf;
    ne = ndims[5];
    npv = ndims[9];
    ncu = ndims[20];
    nch = ndims[23];
    nfe = ndims[2];
    npf = ndims[10];
    ncf = ndims[75];

    Int n1 = npv*ncu;
    Int n5 = nch*npf*nfe;
    Int isEDGelement;
    
    Int nein = sys.elempartpts[0] + sys.elempartpts[1];
    double normRuElem, normRhonly;

    std::fill( sys.r.begin(), sys.r.end(), 0.0 );
    std::fill( sys.Rg.begin(), sys.Rg.end(), 0.0 );

    double elemtimes[32];    
    for (i = 0; i < 32; i++)
        elemtimes[i] = 0.0;
    
    double rNorm = 0.0;
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for private(normRuElem) reduction(+:rNorm)
        for (Int ie=0; ie<ne; ie++) {
            /* Compute Ru and Rh only */
            assembleElementVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie);
            //assembleElementMatrixVector(elems[this_thread], mesh, master, app, sol, temps[this_thread], ie, &elemtimes[0]);
            
            if (ie < nein) {
                normRuElem = DNRM2(&n1, &elems[this_thread].Ru[0], &inc);
                rNorm += normRuElem*normRuElem;
            }

            DCOPY(&n5,&elems[this_thread].Rh[0],&inc,&elems[this_thread].Rhonly[0],&inc);

            isEDGelement = mesh.isEDGelement[ie];

            // Collapse rows in Rhonly (EDG, IEDG and HEDG only)
            if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {
                collapseRowsVector(&elems[this_thread].Rhonly[0], &temps[this_thread].Rh_tmp[0], mesh, app, ie);
            }
            
            /* Map elemental vectors to form the global RHS */
            if (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2)            
                //mapElement2RHSMPI(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie);
                mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie);
            else if (sys.preconditioner == 0)
                mapElement2RHS(sys, elems[this_thread], mesh, app, temps[this_thread], &sys.Rg[0], &sys.r[0], ie);
        }
// // //     }

    /* Compute the residual norm */
    n1 = sys.BJ_nrows*sys.blkSize;
    normRhonly = DNRM2(&n1, &sys.r[0], &inc);
    sol.rNorm = rNorm + normRhonly*normRhonly;
}

#endif

#endif
