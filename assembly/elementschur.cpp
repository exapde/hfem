#ifndef __ELEMENTSCHUR
#define __ELEMENTSCHUR

// Written by: C. Nguyen & P. Fernandez

void schur_primal_matrix(elemstruct &elem, meshstruct &mesh, appstruct &app, tempstruct &temp, Int schurImplementation, Int ie, double* schurtimes, long long * noFlops)
{
    double *D = &elem.BD[0];        /* (npv*ncu)x(npv*ncu) */
    double *F = &elem.F[0];         /* (npv*ncu)x(nch*ndf) */
    double *K = &elem.GK[0];        /* (nch*ndf)x(npv*ncu) */
    double *H = &elem.H[0];         /* (nch*ndf)x(nch*ndf) */
    
    Int inc = 1, i, j, k, l, info;
//     Int nd = ndims[0];
//     Int nfe = ndims[2];
//     Int npf = ndims[10];
//     Int npv = ndims[9];
//     Int nc  = ndims[19];
//     Int ncu = ndims[20];
//     Int nch = ndims[23];
//     Int ncf = ndims[75];
//     Int ndf = npf*nfe;
//     Int na = ncu*npv;
//     Int nb;
//     Int *ipiv = &temp.ipiv[0];    
    Int e = mesh.elementtype[ie];
    Int nfe = mesh.nfes[e];
    Int npv = mesh.npes[e];
    Int ndf = mesh.ndf[e];
    Int ncf = mesh.ncf[e];    
    Int nch = app.nch;
    Int ncu = app.ncu;
    Int na = ncu*npv;
    Int nb;
    Int *ipiv = &temp.ipiv[0];    
    Int *dg2cg = &mesh.dg2cg[e][0];

    char chn = 'N', charu = 'U';
    double one = 1.0, minusone = -1.0, zero = 0.0;

    clock_t t;

    Int isEDGelement = mesh.isEDGelement[ie];

    t = clock();
    if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0)) {
        nb = ndf*nch;
    }
    else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // Collapse rows in K, columns in F, and columns and rows in H
        nb = ncf*nch;

        Int n1 = ncu*npv*ncf*nch;
        for (i=0; i<n1; i++)
            temp.K_tmp[i] = 0.0;
//        memset(&temp.K_tmp[0], 0.0, sizeof(double) * n1);
        Int i_K = 0;
        for (i = 0; i < ncu*npv; i++)
            for (j = 0; j < ndf; j++)
                for (k = 0; k < nch; k++) {
                    temp.K_tmp[k + dg2cg[j]*nch + i*ncf*nch] += K[i_K];
                    i_K++;
                }
        DCOPY(&n1, &temp.K_tmp[0], &inc, &K[0], &inc);

        n1 = ncf*nch*ncu*npv;
        for (i=0; i<n1; i++)
            temp.F_tmp[i] = 0.0;
//        memset(&temp.F_tmp[0], 0.0, sizeof(double)*n1);
        n1 = nch*ncu*npv;
        Int i_F = 0;
        for (i = 0; i < ndf; i++) {
            DAXPY(&n1, &one, &F[i_F], &inc, &temp.F_tmp[dg2cg[i]*n1], &inc);
            i_F += n1;
        }
        n1 = ncf*nch*ncu*npv;
        DCOPY(&n1, &temp.F_tmp[0], &inc, &F[0], &inc);

        n1 = ncf*nch*ncf*nch;
        for (i=0; i<n1; i++)
            temp.H_tmp[i] = 0.0;
//        memset(&temp.H_tmp[0], 0.0, sizeof(double)*n1);
        Int i_H = 0;
        for (i = 0; i < ndf; i++)
            for (j = 0; j < nch; j++)
                for (k = 0; k < ndf; k++)
                    for (l = 0; l < nch; l++) {
                        temp.H_tmp[l + dg2cg[k]*nch + j*nch*ncf + dg2cg[i]*nch*ncf*nch] += H[i_H];
                        i_H ++;
                    }
        DCOPY(&n1, &temp.H_tmp[0], &inc, &H[0], &inc);
    }
    schurtimes[0] += clock() - t;
    
    if (schurImplementation == 0) {
        // Factor D using LU decomposition and store it in D
        t = clock();
        na = npv*ncu;
        DGETRF(&na,&na,&D[0],&na,ipiv,&info);
        schurtimes[1] += clock() - t;

        // Compute inv(D)*F and store it in F
        t = clock();
        DGETRS(&chn,&na,&nb,&D[0],&na,ipiv,&F[0],&na,&info);
        schurtimes[2] += clock() - t;
    }
    else if (schurImplementation == 1) {
        double * DiF = &F[0];
        double * D_LU = &elem.D_LU[0];
        double * D_LU_tmp = &elem.D_LU_tmp[0];
        double * DiF_ij = &elem.DiF_ij[0];
        double * workDii = &elem.workDii[0];
        Int * pivDii = &elem.pivDii[0];

        Int lenWorkDii = npv;
        Int i, j, k, l, m, n, ii, jj, len, i_L, i_Lnew, i_U, i_Unew, i_Diag, i_DiagNew;
        Int blkSzD = npv, blkSzD2 = npv * npv, nBlkRowsD = ncu, nBlkColmsF = nch, LD_D = blkSzD * nBlkRowsD;
        Int colSzF;

        if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0)) {       // HDG
            colSzF = ndf;
        }
        else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // EDG
            colSzF = ncf;
        }

        t = clock();

        // Copy D into D_LU with desired storage format (block format with first L blocks and then U blocks. Uii block is stored at the end of the U row)
        Int i_D_LU = 0;
        for (i = 0; i < nBlkRowsD; i++)
            for (j = 0; j < i; j++)
                for (k = 0; k < blkSzD; k++) {
                    DCOPY(&blkSzD, &D[blkSzD2*nBlkRowsD*j + i*blkSzD + k*blkSzD*nBlkRowsD], &inc, &D_LU[i_D_LU], &inc);
                    i_D_LU += blkSzD;
                }
        for (ii = 0; ii < nBlkRowsD; ii++) {
            i = nBlkRowsD- ii - 1;
            for (j = i+1; j < nBlkRowsD; j++)
                for (k = 0; k < blkSzD; k++) {
                    DCOPY(&blkSzD, &D[blkSzD2*nBlkRowsD*j + i*blkSzD + k*blkSzD*nBlkRowsD], &inc, &D_LU[i_D_LU], &inc);
                    i_D_LU += blkSzD;
                }
            for (k = 0; k < blkSzD; k++) {
                DCOPY(&blkSzD, &D[blkSzD2*nBlkRowsD*i + i*blkSzD + k*blkSzD*nBlkRowsD], &inc, &D_LU[i_D_LU], &inc);
                i_D_LU += blkSzD;
            }
        }

        // In-house block LU factorization of D
        for (i = 0; i < nBlkRowsD; i++) {
            // Compute inverse of block diagonal
            ii = nBlkRowsD- i - 1;
            i_Diag = blkSzD2 * ((nBlkRowsD* (nBlkRowsD- 1)) / 2 + ((ii + 1) * (ii + 2)) / 2 - 1);
            DGETRF(&blkSzD, &blkSzD, &D_LU[i_Diag], &blkSzD, pivDii, &info);
            DGETRI(&blkSzD, &D_LU[i_Diag], &blkSzD, pivDii, workDii, &lenWorkDii, &info);

            // Compute L factors and compute / update U factors:
            for (j = i+1; j < nBlkRowsD; j++) {
                jj = nBlkRowsD- j - 1;

                // Compute L factors
                i_L = ((j * (j-1)) / 2 + i) * blkSzD2;
                DGEMM(&chn, &chn, &blkSzD, &blkSzD, &blkSzD, &one, &D_LU[i_L], &blkSzD, &D_LU[i_Diag], &blkSzD, &zero, &D_LU_tmp[0], &blkSzD);
                DCOPY(&blkSzD2, &D_LU_tmp[0], &inc, &D_LU[i_L], &inc);

                // Compute / update U factors
                i_U = blkSzD2 * ((nBlkRowsD* (nBlkRowsD- 1)) / 2 + (ii * (ii + 1)) / 2);
                i_Lnew = ((j * (j-1)) / 2 + i + 1) * blkSzD2;
                for (k = i+1; k < j; k++) {
                    DGEMM(&chn, &chn, &blkSzD, &blkSzD, &blkSzD, &minusone, &D_LU[i_L], &blkSzD, &D_LU[i_U], &blkSzD, &one, &D_LU[i_Lnew], &blkSzD);
                    i_Lnew += blkSzD2;
                    i_U += blkSzD2;
                }

                i_DiagNew = blkSzD2 * ((nBlkRowsD* (nBlkRowsD- 1)) / 2 + ((jj + 1) * (jj + 2)) / 2 - 1);
                DGEMM(&chn, &chn, &blkSzD, &blkSzD, &blkSzD, &minusone, &D_LU[i_L], &blkSzD, &D_LU[i_U], &blkSzD, &one, &D_LU[i_DiagNew], &blkSzD);
                i_U += blkSzD2;

                i_Unew = blkSzD2 * ((nBlkRowsD* (nBlkRowsD- 1)) / 2 + (jj * (jj + 1)) / 2);
                for (k = j+1; k < nBlkRowsD; k++) {
                    DGEMM(&chn, &chn, &blkSzD, &blkSzD, &blkSzD, &minusone, &D_LU[i_L], &blkSzD, &D_LU[i_U], &blkSzD, &one, &D_LU[i_Unew], &blkSzD);
                    i_Unew += blkSzD2;
                    i_U += blkSzD2;
                }
            }
        }
        schurtimes[1] += clock() - t;

        // In-house forward and backward solves to compute inv(D)*F from block LU factors of D
        t = clock();
        for (j = 0; j < nBlkColmsF; j++) {
            i_D_LU = 0;
            // Forward solve
            for (i = 0; i < nBlkRowsD; i++)
                for (k = 0; k < i; k++) {
                    DGEMM(&chn, &chn, &blkSzD, &colSzF, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &F[blkSzD*colSzF*nBlkRowsD*j+blkSzD*k], &LD_D, &one, &F[blkSzD*colSzF*nBlkRowsD*j+blkSzD*i], &LD_D);
                    i_D_LU += blkSzD2;
                }

            // Backward solve
            for (ii = 0; ii < nBlkRowsD; ii++) {
                i = nBlkRowsD- ii - 1;
                for (k = i+1; k < nBlkRowsD; k++) {
                    DGEMM(&chn, &chn, &blkSzD, &colSzF, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &F[blkSzD*colSzF*nBlkRowsD*j+blkSzD*k], &LD_D, &one, &F[blkSzD*colSzF*nBlkRowsD*j+blkSzD*i], &LD_D);
                    i_D_LU += blkSzD2;
                }
                DGEMM(&chn, &chn, &blkSzD, &colSzF, &blkSzD, &one, &D_LU[i_D_LU], &blkSzD, &F[blkSzD*colSzF*nBlkRowsD*j+blkSzD*i], &LD_D, &zero, &DiF_ij[0], &blkSzD);
                i_D_LU += blkSzD2;
                for (k = 0; k < colSzF; k++) {
                    DCOPY(&blkSzD, &DiF_ij[k*blkSzD], &inc, &DiF[blkSzD*colSzF*nBlkRowsD*j + blkSzD*i + blkSzD*nBlkRowsD*k], &inc);
                }
            }
        }
        schurtimes[2] += clock() - t;
    }

    /* H = H - K*inv(D)*F */
    t = clock();
    na = npv*ncu;
    DGEMM(&chn, &chn, &nb, &nb, &na, &minusone, K, &nb, F, &na, &one, H, &nb);
    schurtimes[3] += clock() - t;
    noFlops[0] += 2*nb*nb*na;
}

void schur_primal_vector(elemstruct &elem, meshstruct &mesh, appstruct &app, tempstruct &temp, Int schurImplementation, Int ie, double* schurtimes)
{
    /* NOTE: This routine assumes that Rq = 0*/
    double *Dinv = &elem.BD[0];
    double *K = &elem.GK[0];
    double *Ru = &elem.Ru[0];
    double *Rh = &elem.Rh[0];
    double *Rhonly = &elem.Rhonly[0];

    Int inc = 1, i, info;
//     Int nfe = ndims[2];
//     Int npf = ndims[10];
//     Int npv = ndims[9];
//     Int ncu = ndims[20];
//     Int nch = ndims[23];
//     Int ncf = ndims[75];
//     Int ndf = npf*nfe;
//     Int na = npv*ncu;
//     Int nb;
//     Int nc = 1;
//     Int *ipiv = &temp.ipiv[0];
    //Int *dg2cg = &mesh.dg2cg[0];

    Int e = mesh.elementtype[ie];
    Int nfe = mesh.nfes[e];
    Int npv = mesh.npes[e];
    Int ndf = mesh.ndf[e];
    Int ncf = mesh.ncf[e];    
    Int nch = app.nch;
    Int ncu = app.ncu;    
    Int na = npv*ncu;
    Int nb;
    Int nc = 1;
    Int *ipiv = &temp.ipiv[0];    
    Int *dg2cg = &mesh.dg2cg[e][0];
    
    char chn = 'N';
    double one = 1.0, minusone = -1.0, zero = 0.0;

    Int isEDGelement = mesh.isEDGelement[ie];

    if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0)) {
        nb = ndf*nch;
    }
    else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // Collapse rows in Rh and Rhonly
        nb = ncf*nch;

        Int n1 = ncf*nch;
        for (i=0; i<n1; i++)
            temp.Rh_tmp[i] = 0.0;
//        memset(&temp.Rh_tmp[0], 0.0, sizeof(double) * n1);
        n1 = nch;
        Int i_Rh = 0;
        for (i = 0; i < ndf; i++) {
            DAXPY(&n1, &one, &Rh[i_Rh], &inc, &temp.Rh_tmp[dg2cg[i]*nch], &inc);
            i_Rh += n1;
        }
        n1 = ncf*nch;
        DCOPY(&n1, &temp.Rh_tmp[0], &inc, &Rh[0], &inc);
        DCOPY(&n1, &temp.Rh_tmp[0], &inc, &Rhonly[0], &inc);
    }

    clock_t t;

    if (schurImplementation == 0) {
        // compute inv(D)*Ru and store it in Ru
        t = clock();
        DGETRS(&chn,&na,&nc,Dinv,&na,ipiv,Ru,&na,&info);
        schurtimes[4] += clock() - t;
    }
    else if (schurImplementation == 1) {
        double * Ru_i = &elem.Ru_i[0];
        double * D_LU = &elem.D_LU[0];
        Int blkSzD = npv, blkSzD2 = npv * npv, nBlkRowsD = ncu;
        Int i, k, ii;

        t = clock();

        Int i_D_LU = 0;
        // Forward solve
        for (i = 0; i < nBlkRowsD; i++)
            for (k = 0; k < i; k++) {
                DGEMV(&chn, &blkSzD, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &Ru[blkSzD*k], &inc, &one, &Ru[blkSzD*i], &inc);
                i_D_LU += blkSzD2;
            }

        // Backward solve
        for (ii = 0; ii < nBlkRowsD; ii++) {
            i = nBlkRowsD- ii - 1;
            for (k = i+1; k < nBlkRowsD; k++) {
                DGEMV(&chn, &blkSzD, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &Ru[blkSzD*k], &inc, &one, &Ru[blkSzD*i], &inc);
                i_D_LU += blkSzD2;
            }
            DGEMV(&chn, &blkSzD, &blkSzD, &one, &D_LU[i_D_LU], &blkSzD, &Ru[i*blkSzD], &inc, &zero, &Ru_i[0], &inc);
            i_D_LU += blkSzD2;
            DCOPY(&blkSzD, &Ru_i[0], &inc, &Ru[i*blkSzD], &inc);
        }
        schurtimes[4] += clock() - t;
    }

    // Rh <- Rh - K*inv(D)*Ru
    t = clock();
    DGEMV(&chn, &nb, &na, &minusone, K, &nb, Ru, &inc, &one, Rh, &inc);
    schurtimes[5] += clock() - t;
}

void schurElementMatrixVector(elemstruct &elem, meshstruct &mesh, appstruct &app, tempstruct &temp, Int schurImplementation, Int ie, double* schurtimes, long long * noFlops)
{
    // Test case with p=3 hexa: 3, 8 and 9 dominate equally. One order of magnitude below are 7 and then 2. Then, everything else

    if (app.adjoint==0) {
        // Note: schur_primal_vector vs. schur_primal_matrix is x10 (p=1) and x50 (p=6)
        schur_primal_matrix(elem, mesh, app, temp, schurImplementation, ie, &schurtimes[0], noFlops);
        schur_primal_vector(elem, mesh, app, temp, schurImplementation, ie, &schurtimes[0]);
    }
    else {
        error("Adjoint solver not implemented yet.");
    }
}

void schurElementVector(sysstruct &sys, elemstruct &elem, meshstruct &mesh, solstruct &sol, appstruct &app, tempstruct &temp, Int ie)
{
    double *Ru = &elem.Ru[0];
    double *Rh = &elem.Rh[0];
    double *Rhonly = &elem.Rhonly[0];

//     Int nfe = ndims[2];
//     Int npf = ndims[10];
//     Int npv = ndims[9];
//     Int ncu = ndims[20];
//     Int nch = ndims[23];
//     Int ncf = ndims[75];
//     Int ndf = npf*nfe;
//     Int na = npv*ncu;
    Int nb1, nb2;
    Int nc = 1;
    Int inc = 1, i, info;
    char chn = 'N';
    double one = 1.0, minusone = -1.0, zero = 0.0;

    double *Dinv;
    double *D_LU;
    double *K;
    //Int *ipiv = &sol.ipivD[ie*na];
    //Int *dg2cg = &mesh.dg2cg[0];

    Int e = mesh.elementtype[ie];
    Int nfe = mesh.nfes[e];
    Int npv = mesh.npes[e];
    Int ndf = mesh.ndf[e];
    Int ncf = mesh.ncf[e];    
    Int nch = app.nch;
    Int ncu = app.ncu;    
    Int na = npv*ncu;
    Int *ipiv = &sol.ipivD[ie*na];
    Int *dg2cg = &mesh.dg2cg[e][0];
    
    Int isEDGelement = mesh.isEDGelement[ie];

    if (app.hybrid == 0 || app.hybrid == 2)
        nb2 = ndf*ncu;
    else if (app.hybrid == 1)
        nb2 = ncf*ncu;
    
    if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0))
        nb1 = ndf*ncu;
    else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {     // Collapse rows in Rh and Rhonly
        nb1 = ncf*ncu;
        Int n1 = ncf*nch;
        for (i=0; i<n1; i++)
            temp.Rh_tmp[i] = 0.0;
//        memset(&temp.Rh_tmp[0], 0.0, sizeof(double) * n1);
        n1 = nch;
        Int i_Rh = 0;
        for (i = 0; i < ndf; i++) {
            DAXPY(&n1, &one, &Rh[i_Rh], &inc, &temp.Rh_tmp[dg2cg[i]*n1], &inc);
            i_Rh += n1;
        }
        n1 = ncf*nch;
        DCOPY(&n1, &temp.Rh_tmp[0], &inc, &Rh[0], &inc);
        DCOPY(&n1, &temp.Rh_tmp[0], &inc, &Rhonly[0], &inc);
    }

    if (sys.quasiNewtonAccuracy == 0) {
        Dinv = &elem.D_inv[0];
        D_LU = &elem.D_inv[0];
        for (Int i_ = 0; i_ < na*na; i_++) {
            Dinv[i_] = (double) sol.DinvFloat[ie*na*na+i_];
            D_LU[i_] = (double) sol.DinvFloat[ie*na*na+i_];
        }
        K = &elem.K[0];
        for (Int i_ = 0; i_ < na*nb2; i_++) {
            K[i_] = (double) sol.Kfloat[ie*na*nb2+i_];
        }
    }
    else if (sys.quasiNewtonAccuracy == 1) {
        Dinv = &sol.Dinv[ie*na*na];
        D_LU = &sol.Dinv[ie*na*na];
        K = &sol.K[ie*na*nb2];
    }

    if (sys.schurImplementation == 0) {
        /* compute inv(D)*Ru and store it in Ru */
        DGETRS(&chn,&na,&nc,Dinv,&na,ipiv,Ru,&na,&info);
    }
    else if (sys.schurImplementation == 1) {
        double * Ru_i = &elem.Ru_i[0];
        Int blkSzD = npv, blkSzD2 = npv * npv, nBlkRowsD = ncu;
        Int i, k, ii;

        Int i_D_LU = 0;
        // Forward solve
        for (i = 0; i < nBlkRowsD; i++)
            for (k = 0; k < i; k++) {
                DGEMV(&chn, &blkSzD, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &Ru[blkSzD*k], &inc, &one, &Ru[blkSzD*i], &inc);
                i_D_LU += blkSzD2;
            }

        // Backward solve
        for (ii = 0; ii < nBlkRowsD; ii++) {
            i = nBlkRowsD- ii - 1;
            for (k = i+1; k < nBlkRowsD; k++) {
                DGEMV(&chn, &blkSzD, &blkSzD, &minusone, &D_LU[i_D_LU], &blkSzD, &Ru[blkSzD*k], &inc, &one, &Ru[blkSzD*i], &inc);
                i_D_LU += blkSzD2;
            }
            DGEMV(&chn, &blkSzD, &blkSzD, &one, &D_LU[i_D_LU], &blkSzD, &Ru[i*blkSzD], &inc, &zero, &Ru_i[0], &inc);
            i_D_LU += blkSzD2;
            DCOPY(&blkSzD, &Ru_i[0], &inc, &Ru[i*blkSzD], &inc);
        }
    }

    /* Rh <- Rh - K*inv(D)*Ru */
    DGEMV(&chn, &nb1, &na, &minusone, K, &nb1, Ru, &inc, &one, Rh, &inc);
}

// void adjoint_u(double* ae, double* fe, double* dudg, double* dudg_duh, double* BD,
//         double* F, double* GK, double* H, double* Ru, double* Rh, Int* ndims, Int flag)
// {
//     Int inc = 1, i, j, m, n, k, is, ks, iA, iB, iC, na, nb, info;
//     char *chn = "N", *charu = "U";
//     double one = 1.0, zero = 0.0, minusone = -1.0, fac, tm;
//
//     /* Get dimensions */
//     Int npv, ncu, nc, ndf;
//     nc  = ndims[0];
//     ncu = ndims[1];
//     npv = ndims[2];
//     ndf = ndims[3];
//
//
//     double *Dt  = (double *) malloc( npv*ncu*npv*ncu * sizeof(double) );
//     double *Kt  = (double *) malloc( npv*ncu*ncu*ndf * sizeof(double) );
//     double *Ft  = (double *) malloc( ncu*ndf*npv*ncu * sizeof(double) );
//     double *Ht  = (double *) malloc( ncu*ndf*ncu*ndf * sizeof(double) );
//
//     Int sz[10];
//
//     // Taking transpose
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*npv;
//     for (i=0; i<ncu; i++)
//         for (j=0; j<npv; j++)
//             for (m=0; m<ncu; m++)
//                 for (k=0; k<npv; k++)
//                     Dt[k+npv*m+sz[0]*j+sz[1]*i] = BD[j+npv*i+sz[0]*k+sz[1]*m];
//
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*ncu;
//     sz[2] = ncu*ndf;
//     sz[3] = ncu*ndf*npv;
//     for (i=0; i<ndf; i++)
//         for (j=0; j<ncu; j++)
//             for (m=0; m<ncu; m++)
//                 for (k=0; k<npv; k++)
//                     Kt[k+npv*m+sz[0]*j+sz[1]*i] = GK[j+ncu*i+sz[2]*k+sz[3]*m];
//
//     sz[0] = ncu*ndf;
//     sz[1] = ncu*ndf*npv;
//     sz[2] = npv*ncu;
//     sz[3] = npv*ncu*ncu;
//     for (i=0; i<ncu; i++)
//         for (j=0; j<npv; j++)
//             for (m=0; m<ndf; m++)
//                 for (k=0; k<ncu; k++)
//                     Ft[k+ncu*m+sz[0]*j+sz[1]*i] = F[j+npv*i+sz[2]*k+sz[3]*m];
//
//
//     sz[0] = ncu*ndf;
//     sz[1] = ncu*ndf*ncu;
//     for (i=0; i<ndf; i++)
//         for (j=0; j<ncu; j++)
//             for (m=0; m<ndf; m++)
//                 for (k=0; k<ncu; k++)
//                     Ht[k+ncu*m+sz[0]*j+sz[1]*i] = H[j+ncu*i+sz[0]*k+sz[1]*m];
//
//
//     double *Rt  = (double *) malloc( npv*ncu*(1+ndf*ncu) * sizeof(double) );
//     for (i=0; i<npv*ncu; i++)
//         Rt[i] = Ru[i];
//
//     for (j=0; j<ndf*ncu; j++)
//         for (i=0; i<npv*ncu; i++)
//             Rt[(j+1)*npv*ncu+i] = -Kt[j*npv*ncu+i];
//
//     Int *ipiv =  (Int *) malloc( npv*ncu * sizeof(Int) );
//     na = npv*ncu;
//     nb = 1+ndf*ncu;
//     DGESV(&na, &nb, Dt, &na, ipiv, Rt, &na, &info);
//
//     double *du, *dlu;
//     du  = &Rt[0];
//     dlu = &Rt[npv*ncu];
//
//     na = npv*ncu;
//     DCOPY(&na,du,&inc,dudg,&inc);
//
//     na = npv*ncu*ncu*ndf;
//     DCOPY(&na,dlu,&inc,dudg_duh,&inc);
//
//     if (flag == 0) {
//         na = ncu*ndf;
//         DCOPY(&na,Rh,&inc,fe,&inc);
//
//         na = ncu*ndf*ncu*ndf;
//         DCOPY(&na,Ht,&inc,ae,&inc);
//
//         na = ncu*ndf;
//         nb = npv*ncu;
//         DGEMV(chn, &na, &nb, &minusone, Ft, &na, du,
//                     &inc, &one, fe, &inc);
//         DGEMM(chn, chn, &na, &na, &nb, &one, Ft, &na, dlu,
//                     &nb, &one, ae, &na);
//     }
//
//     free(Dt); free(Kt); free(Ft); free(Ht);
//     free(Rt); free(ipiv);
// }
//
// void adjoint_uq(double* ae, double* fe, double* dudg, double* dudg_duh, double* M,
//         double* C, double* E, double* BD, double* F, double* GK, double* H,
//         double* Rq, double* Ru, double* Rh, Int* ndims, Int flag)
// {
//     Int inc = 1, i, j, m, n, k, is, ks, iA, iB, iC, na, nb, info;
//     char *chn = "N", *charu = "U";
//     double one = 1.0, zero = 0.0, minusone = -1.0, fac, tm;
//
//     /* Get dimensions */
//     Int nd, npv, ncu, nc, ndf, nd1;
//     nd  = ndims[0];
//     nc  = ndims[1];
//     ncu = ndims[2];
//     npv = ndims[3];
//     ndf = ndims[4];
//     nd1 = nd+1;
//
//     double *Ct  = (double *) malloc( npv*npv*nd * sizeof(double) );
//     double *Et  = (double *) malloc( ndf*npv*nd * sizeof(double) );
//     double *Bt  = (double *) malloc( npv*ncu*nd*npv*ncu * sizeof(double) );
//     double *Dt  = (double *) malloc( npv*ncu*npv*ncu * sizeof(double) );
//     double *Gt  = (double *) malloc( npv*ncu*nd*ncu*ndf * sizeof(double) );
//     double *Kt  = (double *) malloc( npv*ncu*ncu*ndf * sizeof(double) );
//     double *Ft  = (double *) malloc( ncu*ndf*npv*ncu * sizeof(double) );
//     double *Ht  = (double *) malloc( ncu*ndf*ncu*ndf * sizeof(double) );
//
//     Int sz[10];
//     sz[0] = npv*npv;
//     sz[1] = npv*ndf;
//
//     // Taking transpose
//     for (i=0; i<nd; i++)
//         for (j=0; j<npv; j++)
//             for (k=0; k<npv; k++)
//                 Ct[k+npv*j+sz[0]*i] = C[j+npv*k+sz[0]*i];
//
//     for (i=0; i<nd; i++)
//         for (j=0; j<npv; j++)
//             for (k=0; k<ndf; k++)
//                 Et[k+ndf*j+sz[1]*i] = E[j+npv*k+sz[1]*i];
//
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*npv;
//     for (i=0; i<ncu; i++)
//         for (j=0; j<npv; j++)
//             for (m=0; m<ncu; m++)
//                 for (k=0; k<npv; k++)
//                     Dt[k+npv*m+sz[0]*j+sz[1]*i] = BD[j+npv*i+sz[0]*k+sz[1]*m];
//
//     double *B;
//     B = &BD[npv*ncu*npv*ncu];
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*npv;
//     sz[2] = npv*ncu*npv*ncu;
//     for (n=0; n<nd; n++)
//         for (i=0; i<ncu; i++)
//             for (j=0; j<npv; j++)
//                 for (m=0; m<ncu; m++)
//                     for (k=0; k<npv; k++)
//                         Bt[k+npv*m+sz[0]*j+sz[1]*i+sz[2]*n] = B[j+npv*i+sz[0]*k+sz[1]*m+sz[2]*n];
//
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*ncu;
//     sz[2] = ncu*ndf;
//     sz[3] = ncu*ndf*npv;
//     for (i=0; i<ndf; i++)
//         for (j=0; j<ncu; j++)
//             for (m=0; m<ncu; m++)
//                 for (k=0; k<npv; k++)
//                     Kt[k+npv*m+sz[0]*j+sz[1]*i] = GK[j+ncu*i+sz[2]*k+sz[3]*m];
//
//     double *G;
//     G = &GK[npv*ncu*ncu*ndf];
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*ncu;
//     sz[2] = npv*ncu*ncu*ndf;
//     sz[3] = ncu*ndf;
//     sz[4] = ncu*ndf*npv;
//     sz[5] = ncu*ndf*npv*ncu;
//     for (n=0; n<nd; n++)
//         for (i=0; i<ndf; i++)
//             for (j=0; j<ncu; j++)
//                 for (m=0; m<ncu; m++)
//                     for (k=0; k<npv; k++)
//                         Gt[k+npv*m+sz[0]*j+sz[1]*i+sz[2]*n] = G[j+ncu*i+sz[3]*k+sz[4]*m+sz[5]*n];
//
//
//     sz[0] = ncu*ndf;
//     sz[1] = ncu*ndf*npv;
//     sz[2] = npv*ncu;
//     sz[3] = npv*ncu*ncu;
//     for (i=0; i<ncu; i++)
//         for (j=0; j<npv; j++)
//             for (m=0; m<ndf; m++)
//                 for (k=0; k<ncu; k++)
//                     Ft[k+ncu*m+sz[0]*j+sz[1]*i] = F[j+npv*i+sz[2]*k+sz[3]*m];
//
//
//     sz[0] = ncu*ndf;
//     sz[1] = ncu*ndf*ncu;
//     for (i=0; i<ndf; i++)
//         for (j=0; j<ncu; j++)
//             for (m=0; m<ndf; m++)
//                 for (k=0; k<ncu; k++)
//                     Ht[k+ncu*m+sz[0]*j+sz[1]*i] = H[j+ncu*i+sz[0]*k+sz[1]*m];
//
//     sz[0] = ncu*nd;
//     sz[1] = ncu*npv*ncu*nd;
//     sz[2] = ncu*ncu*ndf*nd;
//
//     DPOTRF(charu,&npv,M,&npv,&info);
//
//     DPOTRS(charu,&npv,&sz[0],M,&npv,Rq,&npv,&info);
//
//     DPOTRS(charu,&npv,&sz[1],M,&npv,Bt,&npv,&info);
//
//     DPOTRS(charu,&npv,&sz[2],M,&npv,Gt,&npv,&info);
//
//     na = ncu*ncu*ndf;
//     nb = ncu*npv*ncu;
//     for (i=0; i<nd; i++) {
//         iA = i*npv*npv;
//         iB = i*npv*ncu;
//         DGEMM(chn, chn, &npv, &ncu, &npv, &one, &Ct[iA], &npv, &Rq[iB], &npv,
//             &one, Ru, &npv);
//
//         iB = i*npv*ncu*npv*ncu;
//         DGEMM(chn, chn, &npv, &nb, &npv, &one, &Ct[iA], &npv, &Bt[iB], &npv,
//             &one, Dt, &npv);
//
//         iB = i*npv*ncu*ncu*ndf;
//         DGEMM(chn, chn, &npv, &na, &npv, &one, &Ct[iA], &npv, &Gt[iB], &npv,
//             &one, Kt, &npv);
//     }
//
//     double *Rt  = (double *) malloc( npv*ncu*(1+ndf*ncu) * sizeof(double) );
//     for (i=0; i<npv*ncu; i++)
//         Rt[i] = Ru[i];
//
//     for (j=0; j<ndf*ncu; j++)
//         for (i=0; i<npv*ncu; i++)
//             Rt[(j+1)*npv*ncu+i] = -Kt[j*npv*ncu+i];
//
//     Int *ipiv =  (Int *) malloc( npv*ncu * sizeof(Int) );
//     na = npv*ncu;
//     nb = 1+ndf*ncu;
//     DGESV(&na, &nb, Dt, &na, ipiv, Rt, &na, &info);
//
//     double *du, *dlu;
//     du  = &Rt[0];
//     dlu = &Rt[npv*ncu];
//
//     double *dq  = (double *) malloc( npv*ncu*nd * sizeof(double) );
//     double *dlq  = (double *) malloc( npv*ncu*ncu*ndf*nd * sizeof(double) );
//     na = ncu*ndf;
//     nb = npv*ncu;
//     for (i=0; i<nd; i++) {
//         iA = i*npv*ncu*npv*ncu;
//         iB = i*npv*ncu;
//         DGEMV(chn, &nb, &nb, &minusone, &Bt[iA], &nb, du, &inc,
//             &zero, &dq[iB], &inc);
//
//         iB = i*npv*ncu*ncu*ndf;
//         DGEMM(chn, chn, &nb, &na, &nb, &minusone, &Bt[iA], &nb, dlu, &nb,
//             &zero, &dlq[iB], &nb);
//     }
//
//     na = npv*ncu*nd;
//     DAXPY(&na, &one, Rq, &inc, dq, &inc);
//     na = npv*ncu*ncu*ndf*nd;
//     DAXPY(&na, &minusone, Gt, &inc, dlq, &inc);
//
//     for (i=0; i<npv*ncu; i++)
//         dudg[i] = du[i];
//
//     for (i=0; i<npv*ncu*nd; i++)
//         dudg[npv*ncu+i] = dq[i];
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*ncu;
//     sz[2] = npv*ncu*(nd+1);
//     sz[3] = npv*ncu*(nd+1)*ncu;
//     for (j=0; j<ndf; j++)
//         for (k=0; k<ncu; k++)
//             for (m=0; m<ncu; m++)
//                 for (n=0; n<npv; n++)
//                     dudg_duh[n+npv*m+sz[2]*k+sz[3]*j] = dlu[n+npv*m+sz[0]*k+sz[1]*j];
//
//
//     sz[0] = npv*ncu;
//     sz[1] = npv*ncu*(nd+1);
//     sz[2] = npv*ncu*(nd+1)*ncu;
//     sz[3] = npv*ncu*ncu;
//     sz[4] = npv*ncu*ncu*ndf;
//     sz[5] = npv*ndf;
//     for (j=0; j<ndf; j++)
//         for (k=0; k<ncu; k++)
//             for (i=0; i<nd; i++)
//                 for (m=0; m<ncu; m++)
//                     for (n=0; n<npv; n++)
//                         dudg_duh[n+npv*m+sz[0]*(i+1)+sz[1]*k+sz[2]*j]
//                           = dlq[n+npv*m+sz[0]*k+sz[3]*j+sz[4]*i];
//
//
//     if (flag == 0) {
//         na = ncu*ndf;
//         DCOPY(&na,Rh,&inc,fe,&inc);
//
//         na = ncu*ndf*ncu*ndf;
//         DCOPY(&na,Ht,&inc,ae,&inc);
//
//         na = ncu*ndf;
//         nb = npv*ncu;
//         DGEMV(chn, &na, &nb, &minusone, Ft, &na, du,
//                         &inc, &one, fe, &inc);
//         DGEMM(chn, chn, &na, &na, &nb, &one, Ft, &na, dlu,
//                         &nb, &one, ae, &na);
//
//         double *fet  = (double *) malloc( ndf*ncu * sizeof(double) );
//         double *aet  = (double *) malloc( ndf*ncu*ncu*ndf * sizeof(double) );
//
//         Int colsB = ncu*ncu*ndf;
//         DGEMM(chn, chn, &ndf, &ncu, &npv, &one, &Et[0], &ndf, &dq[0],
//                     &npv, &zero, &fet[0], &ndf);
//         DGEMM(chn, chn, &ndf, &colsB, &npv, &one, &Et[0], &ndf, &dlq[0],
//                     &npv, &zero, &aet[0], &ndf);
//         for (i=1; i<nd; i++) {
//             iA = i*ndf*npv;
//             iB = i*npv*ncu;
//             DGEMM(chn, chn, &ndf, &ncu, &npv, &one, &Et[iA], &ndf, &dq[iB],
//                         &npv, &one, &fet[0], &ndf);
//
//             iB = i*npv*ncu*ncu*ndf;
//             DGEMM(chn, chn, &ndf, &colsB, &npv, &one, &Et[iA], &ndf, &dlq[iB],
//                         &npv, &one, &aet[0], &ndf);
//         }
//
//         for (j=0; j<ndf; j++)
//             for (i=0; i<ncu; i++)
//                 fe[i+ncu*j] = fe[i+ncu*j] - fet[j+ndf*i];
//
//         sz[0] = ncu*ndf;
//         sz[1] = ncu*ndf*ncu;
//         for (n=0; n<ndf; n++)
//             for (m=0; m<ncu; m++)
//                 for (j=0; j<ndf; j++)
//                     for (i=0; i<ncu; i++)
//                         ae[i+ncu*j+sz[0]*m+sz[1]*n] = ae[i+ncu*j+sz[0]*m+sz[1]*n] + aet[j+ndf*i+sz[0]*m+sz[1]*n];
//
//         free(fet); free(aet);
//     }
//
//     free(Et); free(Ct); free(Dt); free(Bt);
//     free(Gt); free(Kt); free(Ft); free(Ht);
//     free(Rt); free(ipiv); free(dq); free(dlq);
// }
//
// void adjoint_localsolve(double* dudg, double* dudg_duh,
//         double* M, double* C, double* E, double* BD, double* F, double* G, double* H,
//         double* Rq, double* Ru, double* Rh, Int* ndims)
// {
//     double *ae, *fe;
//     if (M[0] == 0.0)
//         adjoint_u(ae, fe, dudg, dudg_duh, BD, F, G, H, Ru, Rh, ndims, 1);
//     else
//         adjoint_uq(ae, fe, dudg, dudg_duh, M, C, E, BD, F, G, H, Rq, Ru, Rh, ndims, 1);
// }
//
// void primal_localsolve(double* dudg, double* dudg_duh, double* M, double* C, double* E,
//         double* BD, double* F, double* Rq, double* Ru, Int* ndims)
// {
//     if (M[0] == 0.0)
//         primal_u(dudg, dudg_duh, BD, F, Ru, ndims);
//     else
//         primal_uq(dudg, dudg_duh, M, C, E, BD, F, Rq, Ru, ndims);
// }
//
// void localsolve(double* dudg, double* dudg_duh,
//         double* M, double* C, double* E, double* BD, double* F, double* G, double* H,
//         double* Rq, double* Ru, double* Rh, Int adjoint, Int* ndims)
// {
//     if (adjoint==1)
//         adjoint_localsolve(dudg, dudg_duh, M, C, E, BD, F, G, H, Rq, Ru, Rh, ndims);
//     else
//         primal_localsolve(dudg, dudg_duh, M, C, E, BD, F, Rq, Ru, ndims);
// }
//
//
// void primal_schur(double* ae, double* fe, double* dudg, double* dudg_duh,
//         double* M, double* C, double* E, double* BD, double* F, double* G, double* H,
//         double* Rq, double* Ru, double* Rh, Int* ndims)
// {
//     primal_localsolve(dudg, dudg_duh, M, C, E, BD, F, Rq, Ru, ndims);
//
//     Int inc=1, na, nb;
//     na = npf*nfe*nch;
//     nb = na*na;
//     DCOPY(na,Rh,inc,fe,inc);
//     DCOPY(nb,H,inc,ae,inc);
//
//     Int rowsA = nch*npf*nfe;
//     Int colsA = npv*nc;
//     DGEMM(chn, chn, rowsA, rowsA, colsA, &one, G, rowsA, dudg_duh,
//            colsA, &one, ae, rowsA);
//     DGEMV(chn, rowsA, colsA, &minusone, G, rowsA, dudg,
//                     inc, &one, fe, inc);
// }
//
// void adjoint_schur(double* ae, double* fe, double* dudg, double* dudg_duh,
//         double* M, double* C, double* E, double* BD, double* F, double* G, double* H,
//         double* Rq, double* Ru, double* Rh, Int* ndims)
// {
//     if (M[0] == 0.0)
//         adjoint_u(ae, fe, dudg, dudg_duh, BD, F, G, H, Ru, Rh, ndims, 0);
//     else
//         adjoint_uq(ae, fe, dudg, dudg_duh, M, C, E, BD, F, G, H, Rq, Ru, Rh, ndims, 0);
// }
//
// void schur(double* ae, double* fe, double* dudg, double* dudg_duh,
//         double* M, double* C, double* E, double* BD, double* F, double* G, double* H,
//         double* Rq, double* Ru, double* Rh, Int adjoint, Int* ndims)
//
// {
//     if (adjoint==1)
//         adjoint_schur(ae, fe, dudg, dudg_duh, M, C, E, BD, F, G, H, Rq, Ru, Rh, ndims);
//     else
//         primal_schur(ae, fe, dudg, dudg_duh, M, C, E, BD, F, G, H, Rq, Ru, Rh, ndims);
// }
//

#endif
