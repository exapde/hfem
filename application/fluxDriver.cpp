#ifndef __FLUXDRIVER
#define __FLUXDRIVER

#include "FM/chainRule_f_udg.c"
#include "FM/getViscosity.c"
#include "FM/eulerNew/flux_eulerNEW.c"
#include "FM/nsNew/flux_nsNEW.c"
#include "FM/ransSAnew/flux_ransSANEW.c"
#include "FM/poi/flux_poisson.c"
#include "SM/LE_uq/flux_leuq.c"
#include "SM/ledisp/flux_ledisp.c"
#include "FM/mhd/flux_mhd.c"
#include "../utilities/UDG2udg.cpp"

// Written by: C. Nguyen & P. Fernandez

void chainRuleJacobianFlux(double * f_UDG, double * f_udg, double * pg, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd)
{
    /* This function should be used only if ALEflag == 2 or ALEflag == 3 */
    /* Note: There is little to nothing to gain by using BLAS instead of for loops in this function. */

    int i, j, k, l, m, n, g, len, sz2, sz3;

    Int ALEflag = app.ALEflag;
    Int flag_q = app.flag_q;

    if (ALEflag != 2 && ALEflag != 3) {
        printf("Invalid value of ALEflag = %d in chainRuleJacobianFlux function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get Gg */
    double * Gg = &pg[3 * numPoints * nd];


    /* Get gg and grad_g1 */
    double * gg;
    double * grad_g1;
    if (ALEflag == 3) {
        gg = &pg[(3 * nd + nd * nd) * numPoints];
        grad_g1 = &pg[(3 * nd + nd * nd + 2) * numPoints];
    }
    else {
        gg = new double[numPoints];
        for (i = 0; i < numPoints; i++) {
            gg[i] = 1.0;
        }
    }


    /* Compute Ginv */
    double * Ginv = new double[numPoints * nd * nd];
    double * Gg11, * Gg21, * Gg31, * Gg12, * Gg22, * Gg32, * Gg13, * Gg23, * Gg33;
    double * Ginv11, * Ginv21, * Ginv31, * Ginv12, * Ginv22, * Ginv32, * Ginv13, * Ginv23, * Ginv33;
    if (nd == 2)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg12 = &Gg[2*numPoints];
        Gg22 = &Gg[3*numPoints];
        Ginv11 = &Ginv[0];
        Ginv21 = &Ginv[1*numPoints];
        Ginv12 = &Ginv[2*numPoints];
        Ginv22 = &Ginv[3*numPoints];

        for (i=0; i<numPoints; i++) {
            Ginv11[i] =   Gg22[i] / gg[i];
            Ginv12[i] = - Gg12[i] / gg[i];
            Ginv21[i] = - Gg21[i] / gg[i];
            Ginv22[i] =   Gg11[i] / gg[i];
        }
    }
    else if (nd == 3)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg31 = &Gg[2*numPoints];
        Gg12 = &Gg[3*numPoints];
        Gg22 = &Gg[4*numPoints];
        Gg32 = &Gg[5*numPoints];
        Gg13 = &Gg[6*numPoints];
        Gg23 = &Gg[7*numPoints];
        Gg33 = &Gg[8*numPoints];

        Ginv11 = &Ginv[0];
        Ginv21 = &Ginv[1*numPoints];
        Ginv31 = &Ginv[2*numPoints];
        Ginv12 = &Ginv[3*numPoints];
        Ginv22 = &Ginv[4*numPoints];
        Ginv32 = &Ginv[5*numPoints];
        Ginv13 = &Ginv[6*numPoints];
        Ginv23 = &Ginv[7*numPoints];
        Ginv33 = &Ginv[8*numPoints];

        for (i=0; i<numPoints; i++) {
            Ginv11[i] = (Gg22[i] * Gg33[i] - Gg23[i] * Gg32[i]) / gg[i];
            Ginv12[i] = (Gg23[i] * Gg31[i] - Gg21[i] * Gg33[i]) / gg[i];
            Ginv13[i] = (Gg21[i] * Gg32[i] - Gg22[i] * Gg31[i]) / gg[i];
            Ginv21[i] = (Gg13[i] * Gg32[i] - Gg12[i] * Gg33[i]) / gg[i];
            Ginv22[i] = (Gg11[i] * Gg33[i] - Gg13[i] * Gg31[i]) / gg[i];
            Ginv23[i] = (Gg12[i] * Gg31[i] - Gg11[i] * Gg32[i]) / gg[i];
            Ginv31[i] = (Gg12[i] * Gg23[i] - Gg13[i] * Gg22[i]) / gg[i];
            Ginv32[i] = (Gg13[i] * Gg21[i] - Gg11[i] * Gg23[i]) / gg[i];
            Ginv33[i] = (Gg11[i] * Gg22[i] - Gg12[i] * Gg21[i]) / gg[i];
        }
    }

    sz2 = numPoints * ncu;
    sz3 = numPoints * ncu * nd;
    for (i = 0; i < ncu; i++)
        for (j = 0; j < nd; j++)
            for (k = 0; k < ncu; k++)
                for (g = 0; g < numPoints; g++) {
                    f_UDG[g + k * numPoints + j * sz2 + i * sz3] = f_udg[g + k * numPoints + j * sz2 + i * sz3] / gg[g];

                    if (ALEflag == 3) {
                        for (l = 0; l < nd; l++)
                            for (m = 0; m < nd; m++) {
                                f_UDG[g + k * numPoints + j * sz2 + i * sz3] += grad_g1[g + l * numPoints] * Ginv[g + l * numPoints + m * numPoints * nd] *
                                                                                f_udg[g + k * numPoints + j * sz2 + (ncu + i + m * ncu) * sz3];
                            }
                    }
                }


    if (flag_q == 1) {
        for (l = 0; l < nd; l++)
            for (i = 0; i < ncu; i++)
                for (j = 0; j < nd; j++)
                    for (k = 0; k < ncu; k++)
                        for (g = 0; g < numPoints; g++) {
                            f_UDG[g + k * numPoints + j * sz2 + (ncu + i + l * ncu) * sz3] = Ginv[g + l * numPoints + 0 * numPoints * nd] * f_udg[g + k * numPoints + j * sz2 + (ncu + i + 0 * ncu) * sz3] / gg[g];
                            for (m = 1; m < nd; m++) {
                                f_UDG[g + k * numPoints + j * sz2 + (ncu + i + l * ncu) * sz3] += Ginv[g + l * numPoints + m * numPoints * nd] * f_udg[g + k * numPoints + j * sz2 + (ncu + i + m * ncu) * sz3] / gg[g];
                            }
                        }
    }


    /* Deallocate dynamic memory */
    delete[] Ginv;
    if (ALEflag == 2) {
        delete[] gg;
    }
}



void rigidALEcorrectionFlux(double * f, double * f_UDG, double * pg, double * UDG, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd, int computeJacobian)
{
    int i, j, k, l, m, n, g, sz2, sz3;

    Int ALEflag = app.ALEflag;

    double * Vg = &pg[(2 * nd) * numPoints];

    double * gg;
    if (ALEflag == 3) {
        gg = &pg[(3 * nd + nd * nd) * numPoints];
    }
    else if (ALEflag == 1 || ALEflag == 2) {
        gg = new double[numPoints];
        for (i = 0; i < numPoints; i++) {
            gg[i] = 1.0;
        }
    }
    else {
        printf("Invalid value of ALEflag = %d in rigidALEcorrectionFlux function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }


    /* Compute rigid ALE correction */
    sz2 = numPoints * ncu;
    for (j = 0; j < nd; j++)
        for (i = 0; i < ncu; i++)
            for (k = 0; k < numPoints; k++)
                f[k + i * numPoints + j * sz2] -= UDG[k + i * numPoints] * Vg[k + j * numPoints] / gg[k];

    if (computeJacobian == 1) {
        sz3 = numPoints * ncu * nd;
        for (j = 0; j < nd; j++)
            for (i = 0; i < ncu; i++)
                for (k = 0; k < numPoints; k++)
                    f_UDG[k + i * numPoints + j * sz2 + i * sz3] -= Vg[k + j * numPoints] / gg[k];
    }


    /* Deallocate dynamic memory */
    if (ALEflag == 1 || ALEflag == 2) {
        delete[] gg;
    }
}



void nonRigidALEcorrectionFlux(double * f, double * f_UDG, double * pg, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd, int computeJacobian)
{
    /* Note: There is little to nothing to gain by using BLAS instead of for loops in this function. */

    int i, j, k, l, m, n, g, sz2, sz3;
    Int len, inc = 1;

    Int ALEflag = app.ALEflag;

    double * Gg = &pg[3 * numPoints * nd];

    double * gg;
    if (ALEflag == 3) {
        gg = &pg[(3 * nd + nd * nd) * numPoints];
    }
    else if (ALEflag == 2) {
        gg = new double[numPoints];
        for (i = 0; i < numPoints; i++) {
            gg[i] = 1.0;
        }
    }
    else {
        printf("Invalid value of ALEflag = %d in nonRigidALEcorrectionFlux function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }


    /* Compute GinvT */
    double * GinvT = new double[numPoints * nd * nd];
    double * Gg11, * Gg21, * Gg31, * Gg12, * Gg22, * Gg32, * Gg13, * Gg23, * Gg33;
    double * GinvT11, * GinvT21, * GinvT31, * GinvT12, * GinvT22, * GinvT32, * GinvT13, * GinvT23, * GinvT33;
    if (nd == 2)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg12 = &Gg[2*numPoints];
        Gg22 = &Gg[3*numPoints];
        GinvT11 = &GinvT[0];
        GinvT21 = &GinvT[1*numPoints];
        GinvT12 = &GinvT[2*numPoints];
        GinvT22 = &GinvT[3*numPoints];

        for (i=0; i<numPoints; i++) {
            GinvT11[i] =   Gg22[i] / gg[i];
            GinvT12[i] = - Gg21[i] / gg[i];
            GinvT21[i] = - Gg12[i] / gg[i];
            GinvT22[i] =   Gg11[i] / gg[i];
        }
    }
    else if (nd == 3)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[numPoints];
        Gg31 = &Gg[2*numPoints];
        Gg12 = &Gg[3*numPoints];
        Gg22 = &Gg[4*numPoints];
        Gg32 = &Gg[5*numPoints];
        Gg13 = &Gg[6*numPoints];
        Gg23 = &Gg[7*numPoints];
        Gg33 = &Gg[8*numPoints];

        GinvT11 = &GinvT[0];
        GinvT21 = &GinvT[1*numPoints];
        GinvT31 = &GinvT[2*numPoints];
        GinvT12 = &GinvT[3*numPoints];
        GinvT22 = &GinvT[4*numPoints];
        GinvT32 = &GinvT[5*numPoints];
        GinvT13 = &GinvT[6*numPoints];
        GinvT23 = &GinvT[7*numPoints];
        GinvT33 = &GinvT[8*numPoints];

        for (i=0; i<numPoints; i++) {
            GinvT11[i] = (Gg22[i] * Gg33[i] - Gg23[i] * Gg32[i]) / gg[i];
            GinvT12[i] = (Gg13[i] * Gg32[i] - Gg12[i] * Gg33[i]) / gg[i];
            GinvT13[i] = (Gg12[i] * Gg23[i] - Gg13[i] * Gg22[i]) / gg[i];
            GinvT21[i] = (Gg23[i] * Gg31[i] - Gg21[i] * Gg33[i]) / gg[i];
            GinvT22[i] = (Gg11[i] * Gg33[i] - Gg13[i] * Gg31[i]) / gg[i];
            GinvT23[i] = (Gg13[i] * Gg21[i] - Gg11[i] * Gg23[i]) / gg[i];
            GinvT31[i] = (Gg21[i] * Gg32[i] - Gg22[i] * Gg31[i]) / gg[i];
            GinvT32[i] = (Gg12[i] * Gg31[i] - Gg11[i] * Gg32[i]) / gg[i];
            GinvT33[i] = (Gg11[i] * Gg22[i] - Gg12[i] * Gg21[i]) / gg[i];
        }
    }


    /* Compute non-rigid correction */
    double * fCorr = new double[numPoints * ncu * nd];
    sz2 = numPoints * ncu;
    for (i=0; i<nd; i++)
        for (j=0; j<ncu; j++)
            for (k=0; k<numPoints; k++) {
                fCorr[k + j * numPoints + i * sz2] = gg[k] * f[k + j * numPoints + 0 * sz2] * GinvT[k + 0 * numPoints + i * numPoints * nd];
                for (m=1; m<nd; m++){
                    fCorr[k + j * numPoints + i * sz2] += gg[k] * f[k + j * numPoints + m * sz2] * GinvT[k + m * numPoints + i * numPoints * nd];
                }
            }

    len = numPoints * ncu * nd;
    DCOPY(&len, &fCorr[0], &inc, &f[0], &inc);

    double * fCorr_UDG;
    if (computeJacobian == 1) {
        sz3 = numPoints * ncu * nd;
        fCorr_UDG = new double[numPoints * ncu * nd * nc];
        for (n=0; n<nc; n++)
            for (i=0; i<nd; i++)
                for (j=0; j<ncu; j++)
                    for (k=0; k<numPoints; k++) {
                        fCorr_UDG[k + j * numPoints + i * sz2 + n * sz3] = gg[k] * f_UDG[k + j * numPoints + 0 * sz2 + n * sz3] * GinvT[k + 0 * numPoints + i * numPoints * nd];
                        for (m=1; m<nd; m++) {
                            fCorr_UDG[k + j * numPoints + i * sz2 + n * sz3] += gg[k] * f_UDG[k + j * numPoints + m * sz2 + n * sz3] * GinvT[k + m * numPoints + i * numPoints * nd];
                        }
                    }

        len = numPoints * ncu * nd * nc;
        DCOPY(&len, &fCorr_UDG[0], &inc, &f_UDG[0], &inc);
    }


    /* Deallocate dynamic memory */
    delete[] GinvT;
    if (ALEflag == 2) {
        delete[] gg;
    }
    delete[] fCorr;
    if (computeJacobian == 1) {
        delete[] fCorr_UDG;
    }
}


//fluxDriver(f, f_udg, pg, udgg, odgg, mesh, master, app, sol, temp, ie, 1, computeJacobian);
void fluxDriver(double * f,double * f_UDG, double * pg, double * UDG, double * ODG, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        Int ie, Int quadchoice, Int computeJacobian, Int numPoints)
{
    Int nfe, ncu, nch, ncq, nc, nco, nd, ncd;
    nd = master.nd;
    nfe = master.nfe;
    ncd = app.ncd;    
    nc = app.nc;
    ncu = app.ncu;
    ncq = app.ncq;
    nch = app.nch;
    nco = app.nco;

    Int ALEflag = app.ALEflag;

    double *udg;
    double *odg;
    double *f_udg;

    /* 1. Compute (u,q) from (U, gradX U) */
    if (ALEflag == 2 || ALEflag == 3) {
        udg = new double[numPoints * nc];
        odg = new double[numPoints * nco];
        UDG2udg(udg, UDG, pg, app, numPoints, ncu, ncq, nc, nd);
        UDG2udg(odg, ODG, pg, app, numPoints, ncu, ncq, nco, nd);
        if (computeJacobian == 1)
            f_udg = new double[numPoints * ncu * nd * nc];
    }
    else {
        udg = &UDG[0];
        odg = &ODG[0];
        f_udg = &f_UDG[0];
    }

    double *param = &app.physicsparam[0];
    double time = app.time;
    
    /* 2. Compute physical fluxes */
    switch (app.appname) {
        case 0:
            flux_eulerNEW(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 1:
            flux_nsNEW(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 2:
            flux_poisson(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;                
        case 3:
            flux_ransSANEW(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 4:
            flux_leuq(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 5:
            flux_ledisp(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 6:
            flux_mhd(f, f_udg, pg, udg, odg, mesh, master, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        default: {
            printf("Application not implemented (appname = %d)\n",app.appname);
            exit(-1);
        }
    }

    /* 3. Chain rule for the Jacobian from (u,q) to (U, gradX U) */
    if ((ALEflag == 2 || ALEflag == 3) && computeJacobian == 1) {
        chainRuleJacobianFlux(f_UDG, f_udg, pg, app, numPoints, ncu, ncq, nc, nd);
    }

    /* 4. Introduce rigid ALE correction */
    if (ALEflag == 1 || ALEflag == 2 || ALEflag == 3) {
        rigidALEcorrectionFlux(f, f_UDG, pg, UDG, app, numPoints, ncu, ncq, nc, nd, computeJacobian);
    }

    /* 5. Non-rigid ALE correction */
    if (ALEflag == 2 || ALEflag == 3) {
        nonRigidALEcorrectionFlux(f ,f_UDG, pg, app, numPoints, ncu, ncq, nc, nd, computeJacobian);
    }

    /* Deallocate dynamic memory */
    if (ALEflag == 2 || ALEflag == 3) {
        delete[] udg; delete[] odg;
        if (computeJacobian == 1) {
            delete[] f_udg;
        }
    }
}

#endif
