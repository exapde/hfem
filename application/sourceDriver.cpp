#ifndef __SOURCEDRIVER
#define __SOURCEDRIVER

#include "FM/chainRule_s_udg.c"
#include "FM/eulerNew/source_eulerNEW.c"
#include "FM/nsNew/source_nsNEW.c"
#include "FM/ransSAnew/source_ransSANEW.c"
#include "FM/poi/source_poisson.c"
#include "SM/LE_uq/source_leuq.c"
#include "SM/ledisp/source_ledisp.c"
#include "FM/mhd/source_mhd.c"
#include "../utilities/UDG2udg.cpp"

// Written by: C. Nguyen & P. Fernandez

void chainRuleJacobianSource(double * s_UDG, double * s_udg, double * pg, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd)
{
    /* Note: There is little to nothing to gain by using BLAS instead of for loops in this function */

    int i, j, k, l, m, n, g, len, sz2;

    Int ALEflag = app.ALEflag;
    Int flag_q = app.flag_q;

    double * Gg = &pg[3 * numPoints * nd];

    double * gg;
    double * grad_g1;
    if (ALEflag == 3) {
        gg = &pg[(3 * nd + nd * nd) * numPoints];
        grad_g1 = &pg[(3 * nd + nd * nd + 2) * numPoints];
    }
    else if (ALEflag == 2) {
        gg = new double[numPoints];
        for (i = 0; i < numPoints; i++) {
            gg[i] = 1.0;
        }
    }
    else {
        printf("Invalid value of ALEflag = %d in chainRuleJacobianSource function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
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
    for (i = 0; i < ncu; i++)
        for (k = 0; k < ncu; k++)
            for (g = 0; g < numPoints; g++) {
                s_UDG[g + k * numPoints + i * sz2] = s_udg[g + k * numPoints + i * sz2] / gg[g];

                if (ALEflag == 3) {
                    for (l = 0; l < nd; l++)
                        for (m = 0; m < nd; m++) {
                            s_UDG[g + k * numPoints + i * sz2] += grad_g1[g + l * numPoints] * Ginv[g + l * numPoints + m * numPoints * nd] *
                                                                  s_udg[g + k * numPoints + (ncu + i + m * ncu) * sz2];
                        }
                }
            }


    if (flag_q == 1) {
        for (l = 0; l < nd; l++)
            for (i = 0; i < ncu; i++)
                for (k = 0; k < ncu; k++)
                    for (g = 0; g < numPoints; g++) {
                        s_UDG[g + k * numPoints + (ncu + i + l * ncu) * sz2] = Ginv[g + l * numPoints + 0 * numPoints * nd]
                                                                               * s_udg[g + k * numPoints + (ncu + i + 0 * ncu) * sz2] / gg[g];
                        for (m = 1; m < nd; m++) {
                            s_UDG[g + k * numPoints + (ncu + i + l * ncu) * sz2] += Ginv[g + l * numPoints + m * numPoints*nd]
                                                                                    * s_udg[g + k * numPoints + (ncu + i + m * ncu) * sz2] / gg[g];
                        }
                    }
    }


    /* Deallocate dynamic memory */
    delete[] Ginv;
    if (ALEflag == 2) {
        delete[] gg;
    }
}



void nonRigidALEcorrectionSource(double *s, double *s_UDG, double *pg, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd, int computeJacobian)
{

    int i, j, k, sz2;

    double * gg = &pg[(3 * nd + nd * nd)*numPoints];

    for (i = 0; i < ncu; i++)
        for (k = 0; k < numPoints; k++) {
            s[k + i * numPoints] *= gg[k];
        }

    if (computeJacobian == 1) {
        sz2 = numPoints * ncu;
        for (j = 0; j < nc; j++)
            for (i = 0; i < ncu; i++)
                for (k = 0; k < numPoints; k++)
                    s_UDG[k + i * numPoints + j * sz2] *= gg[k];
    }
}

void sourceDriver(double * s,double * s_UDG, double * pg, double * UDG, double * ODG, 
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

    Int inc = 1;
    int i, len;

    Int ALEflag = app.ALEflag;

    double * udg;
    double * odg;
    double * s_udg;
    
    /* 1. Compute (u,q) from (U, gradX U) */
    if ((app.appname != 0 && app.appname != 1) && (ALEflag == 2 || ALEflag == 3)) {
        udg = new double[numPoints * nc];
        odg = new double[numPoints * nco];
        UDG2udg(udg, UDG, pg, app, numPoints, ncu, ncq, nc, nd);
        UDG2udg(odg, ODG, pg, app, numPoints, ncu, ncq, nco, nd);
        if (computeJacobian == 1)
            s_udg = new double[numPoints * ncu * nc];
    }
    else {
        udg = &UDG[0];
        odg = &ODG[0];
        s_udg = &s_UDG[0];
    }

    double *param = &app.physicsparam[0];
    double time = app.time;
    
    /* 2. Compute physical source */
    switch (app.appname) {
        case 0:
            source_eulerNEW(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 1:
            source_nsNEW(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 2:
            source_poisson(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;        
        case 3:
            source_ransSANEW(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 4:
            source_leuq(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;    
        case 5:
            source_ledisp(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 6:
            source_mhd(s, s_udg, pg, udg, odg, app, param, time, numPoints, nc, ncu, nd, ncd, computeJacobian);
            break;
        default:
            printf("Application not implemented (appname = %d)\n",app.appname);
            exit(-1);
    }

    /* 3. Chain rule for the Jacobian from (u,q) to (U, gradX U) */
    if ((app.appname != 0 && app.appname != 1) && (ALEflag == 2 || ALEflag == 3) && computeJacobian == 1) {
        chainRuleJacobianSource(s_UDG, s_udg, pg, app, numPoints, ncu, ncq, nc, nd);
    }

    /* 4. Non-rigid ALE correction */
    if ((app.appname != 0 && app.appname != 1) && ALEflag == 3) {
        nonRigidALEcorrectionSource(s, s_UDG, pg, app, numPoints, ncu, ncq, nc, nd, computeJacobian);
    }

    /* Deallocate dynamic memory */
    if (((app.appname != 0 && app.appname != 1)) && (ALEflag == 2 || ALEflag == 3)) {
        delete[] udg; delete[] odg;
        if (computeJacobian == 1) {
            delete[] s_udg;
        }
    }
}

#endif
