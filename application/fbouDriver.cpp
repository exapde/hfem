#ifndef __FBOUDRIVER
#define __FBOUDRIVER

#include "FM/eulerNew/fbou_eulerNEW.c"
#include "FM/nsNew/fbou_nsNEW.c"
#include "FM/ransSAnew/fbou_ransSANEW.c"
#include "FM/poi/fbou_poisson.c"
#include "SM/LE_uq/fbou_leuq.c"
#include "SM/ledisp/fbou_ledisp.c"
#include "FM/mhd/fbou_mhd.c"


// Written by: C. Nguyen & P. Fernandez

void chainRuleJacobianFbou(double *fhn_U, double *fhn_UH, double *fhn_u, double *fhn_uh, double *pg, appstruct &app,
                           int numPoints, int ncu, int ncq, int nc, int nd)
{
    /* This function should be used only if ALEflag == 2 or ALEflag == 3 */
    /* Note: There is little to nothing to gain by using BLAS instead of for loops in this function. */

    int i, j, k, l, m, n, g, sz2;
    double ggi;

    Int ALEflag = app.ALEflag;
    Int flag_q = app.flag_q;

    if (ALEflag != 2 && ALEflag != 3) {
        printf("Invalid value of ALEflag = %d in chainRuleJacobianFbou function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get Gg */
    double * Gg = &pg[3 * numPoints * nd];

    /* Get grad_g1 */
    double * grad_g1;
    if (ALEflag == 3) {
        grad_g1 = &pg[(3 * nd + nd * nd + 2) * numPoints];
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
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;

            Ginv11[i] =   Gg22[i] / ggi;
            Ginv12[i] = - Gg12[i] / ggi;
            Ginv21[i] = - Gg21[i] / ggi;
            Ginv22[i] =   Gg11[i] / ggi;
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
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;

            Ginv11[i] = (Gg22[i] * Gg33[i] - Gg23[i] * Gg32[i]) / ggi;
            Ginv12[i] = (Gg23[i] * Gg31[i] - Gg21[i] * Gg33[i]) / ggi;
            Ginv13[i] = (Gg21[i] * Gg32[i] - Gg22[i] * Gg31[i]) / ggi;
            Ginv21[i] = (Gg13[i] * Gg32[i] - Gg12[i] * Gg33[i]) / ggi;
            Ginv22[i] = (Gg11[i] * Gg33[i] - Gg13[i] * Gg31[i]) / ggi;
            Ginv23[i] = (Gg12[i] * Gg31[i] - Gg11[i] * Gg32[i]) / ggi;
            Ginv31[i] = (Gg12[i] * Gg23[i] - Gg13[i] * Gg22[i]) / ggi;
            Ginv32[i] = (Gg13[i] * Gg21[i] - Gg11[i] * Gg23[i]) / ggi;
            Ginv33[i] = (Gg11[i] * Gg22[i] - Gg12[i] * Gg21[i]) / ggi;
        }
    }


    // Chain rule for fhn_U
    sz2 = numPoints * ncu;
    for (i = 0; i < ncu; i++)
        for (k = 0; k < ncu; k++)
            for (g = 0; g < numPoints; g++) {
                if (ALEflag == 3)
                    ggi = pg[(3 * nd + nd * nd) * numPoints + g];
                else
                    ggi = 1.0;

                fhn_U[g + k * numPoints + i * sz2] = fhn_u[g + k * numPoints + i * sz2] / ggi;

                if (ALEflag == 3) {
                    for (l = 0; l < nd; l++)
                        for (m = 0; m < nd; m++) {
                            fhn_U[g + k * numPoints + i * sz2] +=
                                    grad_g1[g + l * numPoints] * Ginv[g + l * numPoints + m * numPoints * nd] *
                                    fhn_u[g + k * numPoints + (ncu + i + m * ncu) * sz2];
                        }
                }
            }


    if (flag_q == 1) {
        for (l = 0; l < nd; l++)
            for (i = 0; i < ncu; i++)
                for (k = 0; k < ncu; k++)
                    for (g = 0; g < numPoints; g++) {
                        if (ALEflag == 3)
                            ggi = pg[(3 * nd + nd * nd) * numPoints + g];
                        else
                            ggi = 1.0;

                        fhn_U[g + k * numPoints + (ncu + i + l * ncu) * sz2] =
                                Ginv[g + l * numPoints + 0 * numPoints * nd] *
                                fhn_u[g + k * numPoints + (ncu + i + 0 * ncu) * sz2] / ggi;
                        for (m = 1; m < nd; m++) {
                            fhn_U[g + k * numPoints + (ncu + i + l * ncu) * sz2] +=
                                    Ginv[g + l * numPoints + m * numPoints * nd] *
                                    fhn_u[g + k * numPoints + (ncu + i + m * ncu) * sz2] / ggi;
                        }
                    }
    }


    // Chain rule for fhn_UH
    for (i = 0; i < ncu; i++)
        for (k = 0; k < ncu; k++)
            for (g = 0; g < numPoints; g++) {
                if (ALEflag == 3)
                    ggi = pg[(3 * nd + nd * nd) * numPoints + g];
                else
                    ggi = 1.0;

                fhn_UH[g + k * numPoints + i * sz2] = fhn_uh[g + k * numPoints + i * sz2] / ggi;
            }


    /* Deallocate dynamic memory */
    delete[] Ginv;
}

//fbouDriver(fb, fb_u, fb_uh, pft, uft, uht, oft, nlt, bcs, mesh, master, app, sol, temp,  bn,  computeJacobian, ngf);
// fbouDriver(fb, fb_u, fb_uh, pft, uft, uht, oft, nlt, bcs, mesh, master, app, sol, temp,  
//                bn, ife, ie, quadchoice, computeJacobian, ngf);
void fbouDriver(double *fh, double *fhn_U, double *fhn_UH, double *pg, double *UDG, double *UHG,
                double *ODG, double *NL, double *ui, meshstruct &mesh, masterstruct &master, appstruct &app, 
        solstruct &sol, tempstruct &temp, Int ib, Int ife, Int ie, Int quadchoice, Int computeJacobian, Int numPoints)
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

    double * udg;
    double * odg;
    double * uhg;
    //double * uhg_ref;
    double * fhn_u;
    double * fhn_uh;
    double * nl;

    /* 1. Compute udg and nl in the deformed (physical) domain */
    if (ALEflag == 2 || ALEflag == 3) {
        // TODO: Need to convert avtg_p1CG as well (note that it is 2*ngf for AVflag == 10 (or maybe another number if we renumber the AV flags)
        udg = new double[numPoints * nc];
        UDG2udg(udg, UDG, pg, app, numPoints, ncu, ncq, nc, nd);
        odg = new double[numPoints * nco];
        UDG2udg(odg, ODG, pg, app, numPoints, ncu, ncq, nco, nd);
        nl = new double[numPoints * nd];
        NL2nl(nl, NL, pg, app, numPoints, nd);
        if (computeJacobian == 1)
            fhn_u = new double[numPoints * ncu * nc];
    }
    else {
        udg = &UDG[0];
        odg = &ODG[0];
        nl = &NL[0];
        fhn_u = &fhn_U[0];
    }

    /* 2. Compute uhg and nl in the deformed (physical) domain */
    if (ALEflag == 3) {
        uhg = new double[numPoints * nch];
        UH2uh(uhg, UHG, pg, app, numPoints, nch, nd);
        //uhg_ref = new double[numPoints * nch];
        //UH2uh(uhg_ref, UHGref, pg, app, numPoints, nch, nd);
        if (computeJacobian == 1)
            fhn_uh = new double[numPoints * ncu * nch];
    }
    else {
        uhg = &UHG[0];
        //uhg_ref = &UHGref[0];
        fhn_uh = &fhn_UH[0];
    }

    /* 3. Compute physical boundary conditions */
    switch (app.appname) {
        case 0:
            fbou_eulerNEW(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 1:
            fbou_nsNEW(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 2:
            fbou_poisson(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 3:
            fbou_ransSANEW(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 4:
            fbou_leuq(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 5:
            fbou_ledisp(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        case 6:
            fbou_mhd(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
            break;
        default: {
            printf("Application not implemented (appname = %d)\n",app.appname);
            exit(-1);
        }
    }

//     if (ie==0) {
//         if (app.my_rank==0) {
//             Int ngf = numPoints;
//             cout<<app.my_rank<<"  "<<ib<<"  "<<ife<<"  "<<ie<<endl;        
//             print2darray(pg,ngf,app.nd);
//             print2darray(udg,ngf,app.nc);
//             print2darray(uhg,ngf,app.nch);        
//             print2darray(nl,ngf,app.nd);
//             //print1darray(ui,app.nch);
//             print2darray(fh,ngf,app.nch);
//             //print3darray(fhn_u,ngf,app.nch,app.nc);
//             //print3darray(fhn_uh,ngf,app.nch,app.nch);
//         }
//     }
    
    /* 4. Chain rule for the Jacobian from (u,q) to (U, gradX U) */
    if ((ALEflag == 2 || ALEflag == 3) && computeJacobian == 1) {
        chainRuleJacobianFbou(fhn_U, fhn_UH, fhn_u, fhn_uh, pg, app, numPoints, ncu, ncq, nc, nd);
    }

    /* Deallocate dynamic memory */
    if (ALEflag == 2 || ALEflag == 3) {
        delete[] udg;
        delete[] odg;
        delete[] nl;
        if (computeJacobian == 1) {
            delete[] fhn_u;
        }
    }
    if (ALEflag == 3) {
        delete[] uhg;
        //delete[] uhg_ref;
        if (computeJacobian == 1) {
            delete[] fhn_uh;
        }
    }
}


#endif
