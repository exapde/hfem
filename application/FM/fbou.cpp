#ifndef __FBOU
#define __FBOU

//#include "../fhatDriver.cpp"
#include "getanNEW.c"
#include "getan.c"

// Written by: C. Nguyen & P. Fernandez


void fbou(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *uhg, double *odg,
          double *nl, double *ui, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
          tempstruct &temp, Int ib, Int ife, Int ie, Int quadchoice, int computeJacobian, int numPoints)
{
    // With ALE, this function assumes udg, uhg and nl are in the physical (deformed) domain.
    
    // TODO: Not all capabilities for ALE are included (with deforming mesh it will not work)
    // TODO: Only SA turbulence model is available
    // TODO: Implement symmetry BC (like slip but with different BC for SA)
    
    
    int ng, ncu, ncq, nch, nc, nd, ncd, nco;
    nd = master.nd;
    ncd = app.ncd;    
    nc = app.nc;
    ncu = app.ncu;
    ncq = app.ncq;
    nch = app.nch;
    nco = app.nco;    
    ng = numPoints;

    double gam, epslm, gam1, Re, Pr, Minf, M2_infty, tau;
    double nx, ny, tm;
    int i, j, k, l, ii, m, n, nm, nm2, nk, nn;
    int sz2 = ng * nch, turbEqs = nch - (2+nd);
    
    double time = app.time;
    double *param = &app.physicsparam[0];
    gam   = param[0];
    epslm = param[1];
    Re    = param[2];
    Pr    = param[3];
    Minf  = param[4];
    tau   = param[5];
    gam1  = gam-1.0;
    M2_infty    = Minf*Minf;
    
    // Get mesh velocity:
    double * Vg;
    if (ib==2 || ib==3 || ib == 7 || ib == 11) {
        if (app.ALEflag == 0) {     // Fixed mesh
            Vg = new double[ng * nd];
            for (i = 0; i < ng * nd; i++)
                Vg[i] = 0.0;
        }
        else {                      // Moving mesh
            Vg = &pg[(2 * nd) * ng];
//            for (j = 0; j < nd; j++)
//                for (i = 0; i < ng; i++)
//                    Vg[i + j * ng] = pg[(2 * nd + j) * ng + i];
        }
    }
    
//     double angle, W_z;
//     if (time < 1.0) {
//         angle = 0.0;
//         W_z = 0.0;
//     }
//     else if (time < 2.0) {
//         angle = 0.2*0.5*(time-1.0)*(time-1.0);
//         W_z = 0.2*(time-1.0);
//     }
//     else {
//         angle = 0.2*(0.5 + (time-2.0));
//         W_z = 0.2;
//     }
//     
//     Vg = new double[ng * nd];
//     for (i = 0; i < ng; i++) {
//         Vg[0*ng+i] =   W_z*pg[1 * ng + i];
//         Vg[1*ng+i] = - W_z*pg[0 * ng + i];
//     }
    
//     double angle, W_z;
//     Vg = new double[ng * nd];
//     for (i = 0; i < ng; i++) {
//         if (pg[0 * ng + i]*pg[0 * ng + i]+pg[1 * ng + i]*pg[1 * ng + i] < 1.0) {     // Inner cylinder
//             if (time < 1.0) {
//                 angle = 0.1*time + 0.5*0.1*time*time;
//                 W_z = 0.1 + 0.1*time;
//             }
//             else {
//                 angle = 0.15 + 0.2*(time-1.0);
//                 W_z = 0.2;
//             }
//         }
//         else {          // Outer cylinder
//             if (time < 1.0) {
//                 angle = 0.1*time + 0.5*0.2*time*time;
//                 W_z = 0.1 + 0.2*time;
//             }
//             else {
//                 angle = 0.2 + 0.3*(time-1.0);
//                 W_z = 0.3;
//             }
//         }
//         Vg[0*ng+i] =   W_z*pg[1 * ng + i];
//         Vg[1*ng+i] = - W_z*pg[0 * ng + i];
//     }
    
//     double alpha, W_z;
//     Vg = new double[ng * nd];
//     for (i = 0; i < ng; i++) {
//         if (time < 1.0) {
//             alpha = 0.1*time + 0.5*0.1*time*time;
//             W_z = 0.1 + 0.1*time;
//         }
//         else {
//             alpha = 0.15 + 0.2*(time-1.0);
//             W_z = 0.2;
//         }
//         Vg[0*ng+i] =   W_z * pg[1 * ng + i];
//         Vg[1*ng+i] = - W_z * pg[0 * ng + i];
//     }
    
    if (ib==1) { // Freestream (regular far-field)
        
        double *an = new double[ng * nch * nch];
        double *An = new double[ng * nch * nch];
        double *anm = new double[ng * nch * nch * nch];
        double *Anm = new double[ng * nch * nch * nch];
        double *ui_g = new double[nch];
        
        double Vn;
        
        //getanNEW(an, anm, uhg, pg, nl, app, param, 0, ng, nd, computeJacobian);
        //getanNEW(An, Anm, uhg, pg, nl, app, param, 1, ng, nd, computeJacobian);
        getan(an, anm, uhg, nl, param, 0, ng, nch, nd);
        getan(An, Anm, uhg, nl, param, 1, ng, nch, nd);       
        
        for (i = 0; i < ng * nch; i++)
            fh[i] = 0.0;

        if (computeJacobian == 1) {
            for (i = 0; i < ng * nch * nc; i++)
                fh_u[i] = 0.0;

            for (i = 0; i < ng * nch * nch; i++)
                fh_uh[i] = 0.0;
        }
        
        // Boundary conditions for mass, momentum and energy equations:
        for (k = 0; k < (2 + nd); k++)
            for (j = 0; j < (2 + nd); j++)
                for (i = 0; i < ng; i++) {
                    
                    if (app.rotatingFrame == 0) {
                        for (n = 0; n < nch; n++)
                            ui_g[n] = ui[n];
                    }
                    else if (app.rotatingFrame == 1) {
                        error("Need to define Euler angles and angular velocity in fbou.cpp.\n");
                        ui_g[0] = ui[0];
                        ui_g[1+nd] = 0.0;
                        for (n = 0; n < nd; n++) {
//                             ui_g[1+n] = cos(alpha) + W_z*pg[1 * ng + i]; / - sin(alpha) - W_z*pg[0 * ng + i];
                            ui_g[1+nd] += ui_g[1+n]*ui_g[1+n];
                        }
                        ui_g[1+nd] = sqrt(ui_g[1+nd]) + 1/(gam*gam1*M2_infty);
                        for (n = 2+nd; n < nch; n++)
                            ui_g[n] = ui[n];
// // // //                         ui[1] =   cos(alpha) + W_z*pg[1 * ng + i];
// // // //                         ui[2] = - sin(alpha) - W_z*pg[0 * ng + i];
// // // //                         ui[3] = sqrt(ui[1]*ui[1] + ui[2]*ui[2]) + 1/(gam*gam1*M2_infty);
                    }
                    
                    nm = k * ng * (nd + 2) + j * ng + i;
                    nm2 = k * sz2 + j * ng + i;
                    nk = k * ng + i;
                    fh[j * ng + i] += (an[nm] + An[nm]) * (udg[nk] - uhg[nk]) -
                                      (an[nm] - An[nm]) * (ui_g[k] - uhg[nk]);

                    if (computeJacobian == 1) {
                        fh_u[nm2] = an[nm] + An[nm];
                        fh_uh[nm2] = -2.0 * An[nm];
                        for (n = 0; n < (nd + 2); n++) {
                            nk = k * (nd + 2) * (nd + 2) * ng + n * (nd + 2) * ng + j * ng + i;
                            nn = n * ng + i;
                            fh_uh[nm2] += (anm[nk] + Anm[nk]) * (udg[nn] - uhg[nn]) -
                                          (anm[nk] - Anm[nk]) * (ui_g[n]   - uhg[nn]);
                        }
                    }
                }
        
        // Boundary conditions for the SA equation:
        for (j=(2+nd); j<nch; j++)
            for (i=0; i<ng; i++) {
                // Compute normal component of linear momentum
                Vn = 0.0;
                for (ii = 0; ii < nd; ii++)
                    Vn += nl[ii * ng + i] * uhg[i + (1 + ii) * ng];

                if (Vn < 0.0) {       // Inflow (Dirichlet for rN)
                    fh[j * ng + i] = ui_g[j] - uhg[i + j * ng];
                    if (computeJacobian == 1)
                        fh_uh[j * sz2 + j * ng + i] = -1.0;
                }
                else {                // Outflow (extrapolation of rN)
                    fh[j * ng + i] = udg[i + j * ng] - uhg[i + j * ng];
                    if (computeJacobian == 1) {
                        fh_u[j * sz2 + j * ng + i] = 1.0;
                        fh_uh[j * sz2 + j * ng + i] = -1.0;
                    }
                }
            }
        
        delete[] an; delete[] An; delete[] anm; delete[] Anm; delete[] ui_g;
    }
    else if (ib==2) { // Non-slip, adiabatic wall
        
        // Boundary condition for energy equation (adiabatic condition):
        Int fhatExpression = 0;
        fhatDriver(fh, fh_u, fh_uh, pg, udg, uhg, odg, nl, mesh, master, app, sol, temp, fhatExpression, ie, quadchoice, computeJacobian, numPoints);       // TODO: Impose energy BC in another way, since this won't work for ALE
        
        for (i=0; i<ng; i++) {
            // Boundary condition for mass equation (extrapolation):
            fh[0*ng+i] = udg[0*ng+i]-uhg[0*ng+i];

            // Boundary condition for momentum equation (non-slip condition):
            for (j=0; j<nd; j++)
                fh[(1+j)*ng+i] = udg[0*ng+i]*Vg[j*ng+i] - uhg[(1+j)*ng+i];

// // // //            // Boundary condition for energy equation (adiabatic wall). TODO: Only valid for tau > 0:
// // // //            fh[(1+nd)*ng+i] = 0.0;
// // // //            for (j=0; j<nd; j++)
// // // //                fh[(1+nd)*ng+i] += (uhg[0*ng+i]*udg[((1+j)*ncu+1+nd)*ng+i] - uhg[(1+nd)*ng+i]*udg[(1+j)*ncu*ng+i]) * nl[j * ng + i];        // (rEx*r - rE*rx)*nx + (rEy*r - rE*ry)*ny ( +(rEz*r - rE*rz)*nz )
// // // //            fh[(1+nd)*ng+i] *= gam/(Pr*Re);       // A factor of the form mu/r^2 has been omitted since it is O(1) and only increases the non-linearity of the BC
// // // //            fh[(1+nd)*ng+i] += tau*(udg[(1+nd)*ng+i]-uhg[(1+nd)*ng+i]);

            // Boundary condition for SA equation (homogeneous Dirichlet):
            for (j=(nd+2); j<nch; j++)
                fh[j*ng+i] = -uhg[j*ng+i];

            // Derivatives of BCs
            if (computeJacobian == 1) {
                for (k=0; k<nch; k++)
                    for (j=0; j<nch; j++)
                        if (j != (1+nd))
                            fh_uh[k * sz2 + j * ng + i] = (j==k) ? -1.0 : 0.0;

//                for (k=0; k<nch; k++)
//                    fh_uh[k * sz2 + (1+nd) * ng + i] = 0.0;
//                for (ii=0; ii<nd; ii++) {
//                    fh_uh[0 * sz2 + (1+nd) * ng + i] += udg[((1+ii)*ncu+1+nd)*ng+i] * nl[ii * ng + i] * gam/(Pr*Re);
//                    fh_uh[(1+nd) * sz2 + (1+nd) * ng + i] -= udg[(1+ii)*ncu*ng+i] * nl[ii * ng + i] * gam/(Pr*Re);
//                }
//                fh_uh[(1+nd) * sz2 + (1+nd) * ng + i] -= tau;

                for (k = 0; k < nc; k++)
                    for (j = 0; j < nch; j++)
                        if (j != (1+nd))
                            fh_u[k * sz2 + j * ng + i] = 0.0;

                fh_u[0 * sz2 + 0 * ng + i] = 1.0;

                for (j = 0; j < nd; j++)
                    fh_u[0 * sz2 + (1+j) * ng + i] = Vg[j*ng+i];

//                for (k=0; k<nc; k++)
//                    fh_u[k * sz2 + (1+nd) * ng + i] = 0.0;
//                for (ii=0; ii<nd; ii++) {
//                    fh_u[((1+ii)*ncu+1+nd) * sz2 + (1+nd) * ng + i] += uhg[0*ng+i] * nl[ii * ng + i] * gam/(Pr*Re);
//                    fh_u[((1+ii)*ncu) * sz2 + (1+nd) * ng + i] -= uhg[(1+nd)*ng+i] * nl[ii * ng + i] * gam/(Pr*Re);
//                }
//                fh_u[(1+nd) * sz2 + (1+nd) * ng + i] += tau;
            }
        }
    }
    else if (ib==3) { // Isothermal (T_wall = ratio * T_inf), non-slip wall
//         double kinEnrgy, kinEnrgyInf = 0.0, T_infty;
//         for (j = 0; j < nd; j++)
//             kinEnrgyInf += 0.5*ui[1+j]*ui[1+j]/ui[0];
//         T_infty = (ui[1+nd] - kinEnrgyInf) / ui[0];
//         
//         double ratio = 500.0/200.0;     // 500.0/200.0: Value for cylinder at M = 18, Re = 378k in Barer's thesis
//         double T_w = ratio * T_infty;
// 
//         double *an = new double[ng * nch * nch];
//         double *An = new double[ng * nch * nch];
//         double *anm = new double[ng * nch * nch * nch];
//         double *Anm = new double[ng * nch * nch * nch];
//         double *ui_g = new double[nch];
//         
//         double Vn;
//         
//         getanNEW(an, anm, uhg, pg, nl, app, param, 0, ng, nd, computeJacobian);
//         getanNEW(An, Anm, uhg, pg, nl, app, param, 1, ng, nd, computeJacobian);
// 
//         for (i = 0; i < ng * nch; i++)
//             fh[i] = 0.0;
// 
//         if (computeJacobian == 1) {
//             for (i = 0; i < ng * nch * nc; i++)
//                 fh_u[i] = 0.0;
//             for (i = 0; i < ng * nch * nch; i++)
//                 fh_uh[i] = 0.0;
//         }
//         
//         // Boundary conditions for mass, momentum and energy equations:
//         for (k = 0; k < (2 + nd); k++)
//             for (j = 0; j < (2 + nd); j++)
//                 for (i = 0; i < ng; i++) {
//                     kinEnrgy = 0.0;
//                     for (n = 0; n < nd; n++)
//                         kinEnrgy += 0.5 * udg[(1+n)*ng+i] * udg[(1+n)*ng+i] / udg[0*ng+i];
//                     
//                     ui_g[0] = (udg[(1+nd)*ng+i] - kinEnrgy) / T_w;//udg[0*ng+i];//
//                     for (n = 0; n < nd; n++)
//                         ui_g[1+n] = 0.0;
//                     ui_g[1+nd] = (udg[(1+nd)*ng+i] - kinEnrgy);//udg[0*ng+i] * T_w;//
//                     
//                     nm = k * ng * (nd + 2) + j * ng + i;
//                     
//                     
//                     
//                     if (j==k) {
// //                         an[nm] += 1.0;
//                         An[nm] += 1.0;
//                     }
// //                     else {
// //                         an[nm] = 0.0;
// //                         An[nm] = 0.0;
// //                     }
//                     
//                     
//                     
//                     nm2 = k * sz2 + j * ng + i;
//                     nk = k * ng + i;
// // // //                     fh[j * ng + i] += An[nm] * (udg[nk] - uhg[nk]) +
// // // //                                       An[nm] * (ui_g[k] - uhg[nk]);
//                     fh[j * ng + i] += (an[nm] + An[nm]) * (udg[nk] - uhg[nk]) -
//                                       (an[nm] - An[nm]) * (ui_g[k] - uhg[nk]);
//                     
//                     if (computeJacobian == 1) {
//                         fh_u[nm2] += an[nm] + An[nm];
// // // //                         fh_u[nm2] += An[nm];
//                         if (k == 0) {
//                             fh_u[0*sz2+j*ng+i] -= (an[nm] - An[nm]) * kinEnrgy / (udg[0*ng+i] * T_w);
//                             for (n = 0; n < nd; n++)
//                                 fh_u[(1+n)*sz2+j*ng+i] += (an[nm] - An[nm]) * udg[(1+n)*ng+i] / (udg[0*ng+i] * T_w);
//                             fh_u[(1+nd)*sz2+j*ng+i] -= (an[nm] - An[nm]) / T_w;
// //                             fh_u[0*sz2+j*ng+i] -= (an[nm] - An[nm]);
//                         }
//                         if (k == 1+nd) {
//                             fh_u[0 * sz2 + j * ng + i] -= (an[nm] - An[nm]) * kinEnrgy / udg[0*ng+i];
//                             for (n = 0; n < nd; n++)
//                                 fh_u[(1+n) * sz2 + j * ng + i] += (an[nm] - An[nm]) * udg[(1+n)*ng+i] / udg[0*ng+i];
//                             fh_u[(1+nd) * sz2 + j * ng + i] -= (an[nm] - An[nm]);
// //                             fh_u[0 * sz2 + j * ng + i] -= (an[nm] - An[nm]) * T_w;
//                         }
//                         
//                         fh_uh[nm2] += -2.0 * An[nm];
// // // // //                         fh_uh[nm2] += - 2.0 * An[nm];
//                         for (n = 0; n < (nd + 2); n++) {
//                             nk = k * (nd + 2) * (nd + 2) * ng + n * (nd + 2) * ng + j * ng + i;
//                             nn = n * ng + i;
//                             fh_uh[nm2] += (anm[nk] + Anm[nk]) * (udg[nn] - uhg[nn]) -
//                                           (anm[nk] - Anm[nk]) * (ui_g[n] - uhg[nn]);
// // // // //                             fh_uh[nm2] += Anm[nk] * (udg[nn] - uhg[nn]) +
// // // // //                                           Anm[nk] * (ui_g[n] - uhg[nn]);
//                         }
//                     }
//                 }
//         
//         // Boundary condition for SA equation (homogeneous Dirichlet):
//         for (j = (nd+2); j < nch; j++)
//             for (i = 0; i < ng; i++) {
//                 fh[j*ng+i] = - uhg[j*ng+i];
//                 if (computeJacobian == 1)
//                     fh_uh[j * sz2 + j * ng + i] = -1.0;
//             }
//         
//         delete[] an; delete[] An; delete[] anm; delete[] Anm; delete[] ui_g;
        
        
        
        double kinEnrgy, kinEnrgyInf = 0.0, T_infty;
        for (j=0; j<nd; j++)
            kinEnrgyInf += 0.5*ui[1+j]*ui[1+j]/ui[0];
        T_infty = (ui[1+nd] - kinEnrgyInf) / ui[0];
        
        double ratio = 500.0/200.0;     // 500.0/200.0: Value for cylinder at M = 18, Re = 378k in Barer's thesis
        double T_w = ratio * T_infty;

        for (i=0; i<ng; i++) {
            // Boundary condition for mass equation (extrapolation):
fh[0*ng+i] =  udg[0*ng+i]-uhg[0*ng+i];

            // Boundary condition for momentum equation (non-slip condition):
            for (j=0; j<nd; j++)
                fh[(1+j)*ng+i] = udg[0*ng+i]*Vg[j*ng+i] - uhg[(1+j)*ng+i];

            // Boundary condition for energy equation (isothermal condition):
            kinEnrgy = 0.0;
            for (j=0; j<nd; j++)
                kinEnrgy += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
fh[(1+nd)*ng+i] = ratio * (ui[1+nd] - kinEnrgyInf) / ui[0] - (uhg[(1+nd)*ng+i] - kinEnrgy) / uhg[0*ng+i];       // This line replaces the following line from the MATLAB implementation
//            fh[(1+nd)*ng+i] = ui[1+nd] - gam*gam1*M2_infty*uhg[(1+nd)*ng+i]/uhg[0*ng+i];          // TODO: While u=v=w=0 on the wall (due to non-slip condition), we may want to include that also here through the kinetic energy
                    
            // Boundary condition for SA equation (homogeneous Dirichlet):
            for (j=(nd+2); j<nch; j++)
                fh[j*ng+i] = -uhg[j*ng+i];

            // Derivatives of BCs
            if (computeJacobian == 1) {
                for (k=0; k<nch; k++)
                    for (j=0; j<nch; j++)
                        if (j != (1+nd))
                            fh_uh[k * sz2 + j * ng + i] = (j==k) ? -1.0 : 0.0;

                for (k=(2+nd); k<nch; k++)      // This piece of code replaces the following commented code from the MATLAB implementation
                    fh_uh[k * sz2 + (1+nd) * ng + i] = 0.0;
fh_uh[0 * sz2 + (1+nd) * ng + i] = uhg[(1+nd)*ng+i] / (uhg[0*ng+i] * uhg[0*ng+i]);

                for (j=0; j<nd; j++) {
fh_uh[0 * sz2 + (1 + nd) * ng + i] -=
uhg[(1 + j) * ng + i] * uhg[(1 + j) * ng + i] / (uhg[0 * ng + i] * uhg[0 * ng + i] * uhg[0 * ng + i]);
fh_uh[(1 + j) * sz2 + (1 + nd) * ng + i] =
uhg[(1 + j) * ng + i] / (uhg[0 * ng + i] * uhg[0 * ng + i]);
                }
fh_uh[(1+nd)*sz2 + (1+nd)*ng + i] = -1.0/uhg[0*ng+i];
//                for (k=0; k<nch; k++)
//                    fh_uh[k * sz2 + (1+nd) * ng + i] = 0.0;
//                fh_uh[0*sz2 + (1+nd) * ng + i] = gam * gam1 * M2_infty * uhg[(1+nd) * ng + i] / (uhg[0 * ng + i] * uhg[0 * ng + i]);
//                fh_uh[(1+nd)*sz2 + (1+nd)*ng + i] = - gam * gam1 * M2_infty / uhg[0 * ng + i];

                for (k = 0; k < nc; k++)
                    for (j = 0; j < nch; j++)
                            fh_u[k * sz2 + j * ng + i] = 0.0;
fh_u[0 * sz2 + 0 * ng + i] = 1.0;
                for (j = 0; j < nd; j++)
                    fh_u[0 * sz2 + (1+j) * ng + i] = Vg[j*ng+i];
            }
        }
    }
    else if (ib==4) { // Pressure far-field

        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *anm  = new double [ng*nch*nch*nch];
        double *Anm  = new double [ng*nch*nch*nch];
        double *uinf  = new double [ng*nch];
        double *uinfu  = new double [ng*nch*nch];

        double pinf, kinEnrgyInf, kinEnrgy, Vn;

        //getanNEW(an, anm, uhg, pg, nl, app, param, 0, ng, nd, computeJacobian);
        //getanNEW(An, Anm, uhg, pg, nl, app, param, 1, ng, nd, computeJacobian);
        getan(an, anm, uhg, nl, param, 0, ng, nch, nd);
        getan(An, Anm, uhg, nl, param, 1, ng, nch, nd);       
        
        // Get desired pressure at infinity:
        kinEnrgyInf = 0.0;
        for (j=0; j<nd; j++)
            kinEnrgyInf += 0.5*ui[1+j]*ui[1+j]/ui[0];
        pinf = (gam-1) * (ui[1+nd]-kinEnrgyInf);

        // Get desired infinity conditions and derivatives w.r.t. u:
//         for (j=0; j<nch; j++)
//             for (i=0; i<ng; i++)
//                 uinf[i+j*ng] = ui[j];  // THIS IS A BUG // This loop replaces the next loop from the MATLAB implementation        
       for (i=0; i<ng*nch; i++)
           uinf[i] = udg[i];
        for (i=0; i<ng; i++) {
            n = (nch-1-turbEqs)*ng+i;
            kinEnrgy = 0;
            for (j=0; j<nd; j++)
                kinEnrgy += 0.5*udg[(j+1)*ng+i]*udg[(j+1)*ng+i]/udg[0*ng+i];
            uinf[n] = pinf/(gam-1) + kinEnrgy;
        }

        if (computeJacobian == 1) {
            for (k=0; k<nch; k++)
                for (j=0; j<nch; j++)
                    for (i=0; i<ng; i++) {
                        nm = k*sz2+j*ng+i;
                        //uinfu[nm] = 0.0;  THIS IS A BUG!     // This line replaces the next line from the MATLAB implementation
                        uinfu[nm] = (j==k) ? 1.0 : 0.0;
                    }

            for (i = 0; i < ng; i++) {
                nm = 0 * sz2 + (1 + nd) * ng + i;
                uinfu[nm] = 0.0;
                for (j = 0; j < nd; j++)
                    uinfu[nm] -= 0.5 * udg[(j + 1) * ng + i] * udg[(j + 1) * ng + i] / (udg[0*ng+i] * udg[0*ng+i]);

                for (j = 0; j < nd; j++) {
                    nm = (j + 1) * sz2 + (1 + nd) * ng + i;
                    uinfu[nm] = udg[(j + 1) * ng + i] / udg[0*ng+i];
                }

                nm = (1 + nd) * sz2 + (1 + nd) * ng + i;
                uinfu[nm] = 0.0;
            }
        }

        // Initialize fh, fh_u and fh_uh to zero:
        for (i=0; i<ng*nch; i++)
            fh[i] = 0.0;

        if (computeJacobian == 1) {
            for (i = 0; i < ng * nch * nc; i++)
                fh_u[i] = 0.0;

            for (i = 0; i < ng * nch * nch; i++)
                fh_uh[i] = 0.0;
        }
        
        // Boundary conditions for mass, momentum and energy equations:
        int na, nb;
        for (k = 0; k < (2 + nd); k++)
            for (j = 0; j < (2 + nd); j++)
                for (i = 0; i < ng; i++) {
                    nm = k * ng * (nd + 2) + j * ng + i;
                    nm2 = k * sz2 + j * ng + i;
                    nk = k * ng + i;
                    fh[j * ng + i] += (an[nm] + An[nm])*(udg[nk]  - uhg[nk]) -
                                      (an[nm] - An[nm])*(uinf[nk] - uhg[nk]);

                    if (computeJacobian == 1) {
                        fh_u[nm2] = an[nm] + An[nm];
                        fh_uh[nm2] = -2.0 * An[nm];
                        for (n = 0; n < (nd + 2); n++) {
                            nk = k * (nd + 2) * (nd + 2) * ng + n * (nd + 2) * ng + j * ng + i;
                            nn = n * ng + i;
                            na = n * (nd + 2) * ng + j * ng + i;
                            nb = k * nch * ng + n * ng + i;
                            fh_u[nm2]  -= (an[na] - An[na]) * uinfu[nb];
                            fh_uh[nm2] += (anm[nk] + Anm[nk]) * (udg[nn]  - uhg[nn]) -
                                          (anm[nk] - Anm[nk]) * (uinf[nn] - uhg[nn]);
                        }
                    }
                }

        // Boundary conditions for SA equation
        for (j=(2+nd); j<nch; j++)
            for (i=0; i<ng; i++) {
                // Compute normal component of linear momentum:
                Vn = 0.0;
                for (ii=0; ii<nd; ii++)
                    Vn += nl[ii*ng+i]*uhg[i+(1+ii)*ng];

                if (Vn < 0.0) {       // Inflow (Dirichlet for rN)
                    fh[j * ng + i] = ui[j] - uhg[i + j * ng];
                    if (computeJacobian == 1)
                        fh_uh[j * sz2 + j * ng + i] = -1.0;
                }
                else {                // Outflow (extrapolation of rN)
                    fh[j * ng + i] = udg[i + j * ng] - uhg[i + j * ng];
                    if (computeJacobian == 1) {
                        fh_u[j * sz2 + j * ng + i] = 1.0;
                        fh_uh[j * sz2 + j * ng + i] = -1.0;
                    }
                }
            }

        delete[] an; delete[] An; delete[] anm; delete[] Anm; delete[] uinf; delete[] uinfu;
    }
    else if (ib == 5) { // Dirichlet. // TODO: Should we use tau here?
        for (j=0; j<nch; j++)
            for (i=0; i<ng; i++)
                fh[i+j*ng] = tau*(ui[j] - uhg[i+j*ng]);

        if (computeJacobian == 1) {
            for (i=0; i<ng*nch*nc; i++)
                fh_u[i] = 0.0;

            for (k=0; k<nch; k++)
                for (j=0; j<nch; j++)
                    for (i=0; i<ng; i++)
                        fh_uh[i+j*ng+k*sz2] = (j == k) ? -tau : 0.0;
        }
    }
    else if (ib == 6) { // Extrapolation. // TODO: Should we use tau here?
        for (j=0; j<nch; j++)
            for (i=0; i<ng; i++)
                fh[i+j*ng] = udg[i+j*ng] - uhg[i+j*ng];

        if (computeJacobian == 1) {
            for (k = 0; k < nc; k++)
                for (j = 0; j < nch; j++)
                    for (i = 0; i < ng; i++)
                        fh_u[i + j * ng + k * sz2] = (j == k) ? 1.0 : 0.0;

            for (k = 0; k < nch; k++)
                for (j = 0; j < nch; j++)
                    for (i = 0; i < ng; i++)
                        fh_uh[i + j * ng + k * sz2] = (j == k) ? -1.0 : 0.0;
        }
    }
    else if (ib==7) {       // Slip wall (with temperature extrapolation)
        double un;
        double Vgn;

        for (i=0; i<ng; i++) {
            // Boundary condition for mass equation (extrapolation):
            fh[0*ng+i] =  udg[0*ng+i] - uhg[0*ng+i];

            // Boundary condition for momentum equation (slip condition):
            un = 0.0;
            Vgn = 0.0;
            for (j=0; j<nd; j++) {
                un += nl[j * ng + i] * udg[i + (1 + j) * ng];
                Vgn += nl[j * ng + i] * Vg[i + j * ng];
            }
            for (j=0; j<nd; j++)
                fh[(1+j)*ng+i] = udg[(j+1)*ng+i] - un*nl[j * ng + i] - (uhg[(1+j)*ng+i] - udg[0*ng+i]*Vgn*nl[j * ng + i]);        // TODO: Should we use udg[0*ng+i] or uhg[0*ng+i] here?

            // Boundary condition for energy equation (extrapolation):          // TODO: What BC should we impose for the energy equation?
            fh[(1+nd)*ng+i] = udg[(1+nd)*ng+i]-uhg[(1+nd)*ng+i];

            // Boundary condition for SA equation (extrapolation):
            for (j=(nd+2); j<nch; j++)
                fh[j*ng+i] = udg[j*ng+i] - uhg[j*ng+i];     // TODO: What BC for the SA variable should we use with slip walls?

            // Derivatives of BCs
            if (computeJacobian == 1) {
                for (k = 0; k < nch; k++)
                    for (j = 0; j < nch; j++)
                        fh_uh[k * sz2 + j * ng + i] = (j == k) ? -1.0 : 0.0;


                for (k = 0; k < nc; k++)
                    fh_u[k * sz2 + 0 * ng + i] = (0 == k) ? 1.0 : 0.0;

                for (j = 0; j < nd; j++) {
                    fh_u[0 * sz2 + (1 + j) * ng + i] = Vgn * nl[j * ng + i];
                    for (k = 0; k < nd; k++)
                        fh_u[(1 + k) * sz2 + (1 + j) * ng + i] = -nl[j * ng + i] * nl[k * ng + i];
                    fh_u[(1 + j) * sz2 + (1 + j) * ng + i] += 1.0;
                    for (k = (nd + 1); k < nc; k++)
                        fh_u[k * sz2 + (1 + j) * ng + i] = 0.0;
                }

                for (k = 0; k < nc; k++)
                    fh_u[k * sz2 + (1+nd) * ng + i] = ((1+nd) == k) ? 1.0 : 0.0;

                for (k = 0; k < nc; k++)
                    for (j = (nd + 2); j < nch; j++)
                        fh_u[k * sz2 + j * ng + i] = (j == k) ? 1.0 : 0.0;
            }
        }
    }
   else if (ib==8) { // Subsonic turbomachinery inlet: Impose p0, T0, alpha, theta. Extrapolate R- = u - 2c/(gamma-1)
       // TODO: Code ALE.
        
        double Vn, p0_target, T0_target, alpha_target, theta_target, aux1, aux2;
        double q, p, M2, p02p, p0;
        double * dqdu = new double[2+nd];
        double * dpdu = new double[2+nd];
        double * dM2du = new double[2+nd];
        double * dp02pdu = new double[2+nd];
        double * dp0du = new double[2+nd];

       p0_target = ui[0];
       T0_target = ui[1];
       alpha_target = ui[2];
       if (nd == 3)
           theta_target = ui[3];

       for (i = 0; i < ng; i++) {
            //////
            q = 0.0;
            for (j = 0; j < nd; j++)
                q += 0.5 * uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i] / uhg[0*ng+i];
            p = gam1*(uhg[(1+nd)*ng+i] - q);
            M2 = (2.0*q) / (gam*p);
            p02p = pow(1.0 + 0.5*gam1*M2, gam/gam1);
            p0 = p*p02p;
            fh[0*ng+i] = p0 - p0_target;
            fh[0*ng+i] *= M2_infty;
            /////
           
//            fh[0*ng+i] = (gam-1.0)*uhg[(1+nd)*ng+i] - p0_target;
           fh[1*ng+i] = uhg[(1+nd)*ng+i]/uhg[0*ng+i] - T0_target;
           for (j = 0; j < nd; j++) {
//                fh[0*ng+i] += 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               fh[1*ng+i] -= 0.5*((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           }
//            fh[0*ng+i] *= M2_infty;
           fh[1*ng+i] *= M2_infty;

           fh[2*ng+i] = uhg[1*ng+i]*sin(alpha_target) - uhg[2*ng+i]*cos(alpha_target);

           if (nd == 3)
               fh[3*ng+i] = sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i])*cos(theta_target) - uhg[3*ng+i];

           // Extrapolate R- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: R- = u_n - 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[(1+nd)*ng+i] = - 2.0*sqrt((gam/(gam-1.0))*aux1) + 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[(1+nd)*ng+i] += -nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate J- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: J- = u_n - sqrt(gam*e/(gam-1))
           }
           
           // Impose turbulent variables (if necessary):
           for (j = (2+nd); j < ncu; j++) {
               fh[j*ng+i] = ui[j] - uhg[i + j * ng];
           }
       }

       if (computeJacobian == 1) {
           for (i = 0; i < ng * nch * nc; i++)
               fh_u[i] = 0.0;

           for (i = 0; i < ng * nch * nch; i++)
               fh_uh[i] = 0.0;
           
           for (i = 0; i < ng; i++) {
              //////
                q = 0.0;
                for (j = 0; j < nd; j++)
                    q += 0.5 * uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i] / uhg[0*ng+i];
                dqdu[0] = - q / uhg[0*ng+i];
                for (j = 0; j < nd; j++)
                    dqdu[1+j] = uhg[(1+j)*ng+i] / uhg[0*ng+i];
                dqdu[1+nd] = 0.0;
                
                p = gam1*(uhg[(1+nd)*ng+i] - q);
                for (j = 0; j < (1+nd); j++)
                    dpdu[j] = - gam1*dqdu[j];
                dpdu[1+nd] = gam1;

                M2 = (2.0*q) / (gam*p);
                for (j = 0; j < (2+nd); j++)
                    dM2du[j] = (2.0 / gam) * ((dqdu[j]*p - q*dpdu[j]) / (p*p));

                p02p = pow(1.0 + 0.5*gam1*M2, gam/gam1);
                for (j = 0; j < (2+nd); j++)
                    dp02pdu[j] = 0.5 * gam * pow(1.0 + 0.5*gam1*M2, 1.0/gam1) * dM2du[j];

                p0 = p*p02p;
                for (j = 0; j < (2+nd); j++)
                    dp0du[j] = dpdu[j] * p02p + p * dp02pdu[j];
               
               for (j = 0; j < (2+nd); j++)
                   fh_uh[j * sz2 + 0 * ng + i] = M2_infty * dp0du[j];
               /////
               
//                fh_uh[(1+nd) * sz2 + 0 * ng + i] = gam-1.0;
               fh_uh[(1+nd) * sz2 + 1 * ng + i] = 1.0/uhg[0*ng+i];
               fh_uh[0 * sz2 + 1 * ng + i] = - uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               for (j = 0; j < nd; j++) {
//                    fh_uh[0 * sz2 + 0 * ng + i] -= 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
//                    fh_uh[(1+j) * sz2 + 0 * ng + i] += (2.0-gam)*uhg[(1+j)*ng+i]/uhg[0*ng+i];
                   fh_uh[0 * sz2 + 1 * ng + i] += ((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[(1+j) * sz2 + 1 * ng + i] -= ((gam-1.0)/gam)*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               }
//                fh_uh[0*sz2 + 0*ng+i] *= M2_infty;
//                fh_uh[(1+nd)*sz2 + 0*ng+i] *= M2_infty;
               fh_uh[0*sz2 + 1*ng+i] *= M2_infty;
               fh_uh[(1+nd)*sz2 + 1*ng+i] *= M2_infty;
               for (j = 0; j < nd; j++) {
//                    fh_uh[(1+j) * sz2 + 0 * ng + i] *= M2_infty;
                   fh_uh[(1+j) * sz2 + 1 * ng + i] *= M2_infty;
               }

               fh_uh[1 * sz2 + 2*ng+i] = sin(alpha_target);
               fh_uh[2 * sz2 + 2*ng+i] = - cos(alpha_target);

               if (nd == 3) {
                   fh_uh[1 * sz2 + 3*ng + i] = (uhg[1*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
                   fh_uh[2 * sz2 + 3*ng + i] = (uhg[2*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
                   fh_uh[3 * sz2 + 3*ng + i] = (uhg[3*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target) - 1.0;
               }

               aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
               aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
               for (j=0; j<nd; j++) {
                   aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
               }

               fh_uh[0*sz2 + (1+nd)*ng+i] = - (- uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[0*sz2 + (1+nd)*ng+i] = (- udg[(1+nd)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);

               for (j=0; j<nd; j++) {
                   fh_uh[0*sz2 + (1+nd)*ng+i] += - (uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
                   fh_u[0*sz2 + (1+nd)*ng+i] += (udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
                   fh_uh[0*sz2 + (1+nd)*ng+i] += nl[j*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_u[0*sz2 + (1+nd)*ng+i] += - nl[j*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
                   fh_uh[(1+j) * sz2 + (1+nd)*ng+i] += - nl[j*ng+i]/uhg[0*ng+i] + (uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
                   fh_u[(1+j) * sz2 + (1+nd)*ng+i] += nl[j*ng+i]/udg[0*ng+i] - (udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
               }

               fh_uh[(1+nd)*sz2 + (1+nd)*ng+i] = - (1.0/uhg[0*ng+i])*sqrt((gam/(gam-1.0))/aux1);
               fh_u[(1+nd)*sz2 + (1+nd)*ng+i] = (1.0/udg[0*ng+i])*sqrt((gam/(gam-1.0))/aux2);
           
               for (j = (2+nd); j < ncu; j++) {
                    fh_uh[j * sz2 + j * ng + i] = -1.0;
               }
           }
       }
        delete[] dqdu; delete[] dpdu; delete[] dM2du;
        delete[] dp02pdu; delete[] dp0du;
   }
   else if (ib==9) { // Subsonic turbomachinery outflow: Impose p. Extrapolate R+ = u + 2c/(gamma-1), s and u_t
       // TODO: Code ALE.

       double Vn, p_target, kinEnergy, kinEnergyHat, aux1, aux2;
       
       p_target = ui[0];

       for (i = 0; i < ng; i++) {
           kinEnergyHat = 0.0;
           kinEnergy = 0.0;
           for (j=0; j<nd; j++) {
               kinEnergyHat += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               kinEnergy += 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/udg[0*ng+i];
           }
           fh[0*ng+i] = M2_infty * ((gam-1.0) * (uhg[(1+nd)*ng+i]-kinEnergyHat) - p_target);

           // Extrapolate R+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: R+ = u_n + 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[1*ng+i] = 2.0*sqrt((gam/(gam-1.0))*aux1) - 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[1*ng+i] += nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate R+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: R+ = u_n + sqrt(gam*e/(gam-1))
           }

           // Extrapolate s = p / rho^gam [Dimensionless entropy: s = p / rho^gam , since we define s_ref := u_ref^2 / rho_ref^(gam-1)]
           // Note: p = (gam-1)*(rE - 0.5*r*u^2 - 0.5*r*v^2 - 0.5*r*w^2) [both dimensional and non-dimensional]
           fh[2*ng+i] = M2_infty * (uhg[(1+nd)*ng+i]-kinEnergyHat) / pow(uhg[0*ng+i],gam) - M2_infty * (udg[(1+nd)*ng+i]-kinEnergy) / pow(udg[0*ng+i],gam);

           // Extrapolate v_t:
           if (nd == 2) {
               double tx, ty, t_norm;
               tx = -nl[1*ng+i];
               ty = nl[0*ng+i];
               t_norm = sqrt(tx*tx+ty*ty);
               tx = tx / t_norm;
               ty = ty / t_norm;
               
               fh[3*ng+i] = tx*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                            ty*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]);
           }
           else if (nd == 3) {
               double t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, t2_norm, nl_norm;
               nl_norm = sqrt(nl[0*ng+i]*nl[0*ng+i] + nl[1*ng+i]*nl[1*ng+i] + nl[2*ng+i]*nl[2*ng+i]);
               if (abs(nl[0*ng+i])/nl_norm < 0.9) {
                   t1x = 0.0;
                   t1y = -nl[2*ng+i];
                   t1z = nl[1*ng+i];
               }
               else {
                   t1x = -nl[2*ng+i];
                   t1y = 0.0;
                   t1z = nl[0*ng+i];
               }
               t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
               t1x = t1x / t1_norm;
               t1y = t1y / t1_norm;
               t1z = t1z / t1_norm;
               
               t2x = t1z*nl[1*ng+i] - t1y*nl[2*ng+i];
               t2y = t1x*nl[2*ng+i] - t1z*nl[0*ng+i];
               t2z = t1y*nl[0*ng+i] - t1x*nl[1*ng+i];
               t2_norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
               t2x = t2x / t2_norm;
               t2y = t2y / t2_norm;
               t2z = t2z / t2_norm;

               fh[3*ng+i] = t1x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                            t1y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                            t1z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);

               fh[4*ng+i] = t2x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                            t2y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                            t2z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);
           }
           
           // Extrapolate turbulent variables (if necessary):
           for (j = (2+nd); j < ncu; j++) {
               fh[j*ng+i] = udg[i + j * ng] - uhg[i + j * ng];
           }
       }

       if (computeJacobian == 1) {
           for (i = 0; i < ng * nch * nc; i++)
               fh_u[i] = 0.0;

           for (i = 0; i < ng * nch * nch; i++)
               fh_uh[i] = 0.0;

           for (i = 0; i < ng; i++) {
               kinEnergyHat = 0.0;
               kinEnergy = 0.0;
               for (j=0; j<nd; j++) {
                   kinEnergyHat += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
                   kinEnergy += 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/udg[0*ng+i];

                   fh_uh[0*sz2     + 0*ng+i] += M2_infty * 0.5*(gam-1.0)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[(1+j)*sz2 + 0*ng+i] = - M2_infty * (gam-1.0)*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               }
               fh_uh[(1+nd)*sz2 + 0*ng+i] = M2_infty * (gam-1.0);


               aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
               aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
               for (j=0; j<nd; j++) {
                   aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
               }

               fh_uh[0*sz2 + 1*ng+i] = (- uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[0*sz2 + 1*ng+i] = - (- udg[(1+nd)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);

               for (j=0; j<nd; j++) {
                   fh_uh[0*sz2 + 1*ng+i] += (uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
                   fh_u[0*sz2 + 1*ng+i] += - (udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
                   fh_uh[0*sz2 + 1*ng+i] += - nl[j*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_u[0*sz2 + 1*ng+i] += nl[j*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
                   fh_uh[(1+j) * sz2 + 1*ng+i] += nl[j*ng+i]/uhg[0*ng+i] - (uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
                   fh_u[(1+j) * sz2 + 1*ng+i] += - nl[j*ng+i]/udg[0*ng+i] + (udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
               }

               fh_uh[(1+nd)*sz2 + 1*ng+i] = (1.0/uhg[0*ng+i])*sqrt((gam/(gam-1.0))/aux1);
               fh_u[(1+nd)*sz2 + 1*ng+i] = - (1.0/udg[0*ng+i])*sqrt((gam/(gam-1.0))/aux2);


               fh_uh[0*sz2 + 2*ng+i] = - M2_infty * gam * uhg[(1+nd)*ng+i] / pow(uhg[0*ng+i],gam+1.0) +
                                         M2_infty * (gam+1.0)*kinEnergyHat / pow(uhg[0*ng+i],gam+1.0);
               fh_u[0*sz2 + 2*ng+i] = M2_infty * gam * udg[(1+nd)*ng+i] / pow(udg[0*ng+i],gam+1.0) -
                                      M2_infty * (gam+1.0)*kinEnergy / pow(udg[0*ng+i],gam+1.0);

               for (k=0; k<nd; k++) {
                   fh_uh[(1+k)*sz2 + 2*ng+i] = - M2_infty * uhg[(1+k)*ng+i] / pow(uhg[0*ng+i],gam+1.0);
                   fh_u[(1+k)*sz2 + 2*ng+i] = M2_infty * udg[(1+k)*ng+i] / pow(udg[0*ng+i],gam+1.0);
               }

               fh_uh[(1+nd)*sz2 + 2*ng+i] = M2_infty / pow(uhg[0*ng+i],gam);
               fh_u[(1+nd)*sz2 + 2*ng+i] = - M2_infty / pow(udg[0*ng+i],gam);


               if (nd == 2) {
                   double tx, ty, t_norm;
                   tx = -nl[1*ng+i];
                   ty = nl[0*ng+i];
                   t_norm = sqrt(tx*tx+ty*ty);
                   tx = tx / t_norm;
                   ty = ty / t_norm;

                   fh_uh[0*sz2 + 3*ng+i] = - tx*uhg[1*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) +
                                           - ty*uhg[2*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[1*sz2 + 3*ng+i] = tx/uhg[0*ng+i];
                   fh_uh[2*sz2 + 3*ng+i] = ty/uhg[0*ng+i];

                   fh_u[0*sz2 + 3*ng+i] = tx*udg[1*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                          ty*udg[2*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
                   fh_u[1*sz2 + 3*ng+i] = - tx/udg[0*ng+i];
                   fh_u[2*sz2 + 3*ng+i] = - ty/udg[0*ng+i];
               }
               else if (nd == 3) {
                   double t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, t2_norm, nl_norm;
                   nl_norm = sqrt(nl[0*ng+i]*nl[0*ng+i] + nl[1*ng+i]*nl[1*ng+i] + nl[2*ng+i]*nl[2*ng+i]);
                   if (abs(nl[0*ng+i])/nl_norm < 0.9) {
                       t1x = 0.0;
                       t1y = -nl[2*ng+i];
                       t1z = nl[1*ng+i];
                   }
                   else {
                       t1x = -nl[2*ng+i];
                       t1y = 0.0;
                       t1z = nl[0*ng+i];
                   }
                   t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
                   t1x = t1x / t1_norm;
                   t1y = t1y / t1_norm;
                   t1z = t1z / t1_norm;

                   t2x = t1z*nl[1*ng+i] - t1y*nl[2*ng+i];
                   t2y = t1x*nl[2*ng+i] - t1z*nl[0*ng+i];
                   t2z = t1y*nl[0*ng+i] - t1x*nl[1*ng+i];
                   t2_norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
                   t2x = t2x / t2_norm;
                   t2y = t2y / t2_norm;
                   t2z = t2z / t2_norm;

                   fh_uh[0*sz2 + 3*ng+i] = - t1x*uhg[1*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                             t1y*uhg[2*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                             t1z*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[1*sz2 + 3*ng+i] = t1x/uhg[0*ng+i];
                   fh_uh[2*sz2 + 3*ng+i] = t1y/uhg[0*ng+i];
                   fh_uh[3*sz2 + 3*ng+i] = t1z/uhg[0*ng+i];

                   fh_u[0*sz2 + 3*ng+i] = t1x*udg[1*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                          t1y*udg[2*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                          t1z*udg[3*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
                   fh_u[1*sz2 + 3*ng+i] = - t1x/udg[0*ng+i];
                   fh_u[2*sz2 + 3*ng+i] = - t1y/udg[0*ng+i];
                   fh_u[3*sz2 + 3*ng+i] = - t1z/udg[0*ng+i];

                   fh_uh[0*sz2 + 4*ng+i] = - t2x*uhg[1*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                             t2y*uhg[2*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                             t2z*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[1*sz2 + 4*ng+i] = t2x/uhg[0*ng+i];
                   fh_uh[2*sz2 + 4*ng+i] = t2y/uhg[0*ng+i];
                   fh_uh[3*sz2 + 4*ng+i] = t2z/uhg[0*ng+i];

                   fh_u[0*sz2 + 4*ng+i] = t2x*udg[1*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                          t2y*udg[2*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                          t2z*udg[3*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
                   fh_u[1*sz2 + 4*ng+i] = - t2x/udg[0*ng+i];
                   fh_u[2*sz2 + 4*ng+i] = - t2y/udg[0*ng+i];
                   fh_u[3*sz2 + 4*ng+i] = - t2z/udg[0*ng+i];
               }
               
               for (j = (2+nd); j < ncu; j++) {
                    fh_u[j * sz2 + j * ng + i] = 1.0;
                    fh_uh[j * sz2 + j * ng + i] = -1.0;
               }
           }
       }
   }
   else if (ib==10) { // Subsonic turbomachinery inlet - Cenaero version: Impose p0, T0, alpha, theta. Extrapolate ||v||
       // TODO: Code ALE.
        double Vn, p0_target, T0_target, alpha_target, theta_target, uh_vel, u_vel;
        double q, p, M2, p02p, p0;
        double * dqdu = new double[2+nd];
        double * dpdu = new double[2+nd];
        double * dM2du = new double[2+nd];
        double * dp02pdu = new double[2+nd];
        double * dp0du = new double[2+nd];
        
        p0_target = ui[0];
        T0_target = ui[1];
        alpha_target = ui[2];
        if (nd == 3) {
            theta_target = ui[3];
        }
        
        for (i = 0; i < ng; i++) {
            // Impose p0:
//             fh[0*ng+i] = (gam-1.0)*uhg[(1+nd)*ng+i] - p0_target;
            //////
            q = 0.0;
            for (j = 0; j < nd; j++)
                q += 0.5 * uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i] / uhg[0*ng+i];
            p = gam1*(uhg[(1+nd)*ng+i] - q);
            M2 = (2.0*q) / (gam*p);
            p02p = pow(1.0 + 0.5*gam1*M2, gam/gam1);
            p0 = p*p02p;
            fh[0*ng+i] = p0 - p0_target;
            fh[0*ng+i] *= M2_infty;
            /////
           
           // Impose T0:
           fh[1*ng+i] = uhg[(1+nd)*ng+i]/uhg[0*ng+i] - T0_target;
           for (j = 0; j < nd; j++) {
//                fh[0*ng+i] += 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               fh[1*ng+i] -= 0.5*((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           }
//            fh[0*ng+i] *= M2_infty;
           fh[1*ng+i] *= M2_infty;

           // Impose alpha:
           fh[2*ng+i] = uhg[1*ng+i]*sin(alpha_target) - uhg[2*ng+i]*cos(alpha_target);

           // Impose theta:
           if (nd == 3)
               fh[3*ng+i] = sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i])*cos(theta_target) - uhg[3*ng+i];
            
           // Extrapolate ||v||:
           uh_vel = 0.0;
           u_vel = 0.0;
           for (j = 0; j < nd; j++) {
               uh_vel += uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i];
               u_vel += udg[(1+j)*ng+i]*udg[(1+j)*ng+i];
           }
           uh_vel = sqrt(uh_vel) / uhg[0*ng+i];
           u_vel = sqrt(u_vel) / udg[0*ng+i];
           
           fh[(1+nd)*ng+i] = uh_vel - u_vel;

           // Impose turbulent variables (if necessary):
           for (j = (2+nd); j < ncu; j++) {
               fh[j*ng+i] = ui[j] - uhg[i + j * ng];
           }
       }

       if (computeJacobian == 1) {
           for (i = 0; i < ng * nch * nc; i++)
               fh_u[i] = 0.0;

           for (i = 0; i < ng * nch * nch; i++)
               fh_uh[i] = 0.0;

           for (i = 0; i < ng; i++) {
              //////
                q = 0.0;
                for (j = 0; j < nd; j++)
                    q += 0.5 * uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i] / uhg[0*ng+i];
                dqdu[0] = - q / uhg[0*ng+i];
                for (j = 0; j < nd; j++)
                    dqdu[1+j] = uhg[(1+j)*ng+i] / uhg[0*ng+i];
                dqdu[1+nd] = 0.0;

                p = gam1*(uhg[(1+nd)*ng+i] - q);
                for (j = 0; j < (1+nd); j++)
                    dpdu[j] = - gam1*dqdu[j];
                dpdu[1+nd] = gam1;

                M2 = (2.0*q) / (gam*p);
                for (j = 0; j < (2+nd); j++)
                    dM2du[j] = (2.0 / gam) * ((dqdu[j]*p - q*dpdu[j]) / (p*p));
                
                p02p = pow(1.0 + 0.5*gam1*M2, gam/gam1);
                for (j = 0; j < (2+nd); j++)
                    dp02pdu[j] = 0.5 * gam * pow(1.0 + 0.5*gam1*M2, 1.0/gam1) * dM2du[j];

                p0 = p*p02p;
                for (j = 0; j < (2+nd); j++)
                    dp0du[j] = dpdu[j] * p02p + p * dp02pdu[j];
               
               for (j = 0; j < (2+nd); j++)
                   fh_uh[j * sz2 + 0 * ng + i] = M2_infty * dp0du[j];
               /////
               
//                fh_uh[(1+nd) * sz2 + 0 * ng + i] = gam-1.0;
               fh_uh[(1+nd) * sz2 + 1 * ng + i] = 1.0/uhg[0*ng+i];
               fh_uh[0 * sz2 + 1 * ng + i] = - uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               for (j = 0; j < nd; j++) {
//                    fh_uh[0 * sz2 + 0 * ng + i] -= 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
//                    fh_uh[(1+j) * sz2 + 0 * ng + i] += (2.0-gam)*uhg[(1+j)*ng+i]/uhg[0*ng+i];
                   fh_uh[0 * sz2 + 1 * ng + i] += ((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_uh[(1+j) * sz2 + 1 * ng + i] -= ((gam-1.0)/gam)*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               }
//                fh_uh[0*sz2 + 0*ng+i] *= M2_infty;
//                fh_uh[(1+nd)*sz2 + 0*ng+i] *= M2_infty;
               fh_uh[0*sz2 + 1*ng+i] *= M2_infty;
               fh_uh[(1+nd)*sz2 + 1*ng+i] *= M2_infty;
               for (j = 0; j < nd; j++) {
//                    fh_uh[(1+j) * sz2 + 0 * ng + i] *= M2_infty;
                   fh_uh[(1+j) * sz2 + 1 * ng + i] *= M2_infty;
               }

               fh_uh[1 * sz2 + 2*ng+i] = sin(alpha_target);
               fh_uh[2 * sz2 + 2*ng+i] = - cos(alpha_target);

               if (nd == 3) {
                   fh_uh[1 * sz2 + 3*ng + i] = (uhg[1*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
                   fh_uh[2 * sz2 + 3*ng + i] = (uhg[2*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
                   fh_uh[3 * sz2 + 3*ng + i] = (uhg[3*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target) - 1.0;
               }

               uh_vel = 0.0;
               u_vel = 0.0;
               for (j = 0; j < nd; j++) {
                   uh_vel += uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i];
                   u_vel += udg[(1+j)*ng+i]*udg[(1+j)*ng+i];
               }
               uh_vel = sqrt(uh_vel) / uhg[0*ng+i];
               u_vel = sqrt(u_vel) / udg[0*ng+i];

               fh_uh[0*sz2 + (1+nd)*ng+i] = - uh_vel / uhg[0*ng+i];
               fh_u[0*sz2 + (1+nd)*ng+i] = u_vel / udg[0*ng+i];
               for (j = 0; j < nd; j++) {
                   fh_uh[(1+j)*sz2 + (1+nd)*ng+i] = uhg[(1+j)*ng+i] / (uh_vel*uhg[0*ng+i]*uhg[0*ng+i]);
                   fh_u[(1+j)*sz2 + (1+nd)*ng+i] = - udg[(1+j)*ng+i] / (u_vel*udg[0*ng+i]*udg[0*ng+i]);
               }
               
               for (j = (2+nd); j < ncu; j++) {
                    fh_uh[j * sz2 + j * ng + i] = -1.0;
               }
           }
       }
        delete[] dqdu; delete[] dpdu; delete[] dM2du;
        delete[] dp02pdu; delete[] dp0du;
    }
//     else if (ib==11) { // Freestream (regular far-field) in rotating reference framce
//         // The boundary conditions are imposed in the absolute (non-rotating) reference frame.
//         // The basis used for velocity is that of the rotating reference frame.
//        
//         double *an = new double[ng * nch * nch];
//         double *An = new double[ng * nch * nch];
//         double *anm = new double[ng * nch * nch * nch];
//         double *Anm = new double[ng * nch * nch * nch];
//         
//         double *udg_abs = new double[ng * nc];      // Note: Although we allocate memory for all nc, only the first ncu components are required and actually computed.
//         double *uhg_abs = new double[ng * nch];
//         double *fh_u_abs = new double[ng * nch * nc];
//         double *fh_uh_abs = new double[ng * nch * nch];
//         double *ui_rotBasis = new double[ng * nch];
//         double *kinEn_rel = new double[ng];
//         double *kinEn_abs = new double[ng];
//         double *kinEnh_rel = new double[ng];
//         double *kinEnh_abs = new double[ng];
//         double *DkinEn_rel_Dug = new double[ng * ncu];
//         double *DkinEn_abs_Dug = new double[ng * ncu];
//         double *DkinEnh_rel_Duhg = new double[ng * nch];
//         double *DkinEnh_abs_Duhg = new double[ng * nch];
//         double *Dug_abs_Dug = new double[ng * ncu * ncu];
//         double *Duhg_abs_Duhg = new double[ng * nch * nch];
//         double *rotVel = new double[ng * nd];
//         double *CBM = new double[nd * nd];
//         
//         double Vn;
//         
//         // Compute change of basis matrix from non-rotating to rotating basis:
//         if (nd == 2) {
//             error("Need to define alpha\n");
// //             CBM[0*nd + 0] =   cos(alpha);
// //             CBM[0*nd + 1] = - sin(alpha);
// //             CBM[1*nd + 0] =   sin(alpha);
// //             CBM[1*nd + 1] =   cos(alpha);
//         }
//         else if (nd == 3) {
//             // Use Euler angles here.
//             error("Need to define Euler angles\n");
// //             CBM[0*nd + 0] = ;
// //             CBM[0*nd + 1] = ;
// //             CBM[0*nd + 2] = ;
// //             CBM[1*nd + 0] = ;
// //             CBM[1*nd + 1] = ;
// //             CBM[1*nd + 2] = ;
// //             CBM[2*nd + 0] = ;
// //             CBM[2*nd + 1] = ;
// //             CBM[2*nd + 2] = ;
//         }
//         
//         // Initialize arrays:
//         for (i = 0; i < ng*ncu*ncu; i++)
//             Dug_abs_Dug[i] = 0.0;
//         for (i = 0; i < ng*nch*nch; i++)
//             Duhg_abs_Duhg[i] = 0.0;
//         
//         for (i = 0; i < ng; i++) {
//             kinEn_rel[i] = 0.0; kinEn_abs[i] = 0.0;
//             kinEnh_rel[i] = 0.0; kinEnh_abs[i] = 0.0;
//         }
//         
//         for (i = 0; i < ng * ncu; i++) {
//             DkinEn_rel_Dug[i] = 0.0;
//             DkinEn_abs_Dug[i] = 0.0;
//         }
//         for (i = 0; i < ng * nch; i++) {
//             DkinEnh_rel_Duhg[i] = 0.0;
//             DkinEnh_abs_Duhg[i] = 0.0;
//         }
//         
//         for (i = 0; i < ng * nch; i++)
//             fh[i] = 0.0;
//         
//         if (computeJacobian == 1) {
//             for (i = 0; i < ng * nch * nc; i++) {
//                 fh_u[i] = 0.0;
//                 fh_u_abs[i] = 0.0;
//             }
// 
//             for (i = 0; i < ng * nch * nch; i++) {
//                 fh_uh[i] = 0.0;
//                 fh_uh_abs[i] = 0.0;
//             }
//         }
//         
//         for (i = 0; i < ng; i++) {
//             // Compute v_abs - v_rel, expressed in the basis of the rotating frame
//             if (nd == 2) {
//                 error("Need to define W_z\n");
// //                 rotVel[0*nd + i] = - W_z*pg[1 * ng + i];
// //                 rotVel[1*nd + i] =   W_z*pg[0 * ng + i];
//             }
//             else if (nd == 3) {
//                 error("Need to define W_x, W_y, W_z\n");
// //                 rotVel[0*nd + i] = W_y*pg[2 * ng + i] - W_z*pg[1 * ng + i];
// //                 rotVel[1*nd + i] = W_z*pg[0 * ng + i] - W_x*pg[2 * ng + i];
// //                 rotVel[2*nd + i] = W_x*pg[1 * ng + i] - W_y*pg[0 * ng + i];
//             }
//             
//             // Compute ui_abs, with velocities expressed in the basis of the rotating frame:
//             ui_rotBasis[0*ng + i] = ui[0];
//             for (j = 0; j < nd; j++) {
//                 ui_rotBasis[(1+j)*ng + i] = 0.0;
//                 for (int jj = 0; jj < nd; jj++)
//                     ui_rotBasis[(1+j)*ng + i] += CBM[jj*nd + j] * ui[(1+jj)*ng + i];
//             }
//             ui_rotBasis[(1+nd)*ng + i] = ui[1+nd];
//             for (j = nd+2; j < nch; j++)
//                 ui_rotBasis[j*ng + i] = ui[j];
//             
//             // Compute udg_abs, uhg_abs
//             // kinEn_rel, kinEn_abs, kinEnh_rel, kinEnh_abs,
//             // DkinEn_rel_Dug, DkinEn_abs_Dug, DkinEnh_rel_Duhg, DkinEnh_abs_Duhg
//             udg_abs[0*ng + i] = udg[0*ng + i];
//             uhg_abs[0*ng + i] = uhg[0*ng + i];
//             for (j = 0; j < nd; j++) {
//                 udg_abs[(1+j)*ng + i] = udg[(1+j)*ng + i] + udg[0*ng + i] * rotVel[j*ng + i];
//                 uhg_abs[(1+j)*ng + i] = uhg[(1+j)*ng + i] + uhg[0*ng + i] * rotVel[j*ng + i];
//                 kinEn_rel[i] += udg[(1+j)*ng + i] * udg[(1+j)*ng + i];
//                 kinEn_abs[i] += udg_abs[(1+j)*ng + i] * udg_abs[(1+j)*ng + i];
//                 kinEnh_rel[i] += uhg[(1+j)*ng + i] * uhg[(1+j)*ng + i];
//                 kinEnh_abs[i] += uhg_abs[(1+j)*ng + i] * uhg_abs[(1+j)*ng + i];
//             }
//             kinEn_rel[i] *= 0.5 / udg[0*ng + i];
//             kinEn_abs[i] *= 0.5 / udg_abs[0*ng + i];
//             kinEnh_rel[i] *= 0.5 / uhg[0*ng + i];
//             kinEnh_abs[i] *= 0.5 / uhg_abs[0*ng + i];
//             DkinEn_rel_Dug[0*ng + i] = - kinEn_rel[i] / udg[0*ng + i];
//             DkinEn_abs_Dug[0*ng + i] = - kinEn_abs[i] / udg_abs[0*ng + i];
//             DkinEnh_rel_Duhg[0*ng + i] = - kinEnh_rel[i] / uhg[0*ng + i];
//             DkinEnh_abs_Duhg[0*ng + i] = - kinEnh_abs[i] / uhg_abs[0*ng + i];
//             for (j = 0; j < nd; j++) {
//                 DkinEn_rel_Dug[(1+j)*ng + i] = udg[(1+j)*ng + i] / udg[0*ng + i];
//                 DkinEn_abs_Dug[0*ng + i] += udg_abs[(1+j)*ng + i] * rotVel[j*ng + i] / udg_abs[0*ng + i];
//                 DkinEn_abs_Dug[(1+j)*ng + i] = udg_abs[(1+j)*ng + i] / udg_abs[0*ng + i];
//                 
//                 DkinEnh_rel_Duhg[(1+j)*ng + i] = uhg[(1+j)*ng + i] / uhg[0*ng + i];
//                 DkinEnh_abs_Duhg[0*ng + i] += uhg_abs[(1+j)*ng + i] * rotVel[j*ng + i] / uhg_abs[0*ng + i];
//                 DkinEnh_abs_Duhg[(1+j)*ng + i] = uhg_abs[(1+j)*ng + i] / uhg_abs[0*ng + i];
//             }
//             udg_abs[(1+nd)*ng + i] = udg[(1+nd)*ng + i] - kinEn_rel[i] + kinEn_abs[i];
//             uhg_abs[(1+nd)*ng + i] = uhg[(1+nd)*ng + i] - kinEnh_rel[i] + kinEnh_abs[i];
//             for (j = 2+nd; j < ncu; j++)
//                 udg_abs[j*ng + i] = udg[j*ng + i];
//             for (j = 2+nd; j < nch; j++)
//                 uhg_abs[j*ng + i] = uhg[j*ng + i];
//             
//             // Compute Dug_abs_Dug, Duhg_abs_Duhg:
//             Dug_abs_Dug[0*ncu*ng + 0*ng + i] = 1.0;
//             Duhg_abs_Duhg[0*nch*ng + 0*ng + i] = 1.0;
//             for (j = 0; j < nd; j++) {
//                 Dug_abs_Dug[0*ncu*ng + (1+j)*ng + i] = rotVel[j*ng + i];
//                 Dug_abs_Dug[(1+j)*ncu*ng + (1+j)*ng + i] = 1.0;
//                 Duhg_abs_Duhg[0*nch*ng + (1+j)*ng + i] = rotVel[j*ng + i];
//                 Duhg_abs_Duhg[(1+j)*nch*ng + (1+j)*ng + i] = 1.0;
//             }
//             Dug_abs_Dug[(1+nd)*ncu*ng + (1+nd)*ng + i] = 1.0;
//             Duhg_abs_Duhg[(1+nd)*nch*ng + (1+nd)*ng + i] = 1.0;
//             Dug_abs_Dug[0*ncu*ng + (1+nd)*ng + i] = - DkinEn_rel_Dug[0*ng + i] + DkinEn_abs_Dug[0*ng + i];
//             Duhg_abs_Duhg[0*nch*ng + (1+nd)*ng + i] = - DkinEnh_rel_Duhg[0*ng + i] + DkinEnh_abs_Duhg[0*ng + i];
//             for (j = 0; j < nd; j++) {
//                 Dug_abs_Dug[(1+j)*ncu*ng + (1+nd)*ng + i] = - DkinEn_rel_Dug[(1+j)*ng + i] + DkinEn_abs_Dug[(1+j)*ng + i];
//                 Duhg_abs_Duhg[(1+j)*nch*ng + (1+nd)*ng + i] = - DkinEnh_rel_Duhg[(1+j)*ng + i] + DkinEnh_abs_Duhg[(1+j)*ng + i];
//             }
//             for (j = 2+nd; j < ncu; j++)
//                 Dug_abs_Dug[j*ncu*ng + j*ng + i] = 1.0;
//             for (j = 2+nd; j < nch; j++)
//                 Duhg_abs_Duhg[j*nch*ng + j*ng + i] = 1.0;
//         }
//         
//         getanNEW(an, anm, uhg_abs, pg, nl, app, param, 0, ng, nd, computeJacobian);
//         getanNEW(An, Anm, uhg_abs, pg, nl, app, param, 1, ng, nd, computeJacobian);
//         
//         // Boundary conditions for mass, momentum and energy equations:
//         for (k = 0; k < (2 + nd); k++)
//             for (j = 0; j < (2 + nd); j++)
//                 for (i = 0; i < ng; i++) {
//                     nm = k * ng * (nd + 2) + j * ng + i;
//                     nm2 = k * sz2 + j * ng + i;
//                     nk = k * ng + i;
//                     fh[j * ng + i] += (an[nm] + An[nm]) * (udg_abs[nk]      - uhg_abs[nk]) -
//                                       (an[nm] - An[nm]) * (ui_rotBasis[nk]  - uhg_abs[nk]);
// 
//                     if (computeJacobian == 1) {
//                         fh_u_abs[nm2] = an[nm] + An[nm];
//                         fh_uh_abs[nm2] = -2.0 * An[nm];
//                         for (n = 0; n < (nd + 2); n++) {
//                             nk = k * (nd + 2) * (nd + 2) * ng + n * (nd + 2) * ng + j * ng + i;
//                             nn = n * ng + i;
//                             fh_uh_abs[nm2] += (anm[nk] + Anm[nk]) * (udg_abs[nn]       - uhg_abs[nn]) -
//                                               (anm[nk] - Anm[nk]) * (ui_rotBasis[nn]   - uhg_abs[nn]);  
//                         }
//                     }
//                 }
//         
//         if (computeJacobian == 1) {
//             for (j = 0; j < (2 + nd); j++)
//                 for (i = 0; i < ng; i++) {
//                     fh_u[0 * sz2 + j * ng + i] += fh_u_abs[0 * sz2 + j * ng + i];
//                     fh_u[0 * sz2 + j * ng + i] += fh_u_abs[(1+nd) * sz2 + j * ng + i] * Dug_abs_Dug[0*ncu*ng + (1+nd)*ng + i];
//                     for (k = 0; k < nd; k++) {
//                         fh_u[0 * sz2 + j * ng + i] += fh_u_abs[(1+k) * sz2 + j * ng + i] * Dug_abs_Dug[0*ncu*ng + (1+k)*ng + i];
//                         fh_u[(1+k) * sz2 + j * ng + i] += fh_u_abs[(1+k) * sz2 + j * ng + i] * Dug_abs_Dug[(1+k)*ncu*ng + (1+k)*ng + i];
//                         fh_u[(1+k) * sz2 + j * ng + i] += fh_u_abs[(1+nd) * sz2 + j * ng + i] * Dug_abs_Dug[(1+k)*ncu*ng + (1+nd)*ng + i];
//                     }
//                     fh_u[(1+nd) * sz2 + j * ng + i] += fh_u_abs[(1+nd) * sz2 + j * ng + i];
//                     
//                     fh_uh[0 * sz2 + j * ng + i] += fh_uh_abs[0 * sz2 + j * ng + i];
//                     fh_uh[0 * sz2 + j * ng + i] += fh_uh_abs[(1+nd) * sz2 + j * ng + i] * Duhg_abs_Duhg[0*nch*ng + (1+nd)*ng + i];
//                     for (k = 0; k < nd; k++) {
//                         fh_uh[0 * sz2 + j * ng + i] += fh_uh_abs[(1+k) * sz2 + j * ng + i] * Duhg_abs_Duhg[0*nch*ng + (1+k)*ng + i];
//                         fh_uh[(1+k) * sz2 + j * ng + i] += fh_uh_abs[(1+k) * sz2 + j * ng + i] * Duhg_abs_Duhg[(1+k)*nch*ng + (1+k)*ng + i];
//                         fh_uh[(1+k) * sz2 + j * ng + i] += fh_uh_abs[(1+nd) * sz2 + j * ng + i] * Duhg_abs_Duhg[(1+k)*nch*ng + (1+nd)*ng + i];
//                     }
//                     fh_uh[(1+nd) * sz2 + j * ng + i] += fh_uh_abs[(1+nd) * sz2 + j * ng + i];
//                 }
//         }
//         
//         // Boundary conditions for the SA equation
//         for (j = 2+nd; j < nch; j++)
//             for (i = 0; i < ng; i++) {
//                 // Compute normal component of linear momentum
//                 error("Error No. 346H8 in fbou.cpp. Not sure about how to treat Vn (relative or absolute?) here.\n");
//                 Vn = 0.0;
//                 for (ii = 0; ii < nd; ii++)
//                     Vn += nl[ii * ng + i] * uhg[i + (1 + ii) * ng];
//                 
//                 if (Vn < 0.0) {       // Inflow (Dirichlet for rN)
//                     fh[j * ng + i] = ui_rotBasis[i + j * ng] - uhg_abs[i + j * ng];
//                     if (computeJacobian == 1)
//                         fh_uh[j * sz2 + j * ng + i] = -1.0;
//                 }
//                 else {                // Outflow (extrapolation of rN)
//                     fh[j * ng + i] = udg_abs[i + j * ng] - uhg_abs[i + j * ng];
//                     if (computeJacobian == 1) {
//                         fh_u[j * sz2 + j * ng + i] = 1.0;
//                         fh_uh[j * sz2 + j * ng + i] = -1.0;
//                     }
//                 }
//             }
//         
//         delete[] an; delete[] An; delete[] anm; delete[] Anm;
//         delete[] udg_abs; delete[] uhg_abs;
//         delete[] ui_rotBasis;
//         delete[] fh_u_abs; delete[] fh_uh_abs;
//         delete[] kinEn_rel; delete[] kinEn_abs;
//         delete[] DkinEn_rel_Dug; delete[] DkinEn_abs_Dug;
//         delete[] kinEnh_rel; delete[] kinEnh_abs;
//         delete[] DkinEnh_rel_Duhg; delete[] DkinEnh_abs_Duhg;
//         delete[] Duhg_abs_Duhg; delete[] Dug_abs_Dug;
//         delete[] rotVel; delete[] CBM;
//     }
    else if (ib==11) { // Adiabatic, isothermal wall based on Roe's approximate Riemann solver
        
        double qRoe, qPlus, qMinus, numeratorH, denominatorH, numerator, denominator, rH_Roe;
        
        double *fPlus = new double[ng * nch];
        double *fPlus_u = new double[ng * nch * nc];
        double *fPlus_uh = new double[ng * nch * nch];
        double *fMinus = new double[ng * nch];
        double *fMinus_u = new double[ng * nch * nc];
        double *fMinus_uh = new double[ng * nch * nch];
        double *An = new double[ng * nch * nch];
        double *An_uRoe = new double[ng * nch * nch * nch];
        double *An_uPlus = new double[ng * nch * nch * nch];
        double *uRoe = new double[ng * ncu];
        double *uPlus = new double[ng * ncu];
        double *uMinus = new double[ng * ncu];
        double *DuRoe_DuPlus = new double[ng * ncu * ncu];
        double *DuRoe_DuMinus = new double[ng * ncu * ncu];
        double *DuMinus_DuPlus = new double[ng * ncu * ncu];
        double *DqDu_Plus = new double[ncu];
        double *DqDu_Minus = new double[ncu];
        double *DqDu_Roe = new double[ncu];
        double *DqRoe_DuPlus = new double[ncu];
        double *DqRoe_DuMinus = new double[ncu];
        double *Dnumerator_DuPlus = new double[ncu];
        double *Dnumerator_DuMinus = new double[ncu];
        double *Ddenominator_DuPlus = new double[ncu];
        double *Ddenominator_DuMinus = new double[ncu];
        double *DnumeratorH_DuPlus = new double[ncu];
        double *DnumeratorH_DuMinus = new double[ncu];
        double *DdenominatorH_DuPlus = new double[ncu];
        double *DdenominatorH_DuMinus = new double[ncu];
        double *DrH_Roe_DuPlus = new double[ncu];
        double *DrH_Roe_DuMinus = new double[ncu];
        
        double kinEnrgy, kinEnrgyInf, T_infty;
        
        kinEnrgyInf = 0.0;
        for (j=0; j<nd; j++)
            kinEnrgyInf += 0.5*ui[1+j]*ui[1+j]/ui[0];
        T_infty = (ui[1+nd] - kinEnrgyInf) / ui[0];
        
        double ratio = 500.0/200.0;     // 500.0/200.0: Value for cylinder at M = 18, Re = 378k in Barer's thesis
        double T_w = ratio * T_infty;
        
        // Get Roe average state:
        for (i = 0; i < ng; i++) {
            // Plus state (i.e. interior state):
            for (j = 0; j < ncu; j++)
                uPlus[j*ng+i] = udg[j*ng+i];
            
            // Minus state (i.e. boundary state):
            uMinus[0*ng+i] = uPlus[(1+nd)*ng+i] / T_w;
            for (j = 0; j < nd; j++)
                uMinus[(1+j)*ng+i] = (uPlus[(1+nd)*ng+i] / T_w) * Vg[j*ng+i];
            uMinus[(1+nd)*ng+i] = uPlus[(1+nd)*ng+i];
            for (j=(nd+2); j<ncu; j++)
                uMinus[j*ng+i] = 0.0;
            
            for (j = 0; j<ncu*ncu; j++)
                DuMinus_DuPlus[j*ng+i] = 0.0;
            DuMinus_DuPlus[(1+nd)*ncu*ng+0*ng+i] = 1.0 / T_w;
            for (j = 0; j < nd; j++)
                DuMinus_DuPlus[(1+nd)*ncu*ng+(1+j)*ng+i] = Vg[j*ng+i] / T_w;
            DuMinus_DuPlus[(1+nd)*ncu*ng+(1+nd)*ng+i] = 1.0;
            
            // Roe average state and derivatives with respect to plus and minus states:
            for (j = 0; j < ncu*ncu; j++) {
                DuRoe_DuPlus[j*ng+i] = 0.0;
                DuRoe_DuMinus[j*ng+i] = 0.0;
            }
            
            // - Density:
            uRoe[0*ng+i] = sqrt(uPlus[0*ng+i]) * sqrt(uMinus[0*ng+i]);
            DuRoe_DuPlus[0*ncu*ng+0*ng+i] = 0.5 * sqrt(uMinus[0*ng+i]) / sqrt(uPlus[0*ng+i]);
            DuRoe_DuMinus[0*ncu*ng+0*ng+i] = 0.5 * sqrt(uPlus[0*ng+i]) / sqrt(uMinus[0*ng+i]);
            
            // - Momentum:
            for (j = 0; j < nd; j++) {
                numerator = sqrt(uMinus[0*ng+i])*uPlus[(1+j)*ng+i] + sqrt(uPlus[0*ng+i])*uMinus[(1+j)*ng+i];
                denominator = sqrt(uMinus[0*ng+i]) + sqrt(uPlus[0*ng+i]);
                
                for (k = 0; k < ncu; k++) {
                    Dnumerator_DuPlus[k] = 0.0;
                    Dnumerator_DuMinus[k] = 0.0;
                    Ddenominator_DuPlus[k] = 0.0;
                    Ddenominator_DuMinus[k] = 0.0;
                }
                Dnumerator_DuPlus[0] = 0.5 * uMinus[(1+j)*ng+i] / sqrt(uPlus[0*ng+i]);
                Dnumerator_DuMinus[0] = 0.5 * uPlus[(1+j)*ng+i] / sqrt(uMinus[0*ng+i]);
                Dnumerator_DuPlus[1+j] = sqrt(uMinus[0*ng+i]);
                Dnumerator_DuMinus[1+j] = sqrt(uPlus[0*ng+i]);
                Ddenominator_DuPlus[0] = 0.5 / sqrt(uPlus[0*ng+i]);
                Ddenominator_DuMinus[0] = 0.5 / sqrt(uMinus[0*ng+i]);
                
                uRoe[(1+j)*ng+i] = numerator / denominator;
                
                for (k = 0; k < ncu; k++) {
                    DuRoe_DuPlus[k*ncu*ng+(1+j)*ng+i] = (Dnumerator_DuPlus[k] * denominator - Ddenominator_DuPlus[k] * numerator) / (denominator * denominator);
                    DuRoe_DuMinus[k*ncu*ng+(1+j)*ng+i] = (Dnumerator_DuMinus[k] * denominator - Ddenominator_DuMinus[k] * numerator) / (denominator * denominator);
                }
            }
            
            // - Kinetic energy:
            qPlus = 0.0; qMinus = 0.0; qRoe = 0.0;
            for (j = 0; j < nd; j++) {
                qPlus += 0.5 * uPlus[(1+j)*ng+i] * uPlus[(1+j)*ng+i] / uPlus[0*ng+i];
                qMinus += 0.5 * uMinus[(1+j)*ng+i] * uMinus[(1+j)*ng+i] / uMinus[0*ng+i];
                qRoe += 0.5 * uRoe[(1+j)*ng+i] * uRoe[(1+j)*ng+i] / uRoe[0*ng+i];
            }
            DqDu_Plus[0] = - qPlus / uPlus[0*ng+i];
            DqDu_Minus[0] = - qMinus / uMinus[0*ng+i];
            DqDu_Roe[0] = - qRoe / uRoe[0*ng+i];
            for (j = 0; j < nd; j++) {
                DqDu_Plus[1+j] = uPlus[(1+j)*ng+i] / uPlus[0*ng+i];
                DqDu_Minus[1+j] = uMinus[(1+j)*ng+i] / uMinus[0*ng+i];
                DqDu_Roe[1+j] = uRoe[(1+j)*ng+i] / uRoe[0*ng+i];
            }
            for (j = (1+nd); j < ncu; j++) {
                DqDu_Plus[j] = 0.0;
                DqDu_Minus[j] = 0.0;
                DqDu_Roe[j] = 0.0;
            }
            
            for (j = 0; j < (1+nd); j++) {
                DqRoe_DuPlus[j] = 0.0;
                DqRoe_DuMinus[j] = 0.0;
                for (k = 0; k < (1+nd); k++) {
                    DqRoe_DuPlus[j] += DqDu_Roe[k] * DuRoe_DuPlus[j*ncu*ng+k*ng+i];         // This line must be after DuRoe_DuPlus for density and momentum are computed
                    DqRoe_DuMinus[j] += DqDu_Roe[k] * DuRoe_DuMinus[j*ncu*ng+k*ng+i];       // This line must be after DuRoe_DuMinus for density and momentum are computed
                }
            }
            for (j = (1+nd); j < ncu; j++) {
                DqRoe_DuPlus[j] = 0.0;
                DqRoe_DuMinus[j] = 0.0;
            }
            
            // - Total enthalpy:
            numeratorH = sqrt(uMinus[0*ng+i])*gam*uPlus[(1+nd)*ng+i] - gam1*sqrt(uMinus[0*ng+i])*qPlus + sqrt(uPlus[0*ng+i])*gam*uMinus[(1+nd)*ng+i] - gam1*sqrt(uPlus[0*ng+i])*qMinus;
            denominatorH = sqrt(uMinus[0*ng+i]) + sqrt(uPlus[0*ng+i]);
            
            for (k=0; k<ncu; k++) {
                DnumeratorH_DuPlus[k] = 0.0;
                DnumeratorH_DuMinus[k] = 0.0;
                DdenominatorH_DuPlus[k] = 0.0;
                DdenominatorH_DuMinus[k] = 0.0;
            }
            DnumeratorH_DuPlus[0] = - gam1*sqrt(uMinus[0*ng+i])*DqDu_Plus[0] + 0.5 * ( gam*uMinus[(1+nd)*ng+i] - gam1*qMinus ) / sqrt(uPlus[0*ng+i]);
            DnumeratorH_DuMinus[0] = - gam1*sqrt(uPlus[0*ng+i])*DqDu_Minus[0] + 0.5 * ( gam*uPlus[(1+nd)*ng+i] - gam1*qPlus ) / sqrt(uMinus[0*ng+i]);
            for (k=0; k<nd; k++) {
                DnumeratorH_DuPlus[1+k] = - gam1*sqrt(uMinus[0*ng+i])*DqDu_Plus[1+k];
                DnumeratorH_DuMinus[1+k] = - gam1*sqrt(uPlus[0*ng+i])*DqDu_Minus[1+k];
            }
            DnumeratorH_DuPlus[1+nd] = sqrt(uMinus[0*ng+i])*gam - gam1*sqrt(uMinus[0*ng+i])*DqDu_Plus[1+nd];
            DnumeratorH_DuMinus[1+nd] = sqrt(uPlus[0*ng+i])*gam - gam1*sqrt(uPlus[0*ng+i])*DqDu_Minus[1+nd];
            DdenominatorH_DuPlus[0] = 0.5 / sqrt(uPlus[0*ng+i]);
            DdenominatorH_DuMinus[0] = 0.5 / sqrt(uMinus[0*ng+i]);
            
            rH_Roe = numeratorH / denominatorH;
            
            for (k=0; k<ncu; k++) {
                DrH_Roe_DuPlus[k] = (DnumeratorH_DuPlus[k] * denominatorH - DdenominatorH_DuPlus[k] * numeratorH) / (denominatorH * denominatorH);
                DrH_Roe_DuMinus[k] = (DnumeratorH_DuMinus[k] * denominatorH - DdenominatorH_DuMinus[k] * numeratorH) / (denominatorH * denominatorH);
            }
            
            // - Total energy:
            uRoe[(1+nd)*ng+i] = (rH_Roe + gam1*qRoe) / gam;
            for (k=0; k<ncu; k++) {
                DuRoe_DuPlus[k*ncu*ng+(1+nd)*ng+i] = (DrH_Roe_DuPlus[k] + gam1 * DqRoe_DuPlus[k]) / gam;
                DuRoe_DuMinus[k*ncu*ng+(1+nd)*ng+i] = (DrH_Roe_DuMinus[k] + gam1 * DqRoe_DuMinus[k]) / gam;
            }
            
            for (j=(nd+2); j<ncu; j++) {
                uRoe[j*ng+i] = 0.5 * (uPlus[j*ng+i] + uMinus[j*ng+i]);     // This is a made up state for the turbulent variables. It is irrelevant in practice since we do not use the Riemann solver for them.
                DuRoe_DuPlus[j*ncu*ng+j*ng+i] = 0.5;
                DuRoe_DuMinus[j*ncu*ng+j*ng+i] = 0.5;
            }
            
            // Update DuRoe_DuPlus from partial to total derivative:
            for (j=0; j<ncu; j++)
                for (k=0; k<ncu; k++)
                    for (l=0; l<ncu; l++)
                        DuRoe_DuPlus[k*ncu*ng+j*ng+i] += DuRoe_DuMinus[l*ncu*ng+j*ng+i] * DuMinus_DuPlus[k*ncu*ng+l*ng+i];
        }
        
        //getanNEW(An, An_uRoe, uRoe, pg, nl, app, param, 1, ng, nd, computeJacobian);
        getan(An, An_uRoe, uRoe, nl, param, 1, ng, nch, nd);       
        
        // Total derivatives of An with respect to uPlus:
        if (computeJacobian == 1) {
            for (l = 0; l < (2 + nd); l++)
                for (k = 0; k < (2 + nd); k++)
                    for (j = 0; j < (2 + nd); j++)
                        for (i = 0; i < ng; i++) {
                            An_uPlus[l*(nd+2)*(nd+2)*ng+k*(nd+2)*ng+j*ng+i] = 0.0;
                            for (n = 0; n < (2 + nd); n++)
                                An_uPlus[l*(nd+2)*(nd+2)*ng+k*(nd+2)*ng+j*ng+i] += An_uRoe[n*(nd+2)*(nd+2)*ng+k*(nd+2)*ng+j*ng+i] * DuRoe_DuPlus[l*ncu*ng+n*ng+i];
                        }
        }
        
        // Boundary conditions for mass, momentum and energy equations:
        Int fhatExpression = 0;
        fhatDriver(fh, fh_u, fh_uh, pg, udg, uhg, odg, nl, mesh, master, app, sol, temp, fhatExpression, ie, quadchoice, computeJacobian, numPoints);       // TODO: This won't work for ALE
        fhatExpression = 2;
        fhatDriver(fPlus, fPlus_u, fPlus_uh, pg, uPlus, uhg, odg, nl, mesh, master, app, sol, temp, fhatExpression, ie, quadchoice, computeJacobian, numPoints);       // TODO: This won't work for ALE
        fhatDriver(fMinus, fMinus_u, fMinus_uh, pg, uMinus, uhg, odg, nl,  mesh, master, app, sol, temp, fhatExpression, ie, quadchoice, computeJacobian, numPoints);       // TODO: This won't work for ALE
        for (j = 0; j < (2 + nd); j++)
            for (i = 0; i < ng; i++)
                fh[j*ng+i] -= 0.5 * (fPlus[j*ng+i] + fMinus[j*ng+i]);
        for (k = 0; k < (2 + nd); k++)
            for (j = 0; j < (2 + nd); j++)
                for (i = 0; i < ng; i++) {
                    nm = k * ng * (nd + 2) + j * ng + i;
                    nm2 = k * sz2 + j * ng + i;
                    nk = k * ng + i;
                    fh[j*ng+i] -= 0.5 * An[nm] * (uPlus[nk] - uMinus[nk]);
                    
                    if (computeJacobian == 1) {
                        fh_u[nm2] -= 0.5 * fPlus_u[nm2];
                        fh_u[nm2] -= 0.5 * An[nm];
                        for (n = 0; n < ncu; n++) {
                            fh_u[nm2] -= 0.5 * fMinus_u[n*sz2+j*ng+i] * DuMinus_DuPlus[k*sz2+n*ng+i];
                        }
                        for (n = 0; n < (nd + 2); n++) {
                            fh_u[nm2] -= 0.5 * An_uPlus[k*(nd+2)*(nd+2)*ng + n*(nd+2)*ng + j*ng + i] * (uPlus[n*ng + i] - uMinus[n*ng + i]);
                            fh_u[nm2] -= 0.5 * An[n*(nd+2)*ng+j*ng+i] * ( - DuMinus_DuPlus[k*sz2+n*ng+i] );
                        }
                    }
                }
        if (computeJacobian == 1) {
            for (k = (2+nd); k < ncu; k++)
                for (j = 0; j < (2 + nd); j++)
                    for (i = 0; i < ng; i++) {
                        nm = k * ng * (nd + 2) + j * ng + i;
                        nm2 = k * sz2 + j * ng + i;
                        nk = k * ng + i;
                            fh_u[nm2] -= 0.5 * fPlus_u[nm2];
                            for (n = 0; n < ncu; n++) {
                                fh_u[nm2] -= 0.5 * fMinus_u[n*sz2+j*ng+i] * DuMinus_DuPlus[k*sz2+n*ng+i];
                            }
                            for (n = 0; n < (nd + 2); n++) {
                                fh_u[nm2] -= 0.5 * An[n*(nd+2)*ng+j*ng+i] * ( - DuMinus_DuPlus[k*sz2+n*ng+i] );
                            }
                    }
        }
        
        // Boundary condition for SA equation (homogeneous Dirichlet):
        for (j=(nd+2); j<nch; j++)
            for (i=0; i<ng; i++) {
                fh[j*ng+i] = -uhg[j*ng+i];
                if (computeJacobian == 1) {
                    for (k=0; k<nch; k++)
                        fh_uh[k*nch*ng+j*ng+i] = 0.0;
                    fh_uh[j*nch*ng+j*ng+i] = -1.0;
                    for (k=0; k<nc; k++)
                        fh_u[k*nch*ng+j*ng+i] = 0.0;
                }
            }
        
        delete[] fPlus; delete[] fPlus_u; delete[] fPlus_uh;
        delete[] fMinus; delete[] fMinus_u; delete[] fMinus_uh;
        delete[] An; delete[] An_uRoe; delete[] An_uPlus;
        delete[] uRoe; delete[] uPlus; delete[] uMinus;
        delete[] DuRoe_DuPlus; delete[] DuRoe_DuMinus; delete[] DuMinus_DuPlus;
        delete[] DqDu_Plus; delete[] DqDu_Minus; delete[] DqDu_Roe;
        delete[] DqRoe_DuPlus; delete[] DqRoe_DuMinus;
        delete[] Dnumerator_DuPlus; delete[] Dnumerator_DuMinus;
        delete[] Ddenominator_DuPlus; delete[] Ddenominator_DuMinus;
        delete[] DnumeratorH_DuPlus; delete[] DnumeratorH_DuMinus;
        delete[] DdenominatorH_DuPlus; delete[] DdenominatorH_DuMinus;
        delete[] DrH_Roe_DuPlus; delete[] DrH_Roe_DuMinus;
    }
    else {
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.\n", ib);
        exit(-1);
    }
    
    if ((ib==2 || ib==3 || ib == 7 || ib == 11) && app.ALEflag == 0) {
        delete[] Vg;
    }
}

#endif
