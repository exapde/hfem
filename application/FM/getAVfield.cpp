#ifndef __GETAVFIELD
#define __GETAVFIELD

#include "../../utilities/DG_2_p1CG.cpp"

// Written by: C. Nguyen & P. Fernandez

void getNumTurbEqs(appstruct &app, Int *nte)
{
    switch (app.appname) {
        case 0:     // Euler
            *nte = 0;
            break;
        case 1:     // Navier-Stokes
            *nte = 0;
            break;
        case 3:     // RANS-SA
            *nte = 1;
            break;
        default: {
            printf("Application not implemented (appname = %d)\n",app.appname);
            exit(-1);
        }
    }
}

void checkAVfield(double *avField, Int *ndims)
{
    Int ne = ndims[5];
    Int npv = ndims[9];
    
    for (Int i = 0; i < ne*npv; i++) {
        if (isnan(avField[i])) {
            printf("\nWARNING: NaN in AV field. It will be made equal to zero.\n");
            avField[i] = 0.0;
        }
        else if (isinf(avField[i])) {
            printf("\nWARNING: Inf in AV field. It will be made equal to zero.\n");
            avField[i] = 0.0;
        }
    }
}

void applySmoothing(double *smoothField, double *nonSmoothField, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int numSmoothing, Int avgMethod)
{
    Int i, j;
    Int ne = ndims[5];
    Int npv = ndims[9];
    
    for (i = 0; i < numSmoothing; i++) {
        DG_2_p1CG(&smoothField[0], &nonSmoothField[0], mesh, master, app, ndims, avgMethod);
        
        if (i+1 < numSmoothing) {
            for (j = 0; j < ne*npv; j++)
                nonSmoothField[j] = smoothField[j];
        }
    }
}

void getAVfield_velDiv_2d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;

    double *param = &app.param[0];

    double alpha = 1.0e30;
    double beta = 1.0e-2;
    double k_h = 1.5;
    double Minf = param[4];
    double pi = 3.141592653589793;
    Int avgMethod = 1;

    double r, ru, rv, rE, rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u, v, ux, vy, h, f, x;

    double *x_p1CG, *x_DG, *f_p1CG, *f_DG;
    if (avgType == 0) {     // Sensor averaging
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 1) {        // Max surrogate averaging
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            rx = UDG[i*nc*npv+4*npv+j];
            rux = UDG[i*nc*npv+5*npv+j];
            rvx = UDG[i*nc*npv+6*npv+j];
            rEx = UDG[i*nc*npv+7*npv+j];
            ry = UDG[i*nc*npv+8*npv+j];
            ruy = UDG[i*nc*npv+9*npv+j];
            rvy = UDG[i*nc*npv+10*npv+j];
            rEy = UDG[i*nc*npv+11*npv+j];

            u = ru/r;
            v = rv/r;
            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;

            h = pn[i*npv*ncd + nd*npv + j];

            if (avgType == 0) {     // Sensor averaging
                x_DG[i*npv+j] = (k_h*h/ ((double) porder))*(ux+vy);
            }
            else if (avgType == 1) { // Max surrogate averaging
                x = (k_h*h/((double) porder))*(ux+vy);
                f_DG[i*npv+j] = (x-beta)*(atan(alpha*(x-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        DG_2_p1CG(x_p1CG, x_DG, mesh, master, app, ndims, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + nd*npv + j];
                f = (x_p1CG[i*npv+j]-beta)*(atan(alpha*(x_p1CG[i*npv+j]-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder))*(1.0+1.0/Minf)*f;
            }
        }
    }
    else if (avgType == 1) { // Max surrogate averaging
        DG_2_p1CG(f_p1CG, f_DG, mesh, master, app, ndims, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + nd*npv + j];
                avField[i*npv+j] = (k_h*h/((double) porder))*(1.0+1.0/Minf)*f_p1CG[i*npv+j];
            }
        }
    }
    
    if (avgType == 0) {     // Sensor averaging
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 1) {        // Max surrogate averaging
        delete[] f_p1CG; delete[] f_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_velDiv_3d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;

    double *param = &app.param[0];

    double alpha = 1.0e30;
    double beta = 1.0e-2;
    double k_h = 1.5;
    double Minf = param[4];
    double pi = 3.141592653589793;
    Int avgMethod = 1;

    double r, ru, rv, rw, rE, rx, rux, rvx, rwx, rEx, ry, ruy, rvy, rwy, rEy, rz, ruz, rvz, rwz, rEz, u, v, w, ux, vy, wz, h, f, x;

    double *x_p1CG, *x_DG, *f_p1CG, *f_DG;
    if (avgType == 0) {     // Sensor averaging
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 1) {        // Max surrogate averaging
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }

    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rw = UDG[i*nc*npv+3*npv+j];
            rE = UDG[i*nc*npv+4*npv+j];
            rx = UDG[i*nc*npv+5*npv+j];
            rux = UDG[i*nc*npv+6*npv+j];
            rvx = UDG[i*nc*npv+7*npv+j];
            rwx = UDG[i*nc*npv+8*npv+j];
            rEx = UDG[i*nc*npv+9*npv+j];
            ry = UDG[i*nc*npv+10*npv+j];
            ruy = UDG[i*nc*npv+11*npv+j];
            rvy = UDG[i*nc*npv+12*npv+j];
            rwy = UDG[i*nc*npv+13*npv+j];
            rEy = UDG[i*nc*npv+14*npv+j];
            rz = UDG[i*nc*npv+15*npv+j];
            ruz = UDG[i*nc*npv+16*npv+j];
            rvz = UDG[i*nc*npv+17*npv+j];
            rwz = UDG[i*nc*npv+18*npv+j];
            rEz = UDG[i*nc*npv+19*npv+j];

            u = ru/r;
            v = rv/r;
            w = rw/r;
            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;
            wz = (rwz - w*rz)/r;

            h = pn[i*npv*ncd + nd*npv + j];

            if (avgType == 0) {     // Sensor averaging
                x_DG[i*npv+j] = (k_h*h/((double) porder))*(ux+vy+wz);
            }
            else if (avgType == 1) { // Max surrogate averaging
                x = (k_h*h/((double) porder))*(ux+vy+wz);
                f_DG[i*npv+j] = (x-beta)*(atan(alpha*(x-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        DG_2_p1CG(x_p1CG, x_DG, mesh, master, app, ndims, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + nd*npv + j];
                f = (x_p1CG[i*npv+j]-beta)*(atan(alpha*(x_p1CG[i*npv+j]-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder))*(1.0+1.0/Minf)*f;
            }
        }
    }
    else if (avgType == 1) { // Max surrogate averaging
        DG_2_p1CG(f_p1CG, f_DG, mesh, master, app, ndims, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + nd*npv + j];
                avField[i*npv+j] = (k_h*h/((double) porder))*(1.0+1.0/Minf)*f_p1CG[i*npv+j];
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 1) {        // Max surrogate averaging
        delete[] f_p1CG; delete[] f_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_velDiv(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int nd = ndims[0];

    if (nd == 2)
        getAVfield_velDiv_2d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else if (nd == 3)
        getAVfield_velDiv_3d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else
        error("Invalid number of dimensions.\n");
}

void getAVfield_rhoSmoothness(double *avField_p1CG, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, tempstruct &temp, Int *ndims)
{
    // Note: This routine assumes the mesh size is the (nd+1)-th argument in pn
    // Note: This approach will introduce no AV for porder == 0, and almost no AV for porder == 1.
    
    Int inc = 1;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int porder = ndims[15];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int nqvR = master.nqvR;
    Int avgMethod = 1;
    
    double zero = 0.0, minusone = -1.0, one = 1.0;
    double pi = 3.141592653589793;
    double av_0;
    
    char chn = 'N';
    
    double *shapvtR = &master.shapvtR[0];
    double *gwvlR = &master.gwvlR[0];
    double *J_R = &temp.J_R[0];
    double *pR = &temp.pR[0];
    double *XxR = &temp.XxR[0];
    double *jacR = &temp.jacR[0];
    double *avField_DG = new double[npv*ne];
    double *rhoProj = new double[npv];
    double *orthComplRhoProj_g = new double[nqvR];
    double *rho;
    double *rho_g = new double[nqvR];
    
    double k_h = 1.5;
    double s_0 = 4.0 * log10 ( 1.0 / pow( (double) porder, 4.0) );     // TODO: Confirm this expression
    double kappa;
    if (porder == 0 || porder == 1)
        kappa = 1.0e-3;
    else
        kappa = abs(s_0) / 3.0; //1.0 + 3.0 * log10(porder);       // TODO: I came up with this value for kappa. Find out what Per is using.
    
    Int numFields2Proj = 1;
    Int *fields2Project = new Int[numFields2Proj];
    fields2Project[0] = 0;
    
    for (Int ie = 0; ie < ne; ie++) {
//         elementgeom(mesh, master, temp, ndims, ie);
        elementgeom(&pR[0], &J_R[0], &XxR[0], &jacR[0], &mesh.dgnodes[ie*npv*ncd], &shapvtR[0], &ndims[0], nqvR);
        
        projectUDGfields(&rhoProj[0], &master.projLowP[0], &UDG[ie*npv*nc], &fields2Project[0], numFields2Proj, npv);
        
        rho = &UDG[ie*npv*nc];
        DAXPY(&npv, &minusone, &rho[0], &inc, &rhoProj[0], &inc);
        double *orthComplRhoProj = &rhoProj[0];     // Contains the negative of the actual orthogonal complement
        
        // Compute rho and rhoProj at Gauss points:
        DGEMV(&chn, &nqvR, &npv, &one, &shapvtR[0], &nqvR, &orthComplRhoProj[0], &inc, &zero, &orthComplRhoProj_g[0], &inc);
        DGEMV(&chn, &nqvR, &npv, &one, &shapvtR[0], &nqvR, &rho[0], &inc, &zero, &rho_g[0], &inc);
        
        double numerator = 0.0, denominator = 0.0;
        for (Int ig = 0; ig < nqvR; ig ++) {
            numerator += jacR[ig] * gwvlR[ig] * orthComplRhoProj_g[ig] * orthComplRhoProj_g[ig];
            denominator += jacR[ig] * gwvlR[ig] * rho_g[ig] * rho_g[ig];
        }
        double s_e = log10 ( numerator / denominator );
        
        for (Int in = 0; in < npv; in++) {
            av_0 = k_h * pn[ie*ncd*npv + nd*npv + in] / ((double) porder);
            
            if (s_e < s_0 - kappa)
                avField_DG[ie*npv+in]  = 0.0;
            else if (s_e <= s_0 + kappa) {
                avField_DG[ie*npv+in]  = 0.5 * av_0 * (1.0 + sin( pi*(s_e - s_0) / (2*kappa) ));
            }
            else {
                avField_DG[ie*npv+in]  = av_0;
            }
        }
    }
    
    if (porder == 0) {
        for (Int i = 0; i < ne*npv; i++)
            avField_DG[i] = 0.0;
    }
    
    delete[] fields2Project;
    delete[] rhoProj; delete[] orthComplRhoProj_g; delete[] rho_g;
    
    DG_2_p1CG(avField_p1CG, avField_DG, mesh, master, app, &ndims[0], avgMethod);
    
    delete[] avField_DG;
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField_p1CG[0], &ndims[0]);
    #endif
}

void getAVfield_bulkViscosity_2d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];
    
    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    Int DucrosSensor = 1;       // Flag to determine whether the Ducros sensor is also applied to remove artificial viscosity from regions where vorticity dominates.
    double hMin = 1.0e-6;
    double hMax = 1.0e6;
    double alpha = 1.0e30;
    double beta = 0.01;
    double x_min = 0.0;
    double x_max = 2.0 * (2.0 / sqrt(gam*gam - 1.0));       // Twice the theoretical maximum in a shock wave for large M_{\infty} when using c^* as velocity scale.
    double k_h = 1.5;
    double eps_Ducros = 1.0e-20;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 1;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rE, rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u, v, ux, uy, vx, vy;
    double q, p, E, H, div_v, vort, p_reg, H_reg, DucrosRatio, c_star, h, f, x;
    double nx, ny, nNorm, h_tmp, hRef = 1.0;
    double *Minv;
    
    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    
    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *avField_DG;
    if (avgType == 0) {     // No smoothing
    }
    else if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3)          // AV smoothing
        avField_DG = new double[ne*npv];
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
            
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            rx = UDG[i*nc*npv+(4+nte)*npv+j];
            rux = UDG[i*nc*npv+(5+nte)*npv+j];
            rvx = UDG[i*nc*npv+(6+nte)*npv+j];
            rEx = UDG[i*nc*npv+(7+nte)*npv+j];
            ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(9+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(10+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(11+2*nte)*npv+j];

            u = ru/r;
            v = rv/r;
            E = rE/r;
            p = gam1*(rE - 0.5*r*(u*u+v*v));
            H = gam*E - 0.5*gam1*(u*u+v*v);
            
            p_reg = max( p_min , p*(atan(alpha*p)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            
            div_v = - (ux + vy);
            vort = - (vx - uy);
            DucrosRatio = div_v*div_v / (div_v*div_v + vort*vort + eps_Ducros);
            
            nx = rx;
            ny = ry;
            nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
            nx = nx / nNorm;
            ny = ny / nNorm;
            
            h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
            h = max(hMin,min(h_tmp,hMax));
            
            // The two following expressions for c_star are equivalent:
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );      // c_star = sqrt( (2 * gam * p_reg / r + (gam-1)*(u*u + v*v)) / (gam+1) );
            
            if (avgType == 0) {      // No smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
                if (DucrosSensor == 1)
                    avField[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = - (k_h*h/ ((double) porder)) * div_v / c_star;
                if (DucrosSensor == 1)
                    x_DG[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (DucrosSensor == 1)
                    f_DG[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 3) {      // AV smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (DucrosSensor == 1)
                    f *= DucrosRatio;
                avField_DG[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }

    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                rx = UDG[i*nc*npv+(4+nte)*npv+j];
                ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                H = gam*E - gam1*q;
                
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                nx = rx;
                ny = ry;
                nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                
                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                rx = UDG[i*nc*npv+(4+nte)*npv+j];
                ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                nx = rx;
                ny = ry;
                nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                
                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3)   // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3) {        // AV smoothing
        delete[] avField_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_bulkViscosity_3d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];

    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    Int DucrosSensor = 1;       // Flag to determine whether the Ducros sensor is also applied to remove artificial viscosity from regions where vorticity dominates.
    double hMin = 1.0e-6;
    double hMax = 0.02;
    double alpha = 1.0e30;
    double beta = 0.01;
    double x_min = 0.0;
    double x_max = 2.0 * (2.0 / sqrt(gam*gam - 1.0));       // Twice the theoretical maximum in a shock wave for large M_{\infty} when using c^* as velocity scale.
    double k_h = 1.5;
    double eps_Ducros = 1.0e-20;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 1;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rw, rE, rx, rux, rvx, rwx, rEx, ry, ruy, rvy, rwy, rEy, rz, ruz, rvz, rwz, rEz;
    double u, v, w, q, ux, uy, uz, vx, vy, vz, wx, wy, wz;
    double p, E, H, div_v, vort_x, vort_y, vort_z, vort, p_reg, H_reg, DucrosRatio, c_star, h, f, x;
    double nx, ny, nz, nNorm, h_tmp, hRef = 1.0;
    double *Minv;
    
    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    
    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *avField_DG;
    if (avgType == 0) {     // No smoothing
    }
    if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3) {        // AV smoothing
        avField_DG = new double[ne*npv];
    }

    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
            
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rw = UDG[i*nc*npv+3*npv+j];
            rE = UDG[i*nc*npv+4*npv+j];
            rx = UDG[i*nc*npv+(5+nte)*npv+j];
            rux = UDG[i*nc*npv+(6+nte)*npv+j];
            rvx = UDG[i*nc*npv+(7+nte)*npv+j];
            rwx = UDG[i*nc*npv+(8+nte)*npv+j];
            rEx = UDG[i*nc*npv+(9+nte)*npv+j];
            ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(11+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(12+2*nte)*npv+j];
            rwy = UDG[i*nc*npv+(13+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(14+2*nte)*npv+j];
            rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
            ruz = UDG[i*nc*npv+(16+3*nte)*npv+j];
            rvz = UDG[i*nc*npv+(17+3*nte)*npv+j];
            rwz = UDG[i*nc*npv+(18+3*nte)*npv+j];
            rEz = UDG[i*nc*npv+(19+3*nte)*npv+j];

            u = ru/r;
            v = rv/r;
            w = rw/r;
            E = rE/r;
            p = gam1*(rE - 0.5*r*(u*u+v*v+w*w));
            H = gam*E - 0.5*gam1*(u*u+v*v+w*w);
            
            p_reg = max( p_min , p*(atan(alpha*p)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            wx = (rwx - w*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            wy = (rwy - w*ry)/r;
            uz = (ruz - u*rz)/r;
            vz = (rvz - v*rz)/r;
            wz = (rwz - w*rz)/r;
            
            div_v = - (ux + vy + wz);
            vort_x = - (wy - vz);
            vort_y = - (uz - wx);
            vort_z = - (vx - uy);
            vort = sqrt(vort_x*vort_x + vort_y*vort_y + vort_z*vort_z);
            DucrosRatio = div_v*div_v / (div_v*div_v + vort*vort + eps_Ducros);
            
            nx = rx;
            ny = ry;
            nz = rz;
            nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
            nx = nx / nNorm;
            ny = ny / nNorm;
            nz = nz / nNorm;
            
            h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
            h = max(hMin,min(h_tmp,hMax));
            
            // The two following expressions for c_star are equivalent:
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );      // c_star = sqrt( (2 * gam * p_reg / r + (gam-1)*(u*u + v*v + w*w)) / (gam+1) );
            
            if(avgType == 0) {      // No smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
                if (DucrosSensor == 1)
                    avField[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = - (k_h*h/ ((double) porder)) * div_v / c_star;
                if (DucrosSensor == 1)
                    x_DG[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (DucrosSensor == 1)
                    f_DG[i*npv+j] *= DucrosRatio;
            }
            else if(avgType == 3) {      // AV smoothing
                x = - (k_h*h/ ((double) porder)) * div_v / c_star;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (DucrosSensor == 1)
                    f *= DucrosRatio;
                avField_DG[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }

    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                rx = UDG[i*nc*npv+(5+nte)*npv+j];
                ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
                rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                nx = rx;
                ny = ry;
                nz = rz;
                nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                nz = nz / nNorm;
                
                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                rx = UDG[i*nc*npv+(5+nte)*npv+j];
                ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
                rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                nx = rx;
                ny = ry;
                nz = rz;
                nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                nz = nz / nNorm;
                
                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3) { // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    }
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3) {        // AV smoothing
        delete[] avField_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_bulkViscosity(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int nd = ndims[0];

    if (nd == 2)
        getAVfield_bulkViscosity_2d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else if (nd == 3)
        getAVfield_bulkViscosity_3d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else
        error("Invalid number of dimensions.\n");
}

void getAVfield_thermalConductivity_2d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];

    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    double h_tmp, hRef = 1.0;
    double hMin = 1.0e-6;
    double hMax = 1.0e6;
    double alpha = 1.0e30;
    double beta = 1.0;
    double x_min = 0.0;
    double x_max = 2.0;
    double k_h = 1.5;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 0;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rE, rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u, v, ux, uy, vx, vy;
    double qx, px, Tx, qy, py, Ty, nx, ny, nNorm;
    double q, p, T, E, H, p_reg, H_reg, c_star, T_ref, norm_gradT_master, h, f, x;
    
    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    double T_0_infty = 1/(2*gam) + 1/(gam*gam1*Minf*Minf);      // Freestream total temperature. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    double T_infty = 1/(gam*gam1*Minf*Minf);                    // Freestream temperature. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    
    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *avField_DG, *M, *Minv;
    if (avgType == 0) {     // No smoothing
    }
    if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3) {        // AV smoothing
        avField_DG = new double[ne*npv];
    }
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            M = &mesh.M[i*npv*nd*nd+j*nd*nd];
            Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
            
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            rx = UDG[i*nc*npv+(4+nte)*npv+j];
            rux = UDG[i*nc*npv+(5+nte)*npv+j];
            rvx = UDG[i*nc*npv+(6+nte)*npv+j];
            rEx = UDG[i*nc*npv+(7+nte)*npv+j];
            ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(9+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(10+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(11+2*nte)*npv+j];

            u = ru/r;
            v = rv/r;
            q = 0.5*(u*u+v*v);
            E = rE/r;
            p = gam1*(rE - r*q);
            T = p /(gam1*r);
            H = gam*E - gam1*q;
            
            p_reg = max( p_min , p*(atan(alpha*p)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            
            qx  = u*ux + v*vx;
            px  = gam1*(rEx - rx*q - r*qx);
            Tx  = (px*r - p*rx)/(gam1*r*r);

            qy  = u*uy + v*vy;
            py  = gam1*(rEy - ry*q - r*qy);
            Ty  = (py*r - p*ry)/(gam1*r*r);
            
            nx = Tx;
            ny = Ty;
            nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
            nx = nx / nNorm;
            ny = ny / nNorm;
            
            h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
            h = max(hMin,min(h_tmp,hMax));
            
            norm_gradT_master = sqrt(M[0]*Tx*Tx + M[1]*Tx*Ty + M[2]*Ty*Tx + M[3]*Ty*Ty);
            
            // The two following expressions for c_star are equivalent:
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );      // c_star = sqrt( (2 * gam * p_reg / r + (gam-1)*(u*u + v*v)) / (gam+1) );
            
            // We use the following approximation to improve the smoothness of the AV field:
            T_ref = T_0_infty;
            
            if(avgType == 0) {      // No smoothing
                x = (k_h * hRef/ ((double) porder)) * norm_gradT_master / T_ref;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f = x_max;
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = (k_h * hRef/ ((double) porder)) * norm_gradT_master / T_ref;
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = (k_h * hRef / ((double) porder)) * norm_gradT_master / T_ref;
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f_DG[i*npv+j] = x_max;
            }
            else if(avgType == 3) {      // AV smoothing
                x = (k_h * hRef / ((double) porder)) * norm_gradT_master / T_ref;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0) {
                    f = x_max;
                }
                avField_DG[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }
    
    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                M = &mesh.M[i*npv*nd*nd+j*nd*nd];
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                rx = UDG[i*nc*npv+(4+nte)*npv+j];
                rux = UDG[i*nc*npv+(5+nte)*npv+j];
                rvx = UDG[i*nc*npv+(6+nte)*npv+j];
                rEx = UDG[i*nc*npv+(7+nte)*npv+j];
                ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
                ruy = UDG[i*nc*npv+(9+2*nte)*npv+j];
                rvy = UDG[i*nc*npv+(10+2*nte)*npv+j];
                rEy = UDG[i*nc*npv+(11+2*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                p  = gam1*(rE - r*q);
                T = p /(gam1*r);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                ux = (rux - u*rx)/r;
                vx = (rvx - v*rx)/r;
                uy = (ruy - u*ry)/r;
                vy = (rvy - v*ry)/r;
            
                qx  = u*ux + v*vx;
                px  = gam1*(rEx - rx*q - r*qx);
                Tx  = (px*r - p*rx)/(gam1*r*r);

                qy  = u*uy + v*vy;
                py  = gam1*(rEy - ry*q - r*qy);
                Ty  = (py*r - p*ry)/(gam1*r*r);

                nx = Tx;
                ny = Ty;
                nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                
                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0.0)
                    f = x_max;
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                M = &mesh.M[i*npv*nd*nd+j*nd*nd];
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                rx = UDG[i*nc*npv+(4+nte)*npv+j];
                rux = UDG[i*nc*npv+(5+nte)*npv+j];
                rvx = UDG[i*nc*npv+(6+nte)*npv+j];
                rEx = UDG[i*nc*npv+(7+nte)*npv+j];
                ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
                ruy = UDG[i*nc*npv+(9+2*nte)*npv+j];
                rvy = UDG[i*nc*npv+(10+2*nte)*npv+j];
                rEy = UDG[i*nc*npv+(11+2*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                p  = gam1*(rE - r*q);
                T = p /(gam1*r);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                ux = (rux - u*rx)/r;
                vx = (rvx - v*rx)/r;
                uy = (ruy - u*ry)/r;
                vy = (rvy - v*ry)/r;
            
                qx  = u*ux + v*vx;
                px  = gam1*(rEx - rx*q - r*qx);
                Tx  = (px*r - p*rx)/(gam1*r*r);

                qy  = u*uy + v*vy;
                py  = gam1*(rEy - ry*q - r*qy);
                Ty  = (py*r - p*ry)/(gam1*r*r);

                nx = Tx;
                ny = Ty;
                nNorm = sqrt(nx*nx + ny*ny + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;

                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*ny*nx + Minv[3]*ny*ny + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3) { // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    }
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3) {        // AV smoothing
        delete[] avField_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_thermalConductivity_3d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];
    
    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    double h_tmp, hRef = 1.0;
    double hMin = 1.0e-6;
    double hMax = 0.02;
    double alpha = 1.0e30;
    double beta = 1.0;
    double x_min = 0.0;
    double x_max = 2.0;
    double k_h = 1.5;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 0;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rw, rE, rx, rux, rvx, rwx, rEx, ry, ruy, rvy, rwy, rEy, rz, ruz, rvz, rwz, rEz;
    double u, v, w, q, p, T, E, H, p_reg, H_reg, c_star, T_ref, norm_gradT_master, h, f, x;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, qx, px, Tx, qy, py, Ty, qz, pz, Tz;
    double nx, ny, nz, nNorm;
    
    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    double T_0_infty = 1/(2*gam) + 1/(gam*gam1*Minf*Minf);      // Freestream total temperature. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    double T_infty = 1/(gam*gam1*Minf*Minf);                    // Freestream temperature. We assume the reference velocity for non-dimensionalization is the freestream velocity.

    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *M, *Minv, *avField_DG;
    if (avgType == 0) {     // No smoothing
    }
    if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3) {        // AV smoothing
        avField_DG = new double[ne*npv];
    }

    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            M = &mesh.M[i*npv*nd*nd+j*nd*nd];
            Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
            
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rw = UDG[i*nc*npv+3*npv+j];
            rE = UDG[i*nc*npv+4*npv+j];
            rx = UDG[i*nc*npv+(5+nte)*npv+j];
            rux = UDG[i*nc*npv+(6+nte)*npv+j];
            rvx = UDG[i*nc*npv+(7+nte)*npv+j];
            rwx = UDG[i*nc*npv+(8+nte)*npv+j];
            rEx = UDG[i*nc*npv+(9+nte)*npv+j];
            ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(11+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(12+2*nte)*npv+j];
            rwy = UDG[i*nc*npv+(13+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(14+2*nte)*npv+j];
            rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
            ruz = UDG[i*nc*npv+(16+3*nte)*npv+j];
            rvz = UDG[i*nc*npv+(17+3*nte)*npv+j];
            rwz = UDG[i*nc*npv+(18+3*nte)*npv+j];
            rEz = UDG[i*nc*npv+(19+3*nte)*npv+j];
            
            u = ru/r;
            v = rv/r;
            w = rw/r;
            E = rE/r;
            q = 0.5*(u*u+v*v+w*w);
            p = gam1*(rE - r*q);
            T = p /(gam1*r);
            H = gam*E - gam1*q;
            
            p_reg = max( p_min , p*(atan(alpha*p)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            wx = (rwx - w*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            wy = (rwy - w*ry)/r;
            uz = (ruz - u*rz)/r;
            vz = (rvz - v*rz)/r;
            wz = (rwz - w*rz)/r;

            qx  = u*ux + v*vx + w*wx;
            px  = gam1*(rEx - rx*q - r*qx);
            Tx  = (px*r - p*rx)/(gam1*r*r);

            qy  = u*uy + v*vy + w*wy;
            py  = gam1*(rEy - ry*q - r*qy);
            Ty  = (py*r - p*ry)/(gam1*r*r);

            qz  = u*uz + v*vz + w*wz;
            pz  = gam1*(rEz - rz*q - r*qz);
            Tz  = (pz*r - p*rz)/(gam1*r*r);
            
            nx = Tx;
            ny = Ty;
            nz = Tz;
            nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
            nx = nx / nNorm;
            ny = ny / nNorm;
            nz = nz / nNorm;
            
            h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
            h = max(hMin,min(h_tmp,hMax));
            
            norm_gradT_master = sqrt(M[0]*Tx*Tx + M[1]*Tx*Ty + M[2]*Tx*Tz + M[3]*Ty*Tx + M[4]*Ty*Ty + M[5]*Ty*Tz + M[6]*Tz*Tx + M[7]*Tz*Ty + M[8]*Tz*Tz);
            
            // The two following expressions for c_star are equivalent:
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );      // c_star = sqrt( (2 * gam * p_reg / r + (gam-1)*(u*u + v*v + w*w)) / (gam+1) );
            
            T_ref = T_0_infty;
            
            if(avgType == 0) {      // No smoothing
                x = (k_h*hRef/ ((double) porder)) * norm_gradT_master / T_ref;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f = x_max;
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = (k_h*hRef/ ((double) porder)) * norm_gradT_master / T_ref;
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = (k_h*hRef/ ((double) porder)) * norm_gradT_master / T_ref;
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f_DG[i*npv+j] = x_max;
            }
            else if(avgType == 3) {      // AV smoothing
                x = (k_h*hRef/ ((double) porder)) * norm_gradT_master / T_ref;
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f = x_max;
                avField_DG[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }

    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                M = &mesh.M[i*npv*nd*nd+j*nd*nd];
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                rx = UDG[i*nc*npv+(5+nte)*npv+j];
                rux = UDG[i*nc*npv+(6+nte)*npv+j];
                rvx = UDG[i*nc*npv+(7+nte)*npv+j];
                rwx = UDG[i*nc*npv+(8+nte)*npv+j];
                rEx = UDG[i*nc*npv+(9+nte)*npv+j];
                ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
                ruy = UDG[i*nc*npv+(11+2*nte)*npv+j];
                rvy = UDG[i*nc*npv+(12+2*nte)*npv+j];
                rwy = UDG[i*nc*npv+(13+2*nte)*npv+j];
                rEy = UDG[i*nc*npv+(14+2*nte)*npv+j];
                rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
                ruz = UDG[i*nc*npv+(16+3*nte)*npv+j];
                rvz = UDG[i*nc*npv+(17+3*nte)*npv+j];
                rwz = UDG[i*nc*npv+(18+3*nte)*npv+j];
                rEz = UDG[i*nc*npv+(19+3*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                p  = gam1*(rE - r*q);
                T = p /(gam1*r);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                ux = (rux - u*rx)/r;
                vx = (rvx - v*rx)/r;
                wx = (rwx - w*rx)/r;
                uy = (ruy - u*ry)/r;
                vy = (rvy - v*ry)/r;
                wy = (rwy - w*ry)/r;
                uz = (ruz - u*rz)/r;
                vz = (rvz - v*rz)/r;
                wz = (rwz - w*rz)/r;
                
                qx  = u*ux + v*vx + w*wx;
                px  = gam1*(rEx - rx*q - r*qx);
                Tx  = (px*r - p*rx)/(gam1*r*r);

                qy  = u*uy + v*vy + w*wy;
                py  = gam1*(rEy - ry*q - r*qy);
                Ty  = (py*r - p*ry)/(gam1*r*r);

                qz  = u*uz + v*vz + w*wz;
                pz  = gam1*(rEz - rz*q - r*qz);
                Tz  = (pz*r - p*rz)/(gam1*r*r);

                nx = Tx;
                ny = Ty;
                nz = Tz;
                nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                nz = nz / nNorm;

                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                if (T <= 0)
                    f = x_max;
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                M = &mesh.M[i*npv*nd*nd+j*nd*nd];
                Minv = &mesh.Minv[i*npv*nd*nd+j*nd*nd];
                
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                rx = UDG[i*nc*npv+(5+nte)*npv+j];
                rux = UDG[i*nc*npv+(6+nte)*npv+j];
                rvx = UDG[i*nc*npv+(7+nte)*npv+j];
                rwx = UDG[i*nc*npv+(8+nte)*npv+j];
                rEx = UDG[i*nc*npv+(9+nte)*npv+j];
                ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
                ruy = UDG[i*nc*npv+(11+2*nte)*npv+j];
                rvy = UDG[i*nc*npv+(12+2*nte)*npv+j];
                rwy = UDG[i*nc*npv+(13+2*nte)*npv+j];
                rEy = UDG[i*nc*npv+(14+2*nte)*npv+j];
                rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
                ruz = UDG[i*nc*npv+(16+3*nte)*npv+j];
                rvz = UDG[i*nc*npv+(17+3*nte)*npv+j];
                rwz = UDG[i*nc*npv+(18+3*nte)*npv+j];
                rEz = UDG[i*nc*npv+(19+3*nte)*npv+j];
                
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                p  = gam1*(rE - r*q);
                T = p /(gam1*r);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                ux = (rux - u*rx)/r;
                vx = (rvx - v*rx)/r;
                wx = (rwx - w*rx)/r;
                uy = (ruy - u*ry)/r;
                vy = (rvy - v*ry)/r;
                wy = (rwy - w*ry)/r;
                uz = (ruz - u*rz)/r;
                vz = (rvz - v*rz)/r;
                wz = (rwz - w*rz)/r;
                
                qx  = u*ux + v*vx + w*wx;
                px  = gam1*(rEx - rx*q - r*qx);
                Tx  = (px*r - p*rx)/(gam1*r*r);

                qy  = u*uy + v*vy + w*wy;
                py  = gam1*(rEy - ry*q - r*qy);
                Ty  = (py*r - p*ry)/(gam1*r*r);

                qz  = u*uz + v*vz + w*wz;
                pz  = gam1*(rEz - rz*q - r*qz);
                Tz  = (pz*r - p*rz)/(gam1*r*r);

                nx = Tx;
                ny = Ty;
                nz = Tz;
                nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-20);
                nx = nx / nNorm;
                ny = ny / nNorm;
                nz = nz / nNorm;

                h_tmp = hRef / sqrt(Minv[0]*nx*nx + Minv[1]*nx*ny + Minv[2]*nx*nz + Minv[3]*ny*nx + Minv[4]*ny*ny + Minv[5]*ny*nz + Minv[6]*nz*nx + Minv[7]*nz*ny + Minv[8]*nz*nz + 1.0e-20);
                h = max(hMin,min(h_tmp,hMax));
                
                avField[i*npv+j] = gam * (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3) { // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    }
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3) {        // AV smoothing
        delete[] avField_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_thermalConductivity(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int nd = ndims[0];

    if (nd == 2)
        getAVfield_thermalConductivity_2d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else if (nd == 3)
        getAVfield_thermalConductivity_3d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else
        error("Invalid number of dimensions.\n");
}

void getAVfield_molecularViscosity_2d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];

    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    double hMin = 1.0e-6;
    double hMax = 1.0e6;
    double alpha = 1.0e30;
    double beta = 1.0;
    double x_min = 0.0;
    double x_max = 2.0;
    double k_h = 1.5;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 0;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rE, rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u, v, ux, uy, vx, vy;
    double E, q, p, H, p_reg, H_reg;
    double S_xx, S_xy, S_yx, S_yy;
    double Q_S, h, f, x, c_star;
    
    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *avField_DG;
    if (avgType == 0) {     // No smoothing
    }
    if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3) {        // AV smoothing
        avField_DG = new double[ne*npv];
    }

    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            rx = UDG[i*nc*npv+(4+nte)*npv+j];
            rux = UDG[i*nc*npv+(5+nte)*npv+j];
            rvx = UDG[i*nc*npv+(6+nte)*npv+j];
            rEx = UDG[i*nc*npv+(7+nte)*npv+j];
            ry = UDG[i*nc*npv+(8+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(9+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(10+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(11+2*nte)*npv+j];

            u = ru/r;
            v = rv/r;
            q = 0.5*(u*u+v*v);
            E = rE/r;
            p = gam1*(rE - r*q);
            H = gam*E - gam1*q;
            
            p_reg = max( p_min , p*(atan(alpha*p)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            
            S_xx = 0.5 * (ux + ux);
            S_xy = 0.5 * (uy + vx);
            S_yx = 0.5 * (vx + uy);
            S_yy = 0.5 * (vy + vy);
            
            Q_S = max( 0.0 , (S_xx*S_xx + S_xy*S_xy + S_yx*S_yx + S_yy*S_yy) - (ux+vy)*(ux+vy));
            
            h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
            
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
            
            if(avgType == 0) {      // No smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
            else if(avgType == 3) {      // AV smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField_DG[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }

    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rE = UDG[i*nc*npv+3*npv+j];
                u = ru/r;
                v = rv/r;
                E = rE/r;
                q = 0.5*(u*u+v*v);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
                
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3) { // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    }
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3) {        // AV smoothing
        delete[] avField_DG;
    }
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_molecularViscosity_3d(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int i, j;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc = ndims[19];
    Int porder = app.porder;
    
    double *param = &app.param[0];

    double gam = param[0];
    double Minf = param[4];
    double gam1 = gam - 1.0;
    double pi = 3.141592653589793;
    
    double hMin = 1.0e-6;
    double hMax = 0.02;
    double alpha = 1.0e30;
    double beta = 1.0;
    double x_min = 0.0;
    double x_max = 2.0;
    double k_h = 1.5;
    double p_min = 1.0e-6;
    double H_min = 1.0e-6;
    Int numSmoothing = 1;
    Int avgMethod = 0;
    
    Int nte;
    getNumTurbEqs(app, &nte);
    
    double r, ru, rv, rw, rE, rx, rux, rvx, rwx, rEx, ry, ruy, rvy, rwy, rEy, rz, ruz, rvz, rwz, rEz;
    double u, v, w, Q_S, h, f, x, c_star;
    double E, q, p, H, p_reg, H_reg;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz;
    double S_xx, S_xy, S_xz, S_yx, S_yy, S_yz, S_zx, S_zy, S_zz;
    
    double *x_p1CG, *x_DG, *f_p1CG, *f_DG, *avField_DG;
    if (avgType == 0) {     // No smoothing
    }
    if (avgType == 1) {     // Sensor smoothing
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }
    else if (avgType == 3) {        // AV smoothing
        avField_DG = new double[ne*npv];
    }
    
    double c_star_infty = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );    // Freestream critical speed of sound. We assume the reference velocity for non-dimensionalization is the freestream velocity.
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rw = UDG[i*nc*npv+3*npv+j];
            rE = UDG[i*nc*npv+4*npv+j];
            rx = UDG[i*nc*npv+(5+nte)*npv+j];
            rux = UDG[i*nc*npv+(6+nte)*npv+j];
            rvx = UDG[i*nc*npv+(7+nte)*npv+j];
            rwx = UDG[i*nc*npv+(8+nte)*npv+j];
            rEx = UDG[i*nc*npv+(9+nte)*npv+j];
            ry = UDG[i*nc*npv+(10+2*nte)*npv+j];
            ruy = UDG[i*nc*npv+(11+2*nte)*npv+j];
            rvy = UDG[i*nc*npv+(12+2*nte)*npv+j];
            rwy = UDG[i*nc*npv+(13+2*nte)*npv+j];
            rEy = UDG[i*nc*npv+(14+2*nte)*npv+j];
            rz = UDG[i*nc*npv+(15+3*nte)*npv+j];
            ruz = UDG[i*nc*npv+(16+3*nte)*npv+j];
            rvz = UDG[i*nc*npv+(17+3*nte)*npv+j];
            rwz = UDG[i*nc*npv+(18+3*nte)*npv+j];
            rEz = UDG[i*nc*npv+(19+3*nte)*npv+j];

            u = ru/r;
            v = rv/r;
            w = rw/r;
            E = rE/r;
            q = 0.5*(u*u+v*v+w*w);
            H = gam*E - gam1*q;
            H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
            
            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            wx = (rwx - w*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            wy = (rwy - w*ry)/r;
            uz = (ruz - u*rz)/r;
            vz = (rvz - v*rz)/r;
            wz = (rwz - w*rz)/r;
            
            S_xx = 0.5 * (ux + ux);
            S_xy = 0.5 * (uy + vx);
            S_xz = 0.5 * (uz + wx);
            S_yx = 0.5 * (vx + uy);
            S_yy = 0.5 * (vy + vy);
            S_yz = 0.5 * (vz + wy);
            S_zx = 0.5 * (wx + uz);
            S_zy = 0.5 * (wy + vz);
            S_zz = 0.5 * (wz + wz);
            
            Q_S = max( 0.0 , (S_xx*S_xx + S_xy*S_xy + S_xz*S_xz + S_yx*S_yx + S_yy*S_yy + S_yz*S_yz + S_zx*S_zx + S_zy*S_zy + S_zz*S_zz) - (ux+vy+wz)*(ux+vy+wz));
            
            h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
            c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
            
            if(avgType == 0) {      // No smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
            else if(avgType == 1) {      // Sensor smoothing
                x_DG[i*npv+j] = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
            }
            else if(avgType == 2) {      // Max surrogate smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f_DG[i*npv+j] = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f_DG[i*npv+j] = f_DG[i*npv+j] - (f_DG[i*npv+j]-x_max)*(atan(alpha*(f_DG[i*npv+j]-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
            else if(avgType == 3) {      // AV smoothing
                x = (k_h*h/ ((double) porder)) * sqrt(Q_S) / sqrt(1 + c_star_infty*c_star_infty);
                f = x_min + (x-beta-x_min)*(atan(alpha*(x-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField_DG[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }

    if (avgType == 1) {     // Sensor smoothing
        applySmoothing(x_p1CG, x_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
                
                f = x_min + (x_p1CG[i*npv+j]-beta-x_min)*(atan(alpha*(x_p1CG[i*npv+j]-beta-x_min))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                f = f - (f-x_max)*(atan(alpha*(f-x_max))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f;
            }
        }
    }
    else if (avgType == 2) { // Max surrogate smoothing
        applySmoothing(f_p1CG, f_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                r = UDG[i*nc*npv+0*npv+j];
                ru = UDG[i*nc*npv+1*npv+j];
                rv = UDG[i*nc*npv+2*npv+j];
                rw = UDG[i*nc*npv+3*npv+j];
                rE = UDG[i*nc*npv+4*npv+j];
                u = ru/r;
                v = rv/r;
                w = rw/r;
                E = rE/r;
                q = 0.5*(u*u+v*v+w*w);
                H = gam*E - gam1*q;
                H_reg = max( H_min , H*(atan(alpha*H)/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0);
                
                c_star = sqrt( (2*gam1*H_reg) / (gam+1) );
                
                h = max(hMin,min(pn[i*npv*ncd + nd*npv + j],hMax));
                
                avField[i*npv+j] = (k_h*h/((double) porder)) * sqrt(u*u + v*v + w*w + c_star*c_star) * f_p1CG[i*npv+j];
            }
        }
    }
    else if (avgType == 3) { // AV smoothing
        applySmoothing(avField, avField_DG, mesh, master, app, ndims, numSmoothing, avgMethod);
    }
    
    if (avgType == 1) {             // Sensor smoothing
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 2) {        // Max surrogate smoothing
        delete[] f_p1CG; delete[] f_DG;
    }
    else if (avgType == 3)        // AV smoothing
        delete[] avField_DG;
    
    #ifdef  HAVE_MPI
    // Communicate AV field among MPI workers
    sendrecvScalarField(sys, &avField[0], &ndims[0]);
    #endif
}

void getAVfield_molecularViscosity(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, Int avgType)
{
    Int nd = ndims[0];

    if (nd == 2)
        getAVfield_molecularViscosity_2d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else if (nd == 3)
        getAVfield_molecularViscosity_3d(avField, UDG, pn, sys, mesh, master, app, ndims, avgType);
    else
        error("Invalid number of dimensions.\n");
}

void getAVfield(double *avField, double *UDG, double *pn, sysstruct &sys, meshstruct &mesh, masterstruct &master, appstruct &app, tempstruct &temp, Int *ndims, Int avgType, Int AVflag)
{
    Int ne = ndims[5];
    Int npv = ndims[9];
    
    if (AVflag == 6)
        getAVfield_velDiv(avField, UDG, pn, sys, mesh, master, app, &ndims[0], avgType);
    else if (AVflag == 8)
        getAVfield_rhoSmoothness(avField, UDG, pn, sys, mesh, master, app, temp, &ndims[0]);
    else if (AVflag == 9)
        getAVfield_bulkViscosity(avField, UDG, pn, sys, mesh, master, app, &ndims[0], avgType);
    else if (AVflag == 10) {
        avgType = 3;
        getAVfield_bulkViscosity(&avField[0*npv*ne], UDG, pn, sys, mesh, master, app, &ndims[0], avgType);
        checkAVfield(&avField[0*npv*ne], &ndims[0]);
        avgType = 3;
        getAVfield_thermalConductivity(&avField[1*npv*ne], UDG, pn, sys, mesh, master, app, &ndims[0], avgType);
        checkAVfield(&avField[1*npv*ne], &ndims[0]);
        avgType = 3;
        getAVfield_molecularViscosity(&avField[2*npv*ne], UDG, pn, sys, mesh, master, app, &ndims[0], avgType);
        checkAVfield(&avField[2*npv*ne], &ndims[0]);
    }
    else
        error("AVflag has invalid value in getAVfield.\n");
}

#endif
