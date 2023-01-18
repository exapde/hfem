#ifndef __COMPUTEDYNAMICSMAGCONSTANT
#define __COMPUTEDYNAMICSMAGCONSTANT

// Written by: C. Nguyen & P. Fernandez

void computeDynamicSmagConstant(double *C_d, sysstruct &sys, meshstruct &mesh, 
        masterstruct &master, appstruct &app, double *UDG, double *Delta, 
        Int projPorder, Int *ndims, Int minimizationTarget)
{
    // C_d: npv / ne
    // Delta: ne
    
    Int ie, j, k, l, numFields2Proj, inc = 1;
    Int nd = ndims[0];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int porder = ndims[15];
    Int nc = ndims[19];
    Int nd2 = nd*nd;
    Int avgMethod = 1;
    
    if (nd != 3)
        error("Smagorinsky SGS model only available for 3D solver.\n");
    
    double alpha = ((double) projPorder) / ((double) porder);
    
    double r, ru, rv, rw, rx, rux, rvx, rwx, ry, ruy, rvy, rwy, rz, ruz, rvz, rwz;
    double u, v, w, ux, vx, wx, uy, vy, wy, uz, vz, wz;
    double uu, uv, uw, vu, vv, vw, wu, wv, ww;
    double Delta_node;
    
    double *projMatrix = new double[npv*npv];
    
    double *ui = new double[nd*npv];
    double *uiuj = new double[nd2*npv];
    double *Sij = new double[nd2*npv];
    double *S_mag = new double[npv];
    double *S_magSij = new double[nd2*npv];
    
    double *uiHat = new double[nd*npv];
    double *uiujHat = new double[nd2*npv];
    double *SijHat = new double[nd2*npv];
    double *S_magHat = new double[npv];
    double *S_magSij_Hat = new double[nd2*npv];
    
    double *L = new double[nd2];
    double *M = new double[nd2];
    double *C_d_DG = new double[ne*npv];
    
    getProjectionMatrix(&ndims[0], projPorder, &projMatrix[0]);
    
    for (ie = 0; ie < ne; ie++) {
        for (j = 0; j < npv; j++) {
            r = UDG[ie*nc*npv+0*npv+j];
            ru = UDG[ie*nc*npv+1*npv+j];
            rv = UDG[ie*nc*npv+2*npv+j];
            rw = UDG[ie*nc*npv+3*npv+j];
            rx = UDG[ie*nc*npv+5*npv+j];
            rux = UDG[ie*nc*npv+6*npv+j];
            rvx = UDG[ie*nc*npv+7*npv+j];
            rwx = UDG[ie*nc*npv+8*npv+j];
            ry = UDG[ie*nc*npv+10*npv+j];
            ruy = UDG[ie*nc*npv+11*npv+j];
            rvy = UDG[ie*nc*npv+12*npv+j];
            rwy = UDG[ie*nc*npv+13*npv+j];
            rz = UDG[ie*nc*npv+15*npv+j];
            ruz = UDG[ie*nc*npv+16*npv+j];
            rvz = UDG[ie*nc*npv+17*npv+j];
            rwz = UDG[ie*nc*npv+18*npv+j];

            u = ru/r;
            v = rv/r;
            w = rw/r;

            uu = u*u;
            uv = u*v;
            uw = u*w;
            vu = v*u;
            vv = v*v;
            vw = v*w;
            wu = w*u;
            wv = w*v;
            ww = w*w;

            ux = (rux - u*rx)/r;
            vx = (rvx - v*rx)/r;
            wx = (rwx - w*rx)/r;
            uy = (ruy - u*ry)/r;
            vy = (rvy - v*ry)/r;
            wy = (rwy - w*ry)/r;
            uz = (ruz - u*rz)/r;
            vz = (rvz - v*rz)/r;
            wz = (rwz - w*rz)/r;

            ui[0*npv+j] = u;
            ui[1*npv+j] = v;
            ui[2*npv+j] = w;

            uiuj[0*npv+j] = u*u;
            uiuj[1*npv+j] = u*v;
            uiuj[2*npv+j] = u*w;
            uiuj[3*npv+j] = v*u;
            uiuj[4*npv+j] = v*v;
            uiuj[5*npv+j] = v*w;
            uiuj[6*npv+j] = w*u;
            uiuj[7*npv+j] = w*v;
            uiuj[8*npv+j] = w*w;

            Sij[0*npv+j] = ux;
            Sij[1*npv+j] = 0.5 * (vx + uy);
            Sij[2*npv+j] = 0.5 * (wx + uz);
            Sij[3*npv+j] = 0.5 * (uy + vx);
            Sij[4*npv+j] = vy;
            Sij[5*npv+j] = 0.5 * (wy + vz);
            Sij[6*npv+j] = 0.5 * (uz + wx);
            Sij[7*npv+j] = 0.5 * (vz + wy);
            Sij[8*npv+j] = wz;
            
            S_mag[j] = 0.0;
            for (k = 0; k < nd2; k++)
                S_mag[j] += Sij[k*npv+j]*Sij[k*npv+j];
            S_mag[j] = sqrt(2.0*S_mag[j]);

            for (k = 0; k < nd2; k++)
                S_magSij[k*npv+j] = S_mag[j]*Sij[k*npv+j];
        }
        
        // Project onto lower p space:
        numFields2Proj = nd;
        projectFields(&uiHat[0], &ui[0], &projMatrix[0], numFields2Proj, npv);
        numFields2Proj = nd2;
        projectFields(&uiujHat[0], &uiuj[0], &projMatrix[0], numFields2Proj, npv);
        projectFields(&S_magSij_Hat[0], &S_magSij[0], &projMatrix[0], numFields2Proj, npv);
        projectFields(&SijHat[0], &Sij[0], &projMatrix[0], numFields2Proj, npv);
        numFields2Proj = 1;
        projectFields(&S_magHat[0], &S_mag[0], &projMatrix[0], numFields2Proj, npv);
        
        for (j = 0; j < npv; j++) {
            Delta_node = Delta[ie];
            
            for (k = 0; k < nd; k++)
                for (l = 0; l < nd; l++)
                    L[k*nd+l] = uiujHat[(k*nd+l)*npv+j] - uiHat[l*npv+j]*uiHat[k*npv+j];
            
            for (k = 0; k < nd; k++)
                for (l = 0; l < nd; l++) 
                    M[k*nd+l] = -2.0 * Delta_node*Delta_node * (alpha*alpha*S_magHat[j]*SijHat[(k*nd+l)*npv+j] - S_magSij_Hat[(k*nd+l)*npv+j]);
            
            if (minimizationTarget == 0)        // Least squares,  Lilly (1992), Ghosal et al. (1995)
                C_d_DG[ie*npv+j] = sqrt( DDOT(&nd2, &L[0], &inc, &M[0], &inc) / DDOT(&nd2, &M[0], &inc, &M[0], &inc) );
            else if (minimizationTarget == 1)   // Contract stresses with the strain-rate tensor, Germano et al. (1991)
                C_d_DG[ie*npv+j] = sqrt( DDOT(&nd2, &L[0], &inc, &SijHat[0], &inc) / DDOT(&nd2, &M[0], &inc, &SijHat[0], &inc) );
        }
    }
    
    // Smooth C_d field:
    DG_2_p1CG(C_d, &C_d_DG[0], mesh, master, app, ndims, avgMethod);
    
    #ifdef HAVE_MPI
    // Communicate C_d field among MPI workers:
    sendrecvScalarField(sys, &C_d[0], &ndims[0]);
    #endif
    
    delete[] projMatrix;
    
    delete[] ui; delete[] uiuj;
    delete[] Sij; delete[] S_mag;
    delete[] S_magSij;
    
    delete[] uiHat; delete[] uiujHat;
    delete[] SijHat; delete[] S_magHat;
    delete[] S_magSij_Hat;
    
    delete[] L; delete[] M;
    delete[] C_d_DG;
}

#endif
