#ifndef __CHECKFEASIBLESOLUTION_FM
#define __CHECKFEASIBLESOLUTION_FM

// Written by: C. Nguyen & P. Fernandez

Int checkFeasibleSolution_FM_2d(double *UDG, double *UH, appstruct &app, Int* ndims)
{
    Int i, j, feasibleSolution = 1;
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ndh = ndims[8];
    Int nch = ndims[23];
    
    double r, ru, rv, rE, u, v, q, p;
    
    double *param = &app.param[0];

    double gam = param[0];
    double gam1 = gam - 1.0;
    
    double r_min = 0;
    double p_min = -1;
    
    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            
            u = ru/r;
            v = rv/r;
            q = 0.5*(u*u+v*v);
            p = gam1*(rE - r*q);
    
            if (r < r_min || p < p_min) {
                feasibleSolution = 0;
                printf("UDG r = %g, p = %g\n", r, p);
                break;
            }
        }
        if (feasibleSolution == 0)
            break;
    }
    
    if (feasibleSolution == 1) {
        for (i = 0; i < ndh; i++) {
            r = UH[i*nch+0];
            ru = UH[i*nch+1];
            rv = UH[i*nch+2];
            rE = UH[i*nch+3];
            
            u = ru/r;
            v = rv/r;
            q = 0.5*(u*u+v*v);
            p = gam1*(rE - r*q);
    
            if (r < r_min || p < p_min) {
                feasibleSolution = 0;
                printf("UH r = %g, p = %g\n", r, p);
                break;
            }
        }
    }
    
    return feasibleSolution;
}

Int checkFeasibleSolution_FM_3d(double *UDG, double *UH, appstruct &app, Int* ndims)
{
    Int i, j, feasibleSolution = 1;
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ndh = ndims[8];
    Int nch = ndims[23];
    
    double r, ru, rv, rw, rE, u, v, w, q, p;
    
    double *param = &app.param[0];

    double gam = param[0];
    double gam1 = gam - 1.0;
    
    double r_min = 1.0e-4;
    double p_min = -0.5;
    
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
            q = 0.5*(u*u+v*v+w*w);
            p = gam1*(rE - r*q);
    
            if (r < r_min || p < p_min) {
                feasibleSolution = 0;
                break;
            }
        }
        if (feasibleSolution == 0)
            break;
    }
    
    if (feasibleSolution == 1) {
        for (i = 0; i < ndh; i++) {
            r = UH[i*nch+0];
            ru = UH[i*nch+1];
            rv = UH[i*nch+2];
            rw = UH[i*nch+3];
            rE = UH[i*nch+4];
            
            u = ru/r;
            v = rv/r;
            w = rw/r;
            q = 0.5*(u*u+v*v+w*w);
            p = gam1*(rE - r*q);
    
            if (r < r_min || p < p_min) {
                feasibleSolution = 0;
                break;
            }
        }
    }
    
    return feasibleSolution;
}

Int checkFeasibleSolution_FM(double *UDG, double *UH, appstruct &app, Int* ndims)
{
    Int nd = ndims[0];
    Int feasibleSolution;
    
    switch (nd) {
        case 2:
            feasibleSolution = checkFeasibleSolution_FM_2d(UDG, UH, app, ndims);
            break;
        case 3:
            feasibleSolution = checkFeasibleSolution_FM_3d(UDG, UH, app, ndims);
            break;
        default:
            exit(-1);
            break;
    }
    
    return feasibleSolution;
}

#endif
