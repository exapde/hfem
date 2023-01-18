
#include "fhat_leuq2d.c"
#include "fhat_leuq3d.c"

void fhat_leuq(double *f, double *f_udg, double * f_uh, double *pg, double *udg, double *odg, 
        double *uh, double *nl, meshstruct &mesh, masterstruct &master, appstruct &app, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{    
    switch (nd) {
        case 2:            
            fhat_leuq2d(f, f_udg, f_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);            
            break;
        case 3:
            fhat_leuq3d(f, f_udg, f_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }    
}


