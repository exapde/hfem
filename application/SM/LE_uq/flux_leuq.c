
#include "flux_leuq2d.c"
#include "flux_leuq3d.c"

void flux_leuq(double *f, double *f_udg, double *pg, double *udg, double *odg, meshstruct &mesh,
        masterstruct &master, appstruct &app, double *param, double time, int ng, 
        int nc, int ncu, int nd, int ncd, int computeJacobian)
{     
    switch (nd) {
        case 2:            
            flux_leuq2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);            
            break;
        case 3:
            flux_leuq3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }    
}


