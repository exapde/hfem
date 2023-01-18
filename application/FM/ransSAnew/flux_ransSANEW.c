
#include "flux_ransSA2dNEW.c"
#include "flux_ransSA3dNEW.c"
#include "AVflux_ransSA.c"

// Written by: C. Nguyen & P. Fernandez

void flux_ransSANEW(double *f, double *f_udg, double *pg, double *udg, double *odg, meshstruct &mesh, masterstruct &master, appstruct &app, double *param,
                 double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            flux_ransSA2dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 3:
            flux_ransSA3dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        default:
            exit(-1);
            break;
    }
    
    // Add artificial viscosity (if necessary)
    if (app.AVflag != 0)
        AVflux_ransSA(f, f_udg, pg, udg, odg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
}
