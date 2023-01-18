
#include "flux_ns2dNEW.c"
#include "flux_ns3dNEW.c"
#include "AVflux_ns.c"
#include "SGSflux_ns.c"

// Written by: C. Nguyen & P. Fernandez

void flux_nsNEW(double *f, double *f_udg, double *pg, double *udg, double *odg, meshstruct &mesh, masterstruct &master, appstruct &app, double *param,
             double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            flux_ns2dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 3:
            flux_ns3dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        default:
            exit(-1);
            break;
    }

    // Add artificial viscosity (if necessary)
    if (app.AVflag != 0)
        AVflux_ns(f, f_udg, pg, udg, odg, mesh, master, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);

// void AVflux_ns(double *f, double *f_udg, double *pg, double *udg, double *avField, meshstruct &mesh, masterstruct &master, appstruct &app, double *param,
//                double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
    
    // Add SGS effects (if necessary)
    if (app.SGSmodel != 0)
        SGSflux_ns(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
    
}
