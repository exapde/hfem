
#include "flux_smagSGS_ns3d.c"
#include "flux_waleSGS_ns3d.c"
#include "flux_VremanSGS_ns3d.c"

// Written by: C. Nguyen & P. Fernandez

void SGSflux_ns(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
               double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            error("No SGS models implemented for 2D flows.");
            break;
        case 3:
            switch (app.SGSmodel) {
                case 1:
                    flux_smagSGS_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 2:
                    flux_waleSGS_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 3:
                    flux_VremanSGS_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                default:
                    error("SGS model not implemented.\n");
            }
            break;
        default:
            error("Number of dimensions not available.");
            break;
    }
}
