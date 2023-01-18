
#include "source_ns2dNEW.c"
#include "source_ns3dNEW.c"
#include "source_inertial_ns2d.c"
#include "source_inertial_ns3d.c"

// Written by: C. Nguyen & P. Fernandez

void source_nsNEW(double *s, double *s_udg, double *pg, double *udg, double *udg_ref, appstruct &app, double *param,
             double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            source_ns2dNEW(s, s_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        case 3:
            source_ns3dNEW(s, s_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            break;
        default:
            exit(-1);
            break;
    }

    // Add inertial forces (if necessary)
    if (app.rotatingFrame == 1) {
        switch (nd) {
            case 2:
                source_inertial_ns2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                break;
            case 3:
                source_inertial_ns3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                break;
            default:
                exit(-1);
                break;
        }
    }
}
