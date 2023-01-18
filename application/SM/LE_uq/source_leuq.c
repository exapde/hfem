
#include "source_leuq2d.c"
#include "source_leuq3d.c"

void source_leuq(double *s, double *s_udg, double *pg, double *udg, double *udg_ref, appstruct &app, double *param,
             double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            source_leuq2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        case 3:
            source_leuq3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }

}
