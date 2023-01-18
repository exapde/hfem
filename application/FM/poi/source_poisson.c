
#include "source_poisson2d.c"
#include "source_poisson3d.c"

void source_poisson(double *s, double *s_udg, double *pg, double *udg, double *udg_ref, appstruct &app, double *param,
             double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            source_poisson2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        case 3:
            source_poisson3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }

}
