
#include "flux_isotropicAV_ransSA2d.c"
#include "flux_isotropicAV_ransSA3d.c"
#include "flux_homogeneousAV_ransSA2d.c"
#include "flux_homogeneousAV_ransSA3d.c"
#include "flux_homogeneous2AV_ransSA2d.c"
#include "flux_homogeneous2AV_ransSA3d.c"

// Written by: C. Nguyen & P. Fernandez

void AVflux_ransSA(double *f, double *f_udg, double *pg, double *udg, double *udg_ref, appstruct &app, double *param,
                   double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    switch (nd) {
        case 2:
            switch (app.AVflag) {
                case 1:
                    flux_homogeneousAV_ransSA2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 2:
                    flux_homogeneous2AV_ransSA2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 3:
                    flux_isotropicAV_ransSA2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                default:
                    printf("Artificial viscosity model not implemented\n");
                    exit(-1);
            }
            break;
        case 3:
            switch (app.AVflag) {
                case 1:
                    flux_homogeneousAV_ransSA3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 2:
                    flux_homogeneous2AV_ransSA3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 3:
                    flux_isotropicAV_ransSA3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                default:
                    printf("Artificial viscosity model not implemented\n");
                    exit(-1);
            }
            break;
        default:
            exit(-1);
            break;
    }
}
