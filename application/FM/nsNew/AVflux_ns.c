
#include "flux_isotropicAV_ns2d.c"
#include "flux_isotropicAV_ns3d.c"
#include "flux_homogeneousAV_ns2d.c"
#include "flux_homogeneousAV_ns3d.c"
#include "flux_homogeneous2AV_ns2d.c"
#include "flux_homogeneous2AV_ns3d.c"
#include "flux_AV_ns2d.c"
#include "flux_AV_ns3d.c"
#include "flux_AV2_ns2d.c"
#include "flux_AVprecomputed_ns2d.c"
#include "flux_AVprecomputed_ns3d.c"
#include "flux_AVfrozen_ns2d.c"

// Written by: C. Nguyen & P. Fernandez

void AVflux_ns(double *f, double *f_udg, double *pg, double *udg, double *avField, meshstruct &mesh, masterstruct &master, appstruct &app, double *param,
               double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    // AVtype:
    //      0: Laplacian model with total enthalpy gradient for energy flux.
    //      1: Physical model with bulk viscosity only.
    //      2: Physical model with bulk viscosity and thermal conductivity (Pr^*=Pr).
    //      3: Physical model with thermal conductivity only.
    
    Int AVtype;
    
    switch (nd) {
        case 2:
            switch (app.AVflag) {
                case 1:
                    flux_homogeneousAV_ns2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 2:
                    flux_homogeneous2AV_ns2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 3:
                    flux_isotropicAV_ns2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 4:     // Sort of David's, but with less nonlinearity
                    flux_AV_ns2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 5:     // Cuong's approach with epsilon0
                    flux_AV2_ns2d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 6:     // p1 CG av field
//                     getAVfield(avField_DG, udg_ref, ng, nc, ncu, nd, ncd, app.AVflag);
//                     DG_2_p1CG(avField_p1CG, avField_DG, mesh, master, ndims);
                    AVtype = 0;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, avField, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                case 7:     // Like case 4, but with solution at previous DIRK stage for AV field
                    flux_AVfrozen_ns2d(f, f_udg, pg, udg, avField, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 8:     // p=1 CG AV with density smoothness based sensor (Per's approach). The solution at previous DIRK stage is used to compute AV field.
//                     getAVfield(avField_p1CG, udg_ref, ng, nc, ncu, nd, ncd, app.AVflag);
//                     flux_AV_rhoSmoothness_ns2d(f, f_udg, pg, udg, udg_ref, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    AVtype = 0;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, avField, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                case 9:     // Approach in note.tex
                    AVtype = 0;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, avField, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                case 10:
                    AVtype = 1;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, &avField[0*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    AVtype = 3;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, &avField[1*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    AVtype = 4;
                    flux_AVprecomputed_ns2d(f, f_udg, pg, &avField[2*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                default:
                    printf("Artificial viscosity model not implemented\n");
                    exit(-1);
            }
            break;
        case 3:
            switch (app.AVflag) {
                case 1:
                    flux_homogeneousAV_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 2:
                    flux_homogeneous2AV_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 3:
                    flux_isotropicAV_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 4:
                    flux_AV_ns3d(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
                    break;
                case 6:
                    error("AV 6 not available in 3D.\n");
                    // AVtype = 0;
                    // flux_AVprecomputed_ns3d(f, f_udg, pg, udg, udg_ref, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                case 9:     // Approach in note.tex
                    AVtype = 0;
                    flux_AVprecomputed_ns3d(f, f_udg, pg, avField, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    break;
                case 10:
                    AVtype = 1;
                    flux_AVprecomputed_ns3d(f, f_udg, pg, &avField[0*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    AVtype = 3;
                    flux_AVprecomputed_ns3d(f, f_udg, pg, &avField[1*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
                    AVtype = 4;
                    flux_AVprecomputed_ns3d(f, f_udg, pg, &avField[2*ng], udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian, AVtype);
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
