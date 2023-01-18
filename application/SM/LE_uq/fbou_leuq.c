
#include "fbou_leuq2d.c"
#include "fbou_leuq3d.c"

// void fbou_nsNEW(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *uhg, double *odg,
//           double *nl, double *ui, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
//           tempstruct &temp, Int ib, Int ife, Int ie, Int quadchoice, int computeJacobian, int numPoints)

void fbou_leuq(double *f, double *f_udg, double * f_uh, double *pg, double *udg, double *uhg, 
        double *odg, double *nl, double *ui, meshstruct &mesh, masterstruct &master, appstruct &app, 
        solstruct &sol, tempstruct &temp, Int ib, Int ife, Int ie, Int quadchoice, 
        int computeJacobian, int numPoints)
{    
    
    Int nfe, ncu, nch, ncq, nc, nco, nd, ncd;
    nd = master.nd;
    nfe = master.nfe;
    ncd = app.ncd;    
    nc = app.nc;
    ncu = app.ncu;
    ncq = app.ncq;
    nch = app.nch;
    nco = app.nco;
    double *param = &app.physicsparam[0];
    double time = app.time;
    
    switch (nd) {
        case 2:            
            fbou_leuq2d(f, f_udg, f_uh, pg, udg, uhg, nl, ui, param, time, ib, numPoints, nc, ncu, nd, ncd);            
            break;
        case 3:
            fbou_leuq3d(f, f_udg, f_uh, pg, udg, uhg, nl, ui, param, time, ib, numPoints, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }    
}


