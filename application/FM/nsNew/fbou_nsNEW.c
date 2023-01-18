
#include "../fbou.cpp"

// Written by: C. Nguyen & P. Fernandez

// void fbou_nsNEW(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *udg_ref, double *uhg,
//                 double *uhg_ref, double *nl, double *avg_p1CG, double *ui, meshstruct &mesh,
//                 masterstruct &master, appstruct &app, double *param, double time, Int ib, Int* ndims,
//                 int computeJacobian, int numPoints)

//fbou_nsNEW(fh, fhn_u, fhn_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
void fbou_nsNEW(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *uhg, double *odg,
          double *nl, double *ui, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
          tempstruct &temp, Int ib, Int ife, Int ie, Int quadchoice, int computeJacobian, int numPoints)
{
    fbou(fh, fh_u, fh_uh, pg, udg, uhg, odg, nl, ui, mesh, master, app, sol, temp, ib, ife, ie, quadchoice, computeJacobian, numPoints);
}
