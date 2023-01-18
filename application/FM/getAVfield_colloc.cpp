#ifndef __GETAVFIELD
#define __GETAVFIELD

#include "../../util/DG_2_p1CG_colloc.cpp"

// Written by: C. Nguyen & P. Fernandez

void getAVfield_2d(double *avField, double *UDG, double *pn, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, int avgType)
{

    int i, j;

    int ncd = (int) ndims[1];
    int ne = (int) ndims[5];
    int npv = (int) ndims[9];
    int nc = (int) ndims[19];

    double *param = &app.param[0];

    double alpha = 100.0;
    double beta = 1.0e-2;
    double k_h = 1.5;

    double Minf = param[4];
    double porder = (double) app.porder;

    double pi = 3.141592653589793;

    double r, ru, rv, rE, rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u, v, ux, vy, h, f, x;

    double *x_p1CG, *x_DG, *f_p1CG, *f_DG;
    if (avgType == 0) {     // Sensor averaging
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 1) {        // Max surrogate averaging
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }

    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rE = UDG[i*nc*npv+3*npv+j];
            rx = UDG[i*nc*npv+4*npv+j];
            rux = UDG[i*nc*npv+5*npv+j];
            rvx = UDG[i*nc*npv+6*npv+j];
            rEx = UDG[i*nc*npv+7*npv+j];
            ry = UDG[i*nc*npv+8*npv+j];
            ruy = UDG[i*nc*npv+9*npv+j];
            rvy = UDG[i*nc*npv+10*npv+j];
            rEy = UDG[i*nc*npv+11*npv+j];

            u = ru/r;
            v = rv/r;
            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;

            h = pn[i*npv*ncd + 2*npv + j];

            if (avgType == 0) {     // Sensor averaging
                x_DG[i*npv+j] = (k_h*h/porder)*(ux+vy);
            }
            else if (avgType == 1) { // Max surrogate averaging
                x = (k_h*h/porder)*(ux+vy);
                f_DG[i*npv+j] = (x-beta)*(atan(alpha*(x-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        DG_2_p1CG(x_p1CG, x_DG, mesh, master, app, ndims);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + 2*npv + j];
                f = (x_p1CG[i*npv+j]-beta)*(atan(alpha*(x_p1CG[i*npv+j]-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/porder)*(1.0+1.0/Minf)*f;
            }
        }
    }
    else if (avgType == 1) { // Max surrogate averaging
        DG_2_p1CG(f_p1CG, f_DG, mesh, master, app, ndims);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + 2*npv + j];
                avField[i*npv+j] = (k_h*h/porder)*(1.0+1.0/Minf)*f_p1CG[i*npv+j];
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 1) {        // Max surrogate averaging
        delete[] f_p1CG; delete[] f_DG;
    }
}

void getAVfield_3d(double *avField, double *UDG, double *pn, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, int avgType)
{

    int i, j;

    int ncd = (int) ndims[1];
    int ne = (int) ndims[5];
    int npv = (int) ndims[9];
    int nc = (int) ndims[19];

    double *param = &app.param[0];

    double alpha = 100.0;
    double beta = 1.0e-2;
    double k_h = 1.5;

    double Minf = param[4];
    double porder = (double) app.porder;

    double pi = 3.141592653589793;

    double r, ru, rv, rw, rE, rx, rux, rvx, rwx, rEx, ry, ruy, rvy, rwy, rEy, rz, ruz, rvz, rwz, rEz, u, v, w, ux, vy, wz, h, f, x;

    double *x_p1CG, *x_DG, *f_p1CG, *f_DG;
    if (avgType == 0) {     // Sensor averaging
        x_p1CG = new double[ne*npv];
        x_DG = new double[ne*npv];
    }
    else if (avgType == 1) {        // Max surrogate averaging
        f_p1CG = new double[ne*npv];
        f_DG = new double[ne*npv];
    }

    for (i = 0; i < ne; i++) {
        for (j = 0; j < npv; j++) {
            r = UDG[i*nc*npv+0*npv+j];
            ru = UDG[i*nc*npv+1*npv+j];
            rv = UDG[i*nc*npv+2*npv+j];
            rw = UDG[i*nc*npv+3*npv+j];
            rE = UDG[i*nc*npv+4*npv+j];
            rx = UDG[i*nc*npv+5*npv+j];
            rux = UDG[i*nc*npv+6*npv+j];
            rvx = UDG[i*nc*npv+7*npv+j];
            rwx = UDG[i*nc*npv+8*npv+j];
            rEx = UDG[i*nc*npv+9*npv+j];
            ry = UDG[i*nc*npv+10*npv+j];
            ruy = UDG[i*nc*npv+11*npv+j];
            rvy = UDG[i*nc*npv+12*npv+j];
            rwy = UDG[i*nc*npv+13*npv+j];
            rEy = UDG[i*nc*npv+14*npv+j];
            rz = UDG[i*nc*npv+15*npv+j];
            ruz = UDG[i*nc*npv+16*npv+j];
            rvz = UDG[i*nc*npv+17*npv+j];
            rwz = UDG[i*nc*npv+18*npv+j];
            rEz = UDG[i*nc*npv+19*npv+j];

            u = ru/r;
            v = rv/r;
            w = rw/r;
            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;
            wz = (rwz - w*rz)/r;

            h = pn[i*npv*ncd + 3*npv + j];

            if (avgType == 0) {     // Sensor averaging
                x_DG[i*npv+j] = (k_h*h/porder)*(ux+vy+wz);
            }
            else if (avgType == 1) { // Max surrogate averaging
                x = (k_h*h/porder)*(ux+vy+wz);
                f_DG[i*npv+j] = (x-beta)*(atan(alpha*(x-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        DG_2_p1CG(x_p1CG, x_DG, mesh, master, app, ndims);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + 3*npv + j];
                f = (x_p1CG[i*npv+j]-beta)*(atan(alpha*(x_p1CG[i*npv+j]-beta))/pi + 1.0/2.0) - atan(alpha)/pi + 1.0/2.0;
                avField[i*npv+j] = (k_h*h/porder)*(1.0+1.0/Minf)*f;
            }
        }
    }
    else if (avgType == 1) { // Max surrogate averaging
        DG_2_p1CG(f_p1CG, f_DG, mesh, master, app, ndims);
        for (i = 0; i < ne; i++) {
            for (j = 0; j < npv; j++) {
                h = pn[i*npv*ncd + 3*npv + j];
                avField[i*npv+j] = (k_h*h/porder)*(1.0+1.0/Minf)*f_p1CG[i*npv+j];
            }
        }
    }

    if (avgType == 0) {     // Sensor averaging
        delete[] x_p1CG; delete[] x_DG;
    }
    else if (avgType == 1) {        // Max surrogate averaging
        delete[] f_p1CG; delete[] f_DG;
    }
}

void getAVfield(double *avField, double *UDG, double *pn, meshstruct &mesh, masterstruct &master, appstruct &app, Int *ndims, int avgType)
{

    int nd = (int) ndims[0];

    if (nd == 2) {
        getAVfield_2d(avField, UDG, pn, mesh, master, app, ndims, avgType);
    }
    else if (nd == 3) {
        getAVfield_3d(avField, UDG, pn, mesh, master, app, ndims, avgType);
    }
    else {
        error("Invalid number of dimensions.\n");
    }

}

#endif
