#ifndef __ELEMENTASSEMBLYSTANDARD
#define __ELEMENTASSEMBLYSTANDARD

#include "../application/fluxDriver.cpp"
#include "../application/fhatDriver.cpp"
#include "../application/sourceDriver.cpp"
#include "../application/fbouDriver.cpp"
// 
// Written by: C. Nguyen & P. Fernandez

void elementgeom(meshstruct &mesh, masterstruct &master, appstruct &app, tempstruct &temp, Int ie)
{
    /* This function computes the element mappings from the master element to the undeformed element */

    /* AT INPUT:
     * pn (npv / ncd): Nodal fields at element nodes.
     * shapvt (ngv / npv / (1+nd) (z.r)): Element shape functions and derivatives at element Gauss points.
     * ndims: Array of integers with information about the element and mesh.
     * */

    /* AT OUTPUT:
     * pg (ngv / ncd): Nodal fields at element Gauss points.
     * jacg (ngv): Determinant of mapping from master element to undeformed element evaluated at element Gauss points.
     * Jg (ngv / nd (z.x) / nd (z.r)): Jacobian matrix of mapping from master element to undeformed element evaluated at element Gauss points.
     * Xxg (ngv / nd (z.r) / nd (z.x)): Minus inverse times determinant of mapping from master element to undeformed element evaluated at element Gauss points.
     * */

    Int i, iA, iB;
    char chn = 'N';
    double one = 1.0, zero = 0.0;

    /* Get dimensions */
    Int nd, npv, ngv, ncd;
    nd = master.nd;    
    npv = master.npv;    
    ncd = app.ncd;

    double *shapvt, *pg, *Jg, *jac, *Xx, *pn;        
    pn = &mesh.dgnodes[ie*npv*ncd];    
    ngv = master.ngv;
    shapvt = &master.shapvt[0];
    pg = &temp.pg[0];
    Jg = &temp.Jg[0];
    jac = &temp.jacg[0];
    Xx = &temp.Xxg[0];                        
    
    DGEMM(&chn, &chn, &ngv, &ncd, &npv, &one, shapvt, &ngv, pn,
          &npv, &zero, pg, &ngv);
    /* shapvt: ngv / npv / (1+nd) (z.r) (we only use ngv / npv here) */
    /* pn: npv / ncd */
    /* pg: ngv / ncd */

    /* Compute Jacobian matrix at element Gauss points */
    for (i=0; i<nd; i++) {
        iA = (i+1)*ngv*npv;
        iB = i*ngv*nd;
        DGEMM(&chn, &chn, &ngv, &nd, &npv, &one, &shapvt[iA], &ngv, pn,
              &npv, &zero, &Jg[iB], &ngv);
    }
    /* shapvt: ngv / npv / (1+nd) (z.r) (we only use ngv / npv / nd (z.r) here) */
    /* pn: npv / ncd (we only use npv / nd (z.x) here) */
    /* Jg: ngv / nd (z.x) / nd (z.r) */

    double *Jg11, *Jg12, *Jg13, *Jg21, *Jg22, *Jg23, *Jg31, *Jg32, *Jg33;
    double *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;

    /* Compute determinant and inverse times Jacobian of mapping from master element to undeformed element evaluated at element Gauss points */
    if (nd==1) {
        for (i=0; i<ngv; i++) {
            jac[i] = Jg[i];
            Xx[i] = 1.0;
        }
    }
    else if (nd==2) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg21 = &Jg[2*ngv];
        Jg22 = &Jg[3*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx12 = &Xx[2*ngv];
        Xx22 = &Xx[3*ngv];
        for (i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
            Xx11[i] = Jg22[i];
            Xx21[i] = -Jg21[i];
            Xx12[i] = -Jg12[i];
            Xx22[i] = Jg11[i];
        }
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg13 = &Jg[2*ngv];
        Jg21 = &Jg[3*ngv];
        Jg22 = &Jg[4*ngv];
        Jg23 = &Jg[5*ngv];
        Jg31 = &Jg[6*ngv];
        Jg32 = &Jg[7*ngv];
        Jg33 = &Jg[8*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx31 = &Xx[2*ngv];
        Xx12 = &Xx[3*ngv];
        Xx22 = &Xx[4*ngv];
        Xx32 = &Xx[5*ngv];
        Xx13 = &Xx[6*ngv];
        Xx23 = &Xx[7*ngv];
        Xx33 = &Xx[8*ngv];
        for (i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i]*Jg33[i] - Jg11[i]*Jg32[i]*Jg23[i] +
                     Jg21[i]*Jg32[i]*Jg13[i] - Jg21[i]*Jg12[i]*Jg33[i] +
                     Jg31[i]*Jg12[i]*Jg23[i] - Jg31[i]*Jg22[i]*Jg13[i];
            Xx11[i] = Jg22[i]*Jg33[i] - Jg23[i]*Jg32[i];
            Xx21[i] = Jg23[i]*Jg31[i] - Jg21[i]*Jg33[i];
            Xx31[i] = Jg21[i]*Jg32[i] - Jg22[i]*Jg31[i];
            Xx12[i] = Jg13[i]*Jg32[i] - Jg12[i]*Jg33[i];
            Xx22[i] = Jg11[i]*Jg33[i] - Jg13[i]*Jg31[i];
            Xx32[i] = Jg12[i]*Jg31[i] - Jg11[i]*Jg32[i];
            Xx13[i] = Jg12[i]*Jg23[i] - Jg13[i]*Jg22[i];
            Xx23[i] = Jg13[i]*Jg21[i] - Jg11[i]*Jg23[i];
            Xx33[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
        }
    }
    /* jacg: ngv */
    /* Xxg: ngv / nd (z.r) / nd (z.x) */
}

void facegeom(meshstruct &mesh, masterstruct &master, appstruct &app, tempstruct &temp, Int ie)
{
    /* This function computes the face mappings from the master face to the undeformed face */

    /* AT INPUT::
     * pn (npv / ncd): Nodal fields at element nodes.
     * shapft (ngf / npf / nd (z.r)): Face shape functions and derivatives at face Gauss points.
     * perm (npf / nfe): Mapping from face nodes to element nodes.
     * ndims: Array of integers with information about the element and mesh */

    /* AT OUTPUT:
     * pg (ngf / nfe / ncd): Nodal fields at face Gauss points.
     * Jg (ngf / nfe / nd (z.x) / nd-1 (z.r). Except in 1D, that is length 2): Jacobian matrix of the face mapping at face Gauss points.
     * jacg (ngf / nfe. Except in 1D, that is length 1): Determinant of the face mapping at face Gauss points.
     * nlg (ngf / nfe / nd (z.x)): Normal vector at face Gauss points.
     * */

    Int inc = 1, i, j, k, iA, iB, is, na, nb;
    char chn = 'N';
    double one = 1.0, zero = 0.0;

    /* Get dimensions */
    Int nd, npv, npf, ngf, nfe, ncd, ndf;
    ncd = app.ncd;
    nd = master.nd;    
    nfe = master.nfe;
    npf = master.npf[0];
    ndf = master.ndf;
    ngf = master.ngf[0];
    npv = master.npv;
    
    double *shapft, *pf, *pg, *Jg, *jacg, *nlg, *pn;
    shapft = &master.shapft[0][0];
    pg = &temp.pgf[0];
    pf = &temp.pf[0];
    Jg = &temp.Jgf[0];
    jacg = &temp.jacgf[0];
    nlg = &temp.nlgf[0];
    pn = &mesh.dgnodes[ie*npv*ncd];    

    /* Get nodal fields at face nodes */
    for (j=0; j<ncd; j++)
        for (is=0; is<nfe; is++)
            for (k=0; k<npf; k++)
                pf[j*ndf+is*npf+k] = pn[j*npv+master.permgeom[is][k]];
    /* pn: npv / ncd */
    /* pf: npf / nfe / ncd */

    /* Compute nodal fields at face Gauss points */
    na = nfe*ncd;
    DGEMM(&chn, &chn, &ngf, &na, &npf, &one, shapft, &ngf, pf,
          &npf, &zero, pg, &ngf);
    /* shapft: ngf / npf / nd (z.r) (only ngf / npf is used here) */
    /* pf: npf / nfe / ncd */
    /* pg: ngf / nfe / ncd */

    /* Compute Jacobian matrix at face Gauss points */
    nb = nfe*nd;
    for (i=0; i<nd-1; i++) {
        iA = (i+1)*ngf*npf;
        iB = i*ngf*nfe*nd;
        DGEMM(&chn, &chn, &ngf, &nb, &npf, &one, &shapft[iA], &ngf, pf,
              &npf, &zero, &Jg[iB], &ngf);
    }
    /* shapft: ngf / npf / nd (z.r) (only ngf / npf / nd-1 is used here) */
    /* pf: npf / nfe / ncd (only npf / nfe / nd (z.x) is used here) /*/
    /* Jg: ngf / nfe / nd (z.x) / nd-1 (z.r) */

    double *Jg11, *Jg12, *Jg21, *Jg22, *Jg31, *Jg32;

    /* Compute determinant and normal vector at face Gauss points */
    na = ngf*nfe;
    if (nd==1) {
        jacg[0] = 1.0;
        nlg[0] = 1.0;
        nlg[1] = -1.0;
    }
    else if (nd==2) {
        for (i=0; i<na; i++) {
            j = i+na;
            jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
            nlg[i] = Jg[j]/jacg[i];
            nlg[j] = -Jg[i]/jacg[i];
        }
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg21 = &Jg[na];
        Jg31 = &Jg[2*na];
        Jg12 = &Jg[3*na];
        Jg22 = &Jg[4*na];
        Jg32 = &Jg[5*na];
        for (i=0; i<na; i++) {
            j = i+na;
            k = i+2*na;
            nlg[i] = Jg21[i]*Jg32[i] - Jg31[i]*Jg22[i];
            nlg[j] = Jg31[i]*Jg12[i] - Jg11[i]*Jg32[i];
            nlg[k] = Jg11[i]*Jg22[i] - Jg21[i]*Jg12[i];
            jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
            nlg[i] = nlg[i]/jacg[i];
            nlg[j] = nlg[j]/jacg[i];
            nlg[k] = nlg[k]/jacg[i];
        }
    }
    /* Jg: ngf / nfe / nd (z.x) / nd-1 (z.r). Except in 1D, that is length 2. */
    /* Nfg: ngf / nfe / nd (z.x) */
    /* jacg: ngf / nfe. Except in 1D, that is length 1. */
    /* pg: ngf / nfe / ncd */
}

void elementflux(meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        Int ie, int computeJacobian)
{

    Int inc = 1, i, j, m, n, k, oneInt = 1;
    char chn = 'N';
    double one = 1.0, zero = 0.0;

    /* Get dimensions */    
    Int npv = master.npv;    
    Int ngv = master.ngv;    
    Int nc  = app.nc;
    Int ncu = app.ncu;
    Int nch = app.nch;   
    Int ncd = app.ncd;
    Int nco = app.nco;

    double *shapvt, *s, *f, *s_udg, *f_udg, *pg, *udg, *odg, *udgg, *odgg, *sh;
    shapvt = &master.shapvt[0];
    s = &temp.s[0];
    f = &temp.f[0];
    s_udg = &temp.s_udg[0];
    f_udg = &temp.f_udg[0];
    pg = &temp.pg[0];
    udgg = &temp.udgg[0];
    udg = &sol.UDG[ie*npv*nc];
    sh = &sol.SH[ie*npv*nc];   
    
    /* Compute solution at Gauss points */
    DGEMM(&chn, &chn, &ngv, &nc, &npv, &one, shapvt, &ngv, udg,
                    &npv, &zero, udgg, &ngv);

    /* Compute other fields at Gauss points */
    if (nco>0) {                 
        odg = &sol.ODG[ie*npv*nco];        
        odgg = &temp.odgg[0];
        DGEMM(&chn, &chn, &ngv, &nco, &npv, &one, shapvt, &ngv, odg,
                        &npv, &zero, odgg, &ngv);
    }
        
    /* Compute flux and source at Gauss points */
    //flux(f, f_udg, pg, udgg, param, time, appname, ndims);
    //source(s, s_udg, pg, udgg, param, time, appname, ndims);
    Int quadchoice = 1;
    fluxDriver(f, f_udg, pg, udgg, odgg, mesh, master, app, sol, temp, ie, quadchoice, computeJacobian, ngv);
    sourceDriver(s, s_udg, pg, udgg, odgg, mesh, master, app, sol, temp, ie, quadchoice, computeJacobian, ngv);
    
// void fluxDriver(double * f,double * f_UDG, double * pg, double * UDG, double * UDG_ref, double *avField_p1CG, 
//         meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
//         int computeJacobian, int numPoints)
    
// void fluxDriver(double * f,double * f_UDG, double * pg, double * UDG, double * UDG_ref, double *avField_p1CG, meshstruct &mesh, masterstruct &master, appstruct &app, double * param, double time,
//                 int computeJacobian, int numPoints)
                
    if (app.tdep==1) {
        /* Compute time source at Gauss points */
        DGEMM(&chn, &chn, &ngv, &ncu, &npv, &one, shapvt, &ngv, sh,
                    &npv, &one, s, &ngv);

        for (j=0; j<ncu; j++)
            for (i=0; i<ngv; i++) {
                s[j*ngv+i] -= udgg[j*ngv+i]*app.fc_u;
                if (computeJacobian == 1)
                    s_udg[(j*ncu+j)*ngv+i] -= app.fc_u;
            }
    }
}

void elementint(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app,
                solstruct &sol, tempstruct &temp, double* elemtimes, Int ie, int computeJacobian)
{
    /* This function computes the contribution to the matrices B and D, and to the residual Ru due to element integrals */

    /* AT INPUT:
     * f (ngv / ncu / nd): Fluxes at element Gauss points.
     * f_UDG (ngv / ncu / nd / nc): Derivatives of fluxes at element Gauss points with respect to conservative variables.
     * s (ngv / ncu): Source term at element nodes in undeformed domain.
     * s_UDG (ngv / ncu / nc): Derivative of source term at element nodes in undeformed domain.
     * jacg (ngv): Determinant of mapping from master element to undeformed element evaluated at element Gauss points.
     * Xxg (ngv / nd (z.r) / nd (z.x)): Minus inverse times determinant of mapping from master element to undeformed element evaluated at element Gauss points.
     * SHAPVG (npv / ngv / (1+nd) (z.r)): Element shape functions and derivatives times Gauss weights at element Gauss points.
     * SHAPVGDOTSHAPVL (npv / npv / ngv / (1+nd) (z.r)): Product of two element shape functions times Gauss weights at element Gauss points.
     * prob:
     * FLAGS: Array of inteters with different flags.
     * NDIMS: Array of integers with information about the element and mesh.
     * /*

    /* AT OUTPUT:
     * Ru (npv / ncu): Residual of u.
     * BD (npv / npv / ncu / nc): Jacobian of residual of u w.r.t. to u (D) and q (B).
     * */

    Int inc = 1, d, i, j, m, n, k, t00, t0, t03, t1, t10, t100, t2, t3, t30, iA, iB, iC, info;
    char chn = 'N', charu = 'U';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    clock_t t;

    /* Get dimensions */    
    Int npv, ngv, ncd, nc, ncu, ncq, nch, nd, nd1, ndf;
    nd = master.nd;        
    npv = master.npv;    
    ngv = master.ngv;
    ndf = master.ndf;
    ncd = app.ncd;
    nc  = app.nc;
    ncu = app.ncu;
    nch = app.nch;    
    ncq = app.ncq;
    nd1 = nd+1;    

    // The cost of elementgeom is non-neglibile for p<=2 (2D tri).
    t = clock();
    if (app.flag_q==0)
        elementgeom(mesh, master, app, temp, ie);
    elemtimes[3] += clock() - t;   
    
    // The cost of elementgeom is non-neglibile for p<=3 (2D tri).
    t = clock();
    elementflux(mesh, master, app, sol, temp, ie, 1);
    elemtimes[4] += clock() - t;
    
    double *shapvt, *shapvg, *shapvgdotshapvl, *Rq, *M, *C;
    double *s, *f, *s_udg, *f_udg, *jac, *Xx, *wrk, *wrl, *Ru, *BDt; // *BD;
    shapvt = &master.shapvt[0];
    shapvg = &master.shapvg[0];
    shapvgdotshapvl = &master.shapvgdotshapvl[0];
    s = &temp.s[0];
    f = &temp.f[0];
    s_udg = &temp.s_udg[0];
    f_udg = &temp.f_udg[0];
    jac = &temp.jacg[0];
    Xx = &temp.Xxg[0];
    wrk = &temp.wrk[0];
    wrl = &temp.wrl[0];
    M = &elem.M[0];
    C  = &elem.C[0];
    Ru = &elem.Ru[0];
    Rq = &elem.Rq[0];
    BDt = &temp.BDt[0];

    Int sz[10];
    sz[0] = ngv*nd;
    sz[1] = ngv*nc;
    sz[2] = ngv*ncu;
    sz[3] = ngv*nd1;
    sz[4] = ngv*ncu*nd;
    sz[5] = ngv*ncu*nc;
    sz[6] = ngv*nd1*ncu;
    sz[7] = ngv*ncu*nd*nc;
    sz[8] = ngv*nd1*ncu*nc;

    /* Compute wrk */
    for (j=0; j<ncu; j++)
        for (i=0; i<ngv; i++)
            wrk[j*sz[3]+0*ngv+i] = s[j*ngv+i]*jac[i];

    // Note: The first implementation is slightly faster for low-mid p, and similar performance for high p (g++ -O3)
    // Since the cost of this operation is only non-negligible w.r.t. Ru computation for low-to-mid p, the first implementation is chosen.
    for (m=0; m<nd; m++) {
        for (j=0; j<ncu; j++)
            for (i=0; i<ngv; i++)
                wrk[j*sz[3]+(m+1)*ngv+i] = f[0*sz[2]+j*ngv+i]*Xx[m*sz[0]+0*ngv+i];
        for (n=1; n<nd; n++)
            for (j=0; j<ncu; j++)
                for (i=0; i<ngv; i++)
                    wrk[j*sz[3]+(m+1)*ngv+i] += f[n*sz[2]+j*ngv+i]*Xx[m*sz[0]+n*ngv+i];
    }

    /* Compute the residual vector */
    DGEMM(&chn, &chn, &npv, &ncu, &sz[3], &one, shapvg, &npv, wrk,
                    &sz[3], &zero, Ru, &npv);

    // wrl: ngv*(nd+1)*ncu*nc
    // Note: No difference between both implementations (g++ -O3)
    // At most, the first is slightly faster for low p. Since the cost of this operation is non-negligible only for low p, the first implementation is used.

    t = clock();
    /* Compute wrl */
    for (k=0; k<nc; k++)
        for (j=0; j<ncu; j++)
            for (i=0; i<ngv; i++)
                wrl[k*sz[6]+j*sz[3]+0*ngv+i] = -s_udg[k*sz[2]+j*ngv+i]*jac[i];

    // Note: For p>=4 the second implementation is faster. For p=1 the first is faster (tested with g++ -O3)
    // Since this operations only is non-negligible w.r.t. BD computation for p=1, the first implementation is chosen.
    for (m=0; m<nd; m++) {
        for (k=0; k<nc; k++)
            for (j=0; j<ncu; j++)
                for (i=0; i<ngv; i++)
                    wrl[k*sz[6]+j*sz[3]+(m+1)*ngv+i] =
                            -f_udg[k*sz[4]+0*sz[2]+j*ngv+i]*Xx[m*sz[0]+0*ngv+i];
        for (n=1; n<nd; n++)
            for (k=0; k<nc; k++)
                for (j=0; j<ncu; j++)
                    for (i=0; i<ngv; i++)
                        wrl[k*sz[6]+j*sz[3]+(m+1)*ngv+i] -= f_udg[k*sz[4]+n*sz[2]+j*ngv+i]*Xx[m*sz[0]+n*ngv+i];
    }
    elemtimes[5] += clock() - t;
    
    //vector<double> wrlshap(ngv*nd1*ncu*npv*ncu);
    t = clock();
    double *wrlshap = &temp.wrlshapMiC[0];
    sz[0] = ngv*nd1;
    sz[1] = sz[0]*ncu;
    sz[2] = sz[1]*npv;        
    for (k=0; k<ncu; k++)
        for (n=0; n<npv; n++)
            for (j=0; j<ncu; j++)
                for (m=0; m<nd1; m++)
                    for (i=0; i<ngv; i++)
                        wrlshap[k*sz[2]+n*sz[1]+j*sz[0]+m*ngv+i] = wrl[k*sz[1]+j*sz[0]+m*ngv+i]*shapvt[n*ngv+i];
    elemtimes[6] += clock() - t;

    if (app.flag_q==1) {
        double *M = &elem.M[0];
        double *MiC = &elem.C[0];
        double *MiE = &elem.E[0];
        double *shapMiC = &temp.shapMiC[0];
        double *shapMiE = &temp.shapMiE[0];
        double *wrlshapMiE = &temp.wrlshapMiE[0];

        t = clock();
        /* Factor M using Cholesky decomposition */        
        DPOTRF(&charu,&npv,&M[0],&npv,&info);
                        
        sz[1] = npv*nd;
        /* compute inv(M)*C and store it in C */        
        DPOTRS(&charu,&npv,&sz[1],&M[0],&npv,&MiC[0],&npv,&info);        
                
        sz[4] = npv*nd;
        DGEMM(&chn, &chn, &ngv, &sz[4], &npv, &one, shapvt, &ngv, &MiC[0],
                    &npv, &zero, &shapMiC[0], &ngv);                
        
        sz[2] = ndf*nd;
        DPOTRS(&charu,&npv,&sz[2],&M[0],&npv,&MiE[0],&npv,&info);
        
        sz[4] = ndf*nd;
        DGEMM(&chn, &chn, &ngv, &sz[4], &npv, &minusone, shapvt, &ngv, &MiE[0],
                    &npv, &zero, &shapMiE[0], &ngv);
        elemtimes[7] += clock() - t;


        t = clock();
        sz[0] = ngv*nd1;
        sz[1] = sz[0]*ncu;
        sz[2] = sz[1]*npv;                        
        sz[3] = sz[1]*ncu;
        sz[4] = ngv*npv;
        sz[5] = ngv*nd1;
        sz[6] = sz[0]*ncu;
        sz[7] = sz[1]*ncu;                        
        sz[8] = sz[1]*ncu;
        sz[9] = ngv*ndf;        
        if (nd==1) {
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshap[k*sz[2]+n*sz[1]+j*sz[0]+m*ngv+i] += 
                                    wrl[sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[n*ngv+i];                                    
            
            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)    
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshapMiE[n*sz[7]+k*sz[6]+j*sz[5]+m*ngv+i] = 
                                    wrl[sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[n*ngv+i];            
        }
        if (nd==2) {                    
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshap[k*sz[2]+n*sz[1]+j*sz[0]+m*ngv+i] += 
                                    wrl[sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[n*ngv+i]+                                        
                                    wrl[2*sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[sz[4]+n*ngv+i];                                    
            
            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)    
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshapMiE[n*sz[7]+k*sz[6]+j*sz[5]+m*ngv+i] = 
                                    wrl[(0+1)*sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[0*sz[9]+n*ngv+i]+
                                    wrl[(1+1)*sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[1*sz[9]+n*ngv+i];                                    
        }                    
        if (nd==3) {   
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshap[k*sz[2]+n*sz[1]+j*sz[0]+m*ngv+i] += 
                                    wrl[sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[n*ngv+i]+                                        
                                    wrl[2*sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[sz[4]+n*ngv+i]+
                                    wrl[3*sz[3]+k*sz[1]+j*sz[0]+m*ngv+i]*shapMiC[2*sz[4]+n*ngv+i];                                    
            
            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)    
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nd1; m++)
                            for (i=0; i<ngv; i++)                                                         
                                wrlshapMiE[n*sz[7]+k*sz[6]+j*sz[5]+m*ngv+i] = 
                                    wrl[(0+1)*sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[0*sz[9]+n*ngv+i]+
                                    wrl[(1+1)*sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[1*sz[9]+n*ngv+i]+
                                    wrl[(2+1)*sz[8]+k*sz[6]+j*sz[5]+m*ngv+i]*shapMiE[2*sz[9]+n*ngv+i];                                    
        }
        elemtimes[8] += clock() - t;

        t = clock();
        sz[3] = ncu*ncu*ndf;
        sz[0] = ngv*nd1;
        /* Compute the matrix F */
        DGEMM(&chn, &chn, &npv, &sz[3], &sz[0], &one, &shapvg[0], &npv, &wrlshapMiE[0],
                    &sz[0], &zero, &elem.F[0], &npv);
        elemtimes[9] += clock() - t;
    }
    else if (app.flag_q==0) {
        for (i=0; i<npv*ncu*ncu*ndf; i++)
            elem.F[i] = 0.0;
    }

    /* Compute the matrix D */
    t = clock();
    sz[3] = ncu*npv*ncu;
    sz[0] = ngv*nd1;
    DGEMM(&chn, &chn, &npv, &sz[3], &sz[0], &one, shapvg, &npv, &wrlshap[0],
                &sz[0], &zero, &elem.BD[0], &npv);
    elemtimes[10] += clock() - t;
}


void faceflux(meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        Int ie, int computeJacobian)
{
    Int inc = 1, i, j, m, n, k, is, ks, na, nb;
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    /* Get dimensions */
    Int npv = master.npv;    
    Int ndf = master.ndf;
    Int npf = master.npf[0];
    Int ngf = master.ngf[0];
    Int nfe = master.nfe;
    Int nc  = app.nc;
    Int nch = app.nch;
    Int nco = app.nco;
        
    double *shapft, *fh, *fh_u, *fh_uh, *udg, *odg, *uf, *of, *ugf, *ogf, *uh, *uhg, *nlg, *pgf;
    shapft = &master.shapft[0][0];
    fh = &temp.fh[0];
    fh_u = &temp.fh_u[0];
    fh_uh = &temp.fh_uh[0];                
    udg = &sol.UDG[ie*npv*nc];        
    uf = &temp.uf[0];
    uh = &temp.uh[0];            
    ugf = &temp.ugf[0];
    uhg = &temp.uhg[0];
    pgf = &temp.pgf[0];
    nlg = &temp.nlgf[0];                
    
    for (j=0; j<ndf; j++) {
        n = mesh.elcon[ie*ndf+j];
        for (k=0; k<nch; k++)
            uh[k*ndf+j] = sol.UH[n*nch+k];
    }
    
//     if (app.my_rank==0) {
//         if (ie==0) {
//             print2darray(uh,ndf,nch);
//             print1iarray(&mesh.elcon[ie*ndf],ndf);
//             error("here");
//         }
//     }

    /* Compute solution on the face */
    for (j=0; j<nc; j++)
        for (is=0; is<nfe; is++)
            for (k=0; k<npf; k++) 
                uf[j*ndf+is*npf+k] = udg[j*npv+master.perm[is][k]];                            

//     cout<<nco;        
//     print2darray(shapft, ngf, npf);
//     print2darray(uf, ngf, nc);
//     error("here");
    
    na = nfe*nc;
    /* Compute solution at Gauss points */
    DGEMM(&chn, &chn, &ngf, &na, &npf, &one, shapft, &ngf, uf,
                    &npf, &zero, ugf, &ngf);        
    
    if (nco>0) {
        odg = &sol.ODG[ie*npv*nco];
        of = &temp.of[0];        
        ogf = &temp.ogf[0];
        for (j=0; j<nco; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<npf; k++) 
                    of[j*ndf+is*npf+k] = odg[j*npv+master.perm[is][k]];            
        
        na = nfe*nco;
        /* Compute solution at Gauss points */
        DGEMM(&chn, &chn, &ngf, &na, &npf, &one, shapft, &ngf, of,
                        &npf, &zero, ogf, &ngf);        
    }
    
        
    nb = nfe*nch;
    /* Compute uh at Gauss points */
    DGEMM(&chn, &chn, &ngf, &nb, &npf, &one, shapft, &ngf, uh,
                    &npf, &zero, uhg, &ngf);

    /* Compute flux and source at Gauss points */
    //fhat(fh, fh_u, fh_uh, pgf, ugf, uhg, nlg, param, time, appname, ndims);
    Int fhatExpression = 0;
    Int quadchoice = 1;
    fhatDriver(fh, fh_u, fh_uh, pgf, ugf, uhg, ogf, nlg, mesh, master, app, sol, temp, 
            fhatExpression, ie, quadchoice, computeJacobian, ngf*nfe);     
    
// void fhatDriver(double * fhn, double * fhn_UDG, double * fhn_UH, double * pg, double * UDG, double * UH, double * ODG,
//                 double * NL,  meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
//                 tempstruct &temp, Int fhatExpression, Int ie, Int quadchoice, Int computeJacobian, Int numPoints)
    
//     if (app.my_rank==0) {         
//         Int nd = app.nd;
//         print2darray(pgf,ngf*nfe,nd);
//         print2darray(nlg,ngf*nfe,nd);
//         print2darray(ugf,ngf*nfe,nc);
//         print2darray(uhg,ngf*nfe,nch);
//         print2darray(fh,ngf*nfe,nch);
//         error("here");
//     }
    
}

void bouflux(meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp,
         Int ib, Int bn, Int ife, Int ie, int computeJacobian)
{    
    double *fb, *fb_u, *fb_uh, *uft, *oft, *uht, *nlt, *pft;    
    fb = &temp.fb[0];
    fb_u = &temp.fb_u[0];
    fb_uh = &temp.fb_uh[0];
    uft = &temp.uft[0];
    oft = &temp.oft[0];
    uht = &temp.uht[0];
    pft = &temp.pft[0];
    nlt = &temp.nlt[0];                
    
    Int ngf = master.ngf[0];
    double *bcs   = &app.bcs[ib*app.nch];
    
    // bn: boundary condition index for the ib_th boundary
    // ib: boundary    
    Int quadchoice = 1;
    fbouDriver(fb, fb_u, fb_uh, pft, uft, uht, oft, nlt, bcs, mesh, master, app, sol, temp,  
               bn, ife, ie, quadchoice, computeJacobian, ngf);
    
//     if (ie==43) {
//         if (app.my_rank==1) {
//             cout<<app.my_rank<<" . "<<ib<<"  "<<bn<<"  "<<ife<<"  "<<ie<<endl;        
//             print2darray(pft,ngf,app.nd);
//             print2darray(uft,ngf,app.nc);
//             print2darray(uht,ngf,app.nch);        
//             print2darray(nlt,ngf,app.nd);
//             print1darray(bcs,app.nch);
//             print2darray(fb,ngf,app.nch);
//         }
//     }
    
    //fbouDriver(fb, fb_u, fb_uh, pft, uft, uft_ref, uht, uht_ref, nlt, avft, bcs, mesh, master, app, param, time, bn, ndims, computeJacobian, (int) numPoints);
    //fhatDriver(fh, fh_u, fh_uh, pgf, ugf, uhg, ogf, nlg, mesh, master, app, sol, temp, fhatExpression, computeJacobian, ngf);       
}

void faceint(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app,
             solstruct &sol, tempstruct &temp, double* elemtimes, Int ie, int computeJacobian)
{
    /* This function computes the contribution to the residual and Jacobian of u and uh due to face integrals */
    /* It takes Ru and BD due to the element integrals and yields the final expression for Ru, BD, F, Rh, GK and H */

    /* AT INPUT:
     * Ru (npv / ncu): Residuals of u.
     * BD (npv / npv / ncu / nc): Jacobian of residual of u w.r.t. to u (D) and q (B).
     * fhn (ngf / nfe / ncu): Normal component of flux trace at face Gauss points.
     * fhn_U (ngf / nfe / ncu / nc): Derivative of normal component of flux trace w.r.t. conserved variables at face Gauss points.
     * fhn_uh (ngf / nfe / ncu / nch): Derivative of normal component of flux trace w.r.t. uh at face Gauss points.
     * jacgf (ngf / nfe. Except in 1D, that is length 1): Determinant of the face mapping at face Gauss points.
     * Ngf (ngf / nfe / nd (z.x)): Normal vector at face Gauss points.
     * Pgf (ngf / nfe / ncd): Nodal fields at face Gauss points.
     * SHAPFG (npf / ngf / nd (z.r)): Face shape functions and derivatives times Gauss weights at face Gauss points.
     * SHAPFGDOTSHAPFC (npf / npf / ngf / nd (z.r)): Product of two face shape functions times Gauss weights at face Gauss points.
     * PERM (npf / nfe): Mapping from face nodes to element nodes.
     * BF (nfe):
     * PARAM:
     * BCM (numBoundaries): Type of boundary condition to be imposed.
     * BCS (ncu / numBoundaries): Conservative variables (state) associated to the boundary.
     * time:
     * prob:
     * FLAGS: Array of inteters with different flags.
     * NDIMS: Array of integers with information about the element and mesh.
     * /*

    /* AT OUTPUT:
     * Ru (npv / ncu): Residuals of u.
     * BD (npv / ncu / ncu / npv / (1+nd) [uq solver] OR npv / ncu / ncu / npv  [u solver]): Jacobian of residual of u w.r.t. to u (D) and q (B).
     * F (npv / ncu / nch / ndf): Jacobian of residual of u w.r.t. to uh.
     * Rh (nch / npf / nfe): Residuals of uh.
     * GK (nch / npf / nfe / npv / nc): Jacobian of residual of uh w.r.t. to u (K) and q (G).
     * H (nch / npf / nfe / nch / npf / nfe): Jacobian of residual of uh w.r.t. to uh.
     * */

    Int inc = 1, d, i, j, m, n, k, is, ks, t00, t0, t1, t2, na, nb, iA, iB, iC, info;
    Int ifhujac, iBDtmp, ifhuhjac, iFtmp, ifhjac, iRutmp;
    Int sz[10];
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    clock_t t;

    Int npv, ngv, ncd, nc, ncu, nch, nco, npf, nfe, ngf, nd, ndf, nd1, sza, szb;
    nd = app.nd;
    ncd = app.ncd;
    nfe = master.nfe;
    npv = master.npv;
    npf = master.npf[0];
    ngf = master.ngf[0];
    ngv = master.ngv;
    nc  = app.nc;
    ncu = app.ncu;
    nch = app.nch;
    nco = app.nco;
    ndf = npf*nfe;
    nd1 = nd+1;
    sza = ngf*nfe;
    szb = npf*npf;

    t = clock();
    if (app.flag_q==0)
        facegeom(mesh, master, app, temp, ie);
    elemtimes[11] += clock() - t;
    
    t = clock();
    faceflux(mesh, master, app, sol, temp, ie, 1);
    elemtimes[12] += clock() - t;

    double *shapft, *shapfg, *fh, *fh_u, *fh_uh, *ugf, *ogf, *uhg, *nlg, *pgf, *jacf, *nlgjac, *param;
    double* shapfgdotshapfc = &master.shapfgdotshapfc[0][0];
    shapfg = &master.shapfg[0][0];
    shapft = &master.shapft[0][0];
    fh = &temp.fh[0];
    fh_u = &temp.fh_u[0];
    fh_uh = &temp.fh_uh[0];
    uhg = &temp.uhg[0];
    ugf = &temp.ugf[0];
    ogf = &temp.ogf[0];
    pgf = &temp.pgf[0];
    nlg = &temp.nlgf[0];
    jacf = &temp.jacgf[0];
    nlgjac = &temp.nlgjac[0];

    double *E, *Rq, *Ru, *BD, *F, *GK, *H, *Rh, *Rutmp, *BDtmp, *Ftmp, *Etmp, *BDt;
    double *fhjac, *fhujac, *fhuhjac;
    E  = &elem.E[0];
    Rq = &elem.Rq[0];
    Ru = &elem.Ru[0];
    BD = &elem.BD[0];
    F = &elem.F[0];
    GK = &elem.GK[0];
    H = &elem.H[0];
    Rh = &elem.Rh[0];
    Rutmp = &temp.Rutmp[0];
    BDtmp = &temp.BDtmp[0];
    BDt = &temp.BDt[0];
    Ftmp = &temp.Ftmp[0];
    Etmp = &temp.Etmp[0];
    fhjac = &fh[0];
    fhujac = &fh_u[0];
    fhuhjac = &fh_uh[0];

    double *fb, *fb_u, *fb_uh, *uft, *oft, *uht, *nlt, *pft;
    fb = &temp.fb[0];
    fb_u = &temp.fb_u[0];
    fb_uh = &temp.fb_uh[0];
    uft = &temp.uft[0];
    uht = &temp.uht[0];
    pft = &temp.pft[0];
    nlt = &temp.nlt[0];

    //udg = &sol.UDG[ie*npv*nc];
    Int *bf = &mesh.bf[ie*nfe];
    //Int *perm = &mesh.perm[0];
    Int *bcm = &app.bcm[0];    
    
    sz[0] = nfe*ngf;
    for (j=0; j<ncu; j++)
        for (is=0; is<nfe; is++)
            for (k=0; k<ngf; k++)
                fhjac[j*sz[0]+is*ngf+k] = fh[j*sz[0]+is*ngf+k]*jacf[is*ngf+k];

    nb = nfe*nch;
    DGEMM(&chn, &chn, &npf, &nb, &ngf, &one, shapfg, &npf, fhjac, &ngf,
          &zero, Rutmp, &npf);

//     if (app.my_rank==0) {
//         print2darray(jacf,ngf,nfe);        
//         print2darray(fhjac,ngf*nfe,nch);
//         print2darray(Rutmp,npf*nfe,nch);
//     }
        
    // Note: No differences observed between both implementations (tested with g++ -O3). In any case, the time spent in this loop is negligible for any p
    for (j=0; j<ncu; j++)
        for (is=0; is<nfe; is++)
            for (k=0; k<npf; k++)
                Ru[j*npv+master.perm[is][k]] = Ru[j*npv+master.perm[is][k]] - Rutmp[j*ndf+is*npf+k];

    // Note: The first implementation is faster (tested with g++ -O) [No need to "help" the compiler since perm does not appear in the loop]
    t = clock();
    sz[0] = nfe*ngf;
    sz[1] = sz[0]*nch;
    for (m=0; m<nc; m++)
        for (j=0; j<ncu; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<ngf; k++)
                    fhujac[m*sz[1]+j*sz[0]+is*ngf+k] = fh_u[m*sz[1]+j*sz[0]+is*ngf+k]*jacf[is*ngf+k];
    elemtimes[13] += clock() - t;

    t = clock();
    na = npf*npf;
    nb = nfe*nch*ncu;
    DGEMM(&chn, &chn, &na, &nb, &ngf, &one, shapfgdotshapfc, &na, fhujac, &ngf,
          &zero, BDtmp, &na);
    elemtimes[14] += clock() - t;

    t = clock();
    sz[0] = nfe*ngf;
    sz[1] = sz[0]*nch;
    sz[2] = nfe*ngf;
    sz[3] = nch*sz[2];
    for (m=0; m<nch; m++)
        for (j=0; j<nch; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<ngf; k++)
                    fhuhjac[m*sz[1]+j*sz[0]+is*ngf+k] = fh_uh[m*sz[3]+j*sz[2]+is*ngf+k]*jacf[is*ngf+k];
    elemtimes[15] += clock() - t;

    t = clock();
    na = npf*npf;
    nb = nfe*nch*nch;
    DGEMM(&chn, &chn, &na, &nb, &ngf, &one, shapfgdotshapfc, &na, fhuhjac, &ngf,
          &zero, Ftmp, &na);
    elemtimes[16] += clock() - t;

    // Note: Not differences observed between both implementations. Tested with g++ -O3.
    t = clock();
    sz[0] = npv*ncu;
    sz[1] = nch*sz[0];
    sz[2] = npf*sz[1];
    sz[3] = npf*npf;
    sz[4] = ndf*npf;
    sz[5] = nch*ndf*npf;
    for (m=0; m<nch; m++)
        for (j=0; j<ncu; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<npf; k++)
                    for (ks=0; ks<npf; ks++)
                        elem.F[is*sz[2]+k*sz[1]+m*sz[0]+j*npv+master.perm[is][ks]] +=
                                Ftmp[m*sz[5]+j*sz[4]+is*sz[3]+k*npf+ks];
    elemtimes[17] += clock() - t;

    t = clock();
    t2 = 0;
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*npv;
    for (m=0; m<ncu; m++)
        for (j=0; j<ncu; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<npf; k++)
                    for (ks=0; ks<npf; ks++) {
                        t1 = m*sz[1]+master.perm[is][k]*sz[0] + j*npv+master.perm[is][ks];
                        elem.BD[t1] += BDtmp[t2];
                        t2++;
                    }
    elemtimes[18] += clock() - t;

    if (app.flag_q==1) {
        double *M = &elem.M[0];
        double *MiC = &elem.C[0];
        double *MiE = &elem.E[0];

        double *MiCf = &temp.MiCf[0];
        double *shapMiCf = &temp.shapMiCf[0];
        double *fhushapMiCf = &temp.fhushapMiCf[0];
        double *Df = &temp.Df[0];
        double *MiEf = &temp.MiEf[0];
        double *shapMiEf = &temp.shapMiEf[0];
        double *fhushapMiEf = &temp.fhushapMiEf[0];
        double *Ff = &temp.Ff[0];

        t = clock();
        sz[0] = npv*ndf;
        sz[1] = npv*npv;
        for (i=0; i<nd; i++)
            for (j=0; j<npv; j++)
                for (is=0; is<nfe; is++)
                    for (k=0; k<npf; k++)
                        MiCf[i*sz[0]+j*ndf+is*npf+k] = MiC[i*sz[1]+j*npv+master.perm[is][k]];

        sz[0] = ndf*ndf;
        sz[1] = ndf*npv;
        for (i=0; i<nd; i++)
            for (j=0; j<ndf; j++)
                for (is=0; is<nfe; is++)
                    for (k=0; k<npf; k++)
                        MiEf[i*sz[0]+j*ndf+is*npf+k] = MiE[i*sz[1]+j*npv+master.perm[is][k]];
        elemtimes[19] += clock() - t;

        t = clock();
        sz[4] = nfe*npv*nd;
        DGEMM(&chn, &chn, &ngf, &sz[4], &npf, &one, shapft, &ngf, &MiCf[0],
              &npf, &zero, &shapMiCf[0], &ngf);

        sz[4] = nfe*ndf*nd;
        DGEMM(&chn, &chn, &ngf, &sz[4], &npf, &one, shapft, &ngf, &MiEf[0],
              &npf, &zero, &shapMiEf[0], &ngf);
        elemtimes[20] += clock() - t;

        t = clock();
        sz[0] = nfe*ngf;
        sz[1] = sz[0]*ncu;
        sz[2] = sz[1]*npv;
        sz[3] = ncu*ncu*nfe*ngf;
        sz[4] = npv*nfe*ngf;
        sz[5] = ndf*nfe*ngf;
        if (nd==1) {
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[n*sz[0]+m*ngf+i];

            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[n*sz[0]+m*ngf+i];
        }
        if (nd==2) {
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[n*sz[0]+m*ngf+i]+
                                        fhujac[2*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[sz[4]+n*sz[0]+m*ngf+i];

            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[n*sz[0]+m*ngf+i]+
                                        fhujac[2*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[sz[5]+n*sz[0]+m*ngf+i];
        }
        if (nd==3) {
            for (k=0; k<ncu; k++)
                for (n=0; n<npv; n++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[n*sz[0]+m*ngf+i]+
                                        fhujac[2*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[sz[4]+n*sz[0]+m*ngf+i]+
                                        fhujac[3*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiCf[2*sz[4]+n*sz[0]+m*ngf+i];

            for (n=0; n<ndf; n++)
                for (k=0; k<ncu; k++)
                    for (j=0; j<ncu; j++)
                        for (m=0; m<nfe; m++)
                            for (i=0; i<ngf; i++)
                                fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i] =
                                        fhujac[sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[n*sz[0]+m*ngf+i]+
                                        fhujac[2*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[sz[5]+n*sz[0]+m*ngf+i]+
                                        fhujac[3*sz[3]+k*sz[1]+j*sz[0]+m*ngf+i]*shapMiEf[2*sz[5]+n*sz[0]+m*ngf+i];
        }
        elemtimes[21] += clock() - t;

        t = clock();
        sz[3] = nfe*ncu*npv*ncu;
        DGEMM(&chn, &chn, &npf, &sz[3], &ngf, &one, shapfg, &npf, &fhushapMiCf[0],
              &ngf, &zero, &Df[0], &npf);
        elemtimes[22] += clock() - t;

        t = clock();
        sz[0] = ncu*npv;
        sz[1] = npv*ncu*npv;
        sz[2] = ncu*ndf;
        sz[3] = npv*ncu*ndf;
        for (k=0; k<ncu; k++)
            for (n=0; n<npv; n++)
                for (j=0; j<ncu; j++)
                    for (m=0; m<nfe; m++)
                        for (i=0; i<npf; i++)
                            elem.BD[k*sz[1]+n*sz[0]+j*npv+master.perm[m][i]] += Df[k*sz[3]+n*sz[2]+j*ndf+m*npf+i];
        elemtimes[23] += clock() - t;

        t = clock();
        sz[3] = nfe*ncu*ncu*ndf;
        DGEMM(&chn, &chn, &npf, &sz[3], &ngf, &one, shapfg, &npf, &fhushapMiEf[0],
              &ngf, &zero, &Ff[0], &npf);
        elemtimes[24] += clock() - t;

        t = clock();
        sz[0] = ncu*npv;
        sz[1] = ncu*ncu*npv;
        sz[2] = ncu*ndf;
        sz[3] = ncu*ncu*ndf;
        for (j=0; j<ndf; j++)
            for (m=0; m<ncu; m++)
                for (i=0; i<ncu; i++)
                    for (is=0; is<nfe; is++)
                        for (k=0; k<npf; k++)
                            elem.F[j*sz[1]+m*sz[0]+i*npv+master.perm[is][k]] -= Ff[j*sz[3]+m*sz[2]+i*ndf+is*npf+k];
        elemtimes[25] += clock() - t;
    }

    t = clock();
    Int ib, bn;
    na = npf*npf;
    for (is=0; is<nfe; is++)
        if (bf[is] < 0) {
            sz[0] = ngf*nfe;
            for (j=0; j<ncd; j++)
                for (k=0; k<ngf; k++)
                    pft[j*ngf+k] = pgf[j*sz[0]+is*ngf+k];

            for (j=0; j<nc; j++)
                for (k=0; k<ngf; k++)
                    uft[j*ngf+k] = ugf[j*sz[0]+is*ngf+k];
            
            if (nco>0) {
                for (j=0; j<nco; j++)
                    for (k=0; k<ngf; k++)
                        oft[j*ngf+k] = ogf[j*sz[0]+is*ngf+k];
            }                

            for (j=0; j<nch; j++)
                for (k=0; k<ngf; k++)
                    uht[j*ngf+k] = uhg[j*sz[0]+is*ngf+k];

            for (j=0; j<nd; j++)
                for (k=0; k<ngf; k++)
                    nlt[j*ngf+k] = nlg[j*sz[0]+is*ngf+k];

            ib = -bf[is]-1;
            bn = bcm[ib];            
            bouflux(mesh, master, app, sol, temp, ib, bn, is, ie, 1);            

            sz[0] = ngf*nfe;
            sz[1] = nch*sz[0];
            for (j=0; j<nch; j++)
                for (k=0; k<ngf; k++)
                    fhjac[j*sz[0]+is*ngf+k] = fb[j*ngf+k]*jacf[is*ngf+k];

            for (m=0; m<nc; m++)
                for (j=0; j<nch; j++)
                    for (k=0; k<ngf; k++)
                        fhujac[m*sz[1]+j*sz[0]+is*ngf+k] = fb_u[m*nch*ngf+j*ngf+k]*jacf[is*ngf+k];

            for (m=0; m<nch; m++)
                for (j=0; j<nch; j++)
                    for (k=0; k<ngf; k++)
                        fhuhjac[m*sz[1]+j*sz[0]+is*ngf+k] = fb_uh[m*nch*ngf+j*ngf+k]*jacf[is*ngf+k];

            // Note: No significant difference between for + DGEMV vs. DGEMM taking advantage of LD in BLAS
//            for (i=0; i<nch; i++) {
//                ifhjac = i*sz[0]+is*ngf;
//                iRutmp = i*ndf+is*npf;
//                DGEMV(&chn, &npf, &ngf, &one, shapfg, &npf, &fhjac[ifhjac], &inc,
//                      &zero, &Rutmp[iRutmp], &inc);
//            }
            sz[2] = ngf*nfe;
            sz[3] = nch;
            sz[4] = npf*nfe;
            DGEMM(&chn, &chn, &npf, &sz[3], &ngf, &one, shapfg, &npf, &fhjac[is*ngf],
                  &sz[2], &zero, &Rutmp[is*npf], &sz[4]);

//            for (j=0; j<nc; j++)
//                for (i=0; i<nch; i++) {
//                    ifhujac = j*sz[1] + i*sz[0] + is*ngf;
//                    iBDtmp = j*npf*npf*nfe*nch + i*npf*npf*nfe + is*npf*npf;
//                    DGEMV(&chn, &na, &ngf, &one, shapfgdotshapfc, &na, &fhujac[ifhujac], &inc,
//                          &zero, &BDtmp[iBDtmp], &inc);
//                }
            sz[2] = ngf*nfe;
            sz[3] = nch*nc;
            sz[4] = npf*npf*nfe;
            DGEMM(&chn, &chn, &na, &sz[3], &ngf, &one, shapfgdotshapfc, &na, &fhujac[is*ngf],
                  &sz[2], &zero, &BDtmp[is*npf*npf], &sz[4]);

//            for (j=0; j<nch; j++)
//                for (i=0; i<nch; i++) {
//                    ifhuhjac = j*sz[1] + i*sz[0] + is*ngf;
//                    iFtmp = j*npf*npf*nfe*nch + i*npf*npf*nfe + is*npf*npf;
//                    DGEMV(&chn, &na, &ngf, &one, shapfgdotshapfc, &na, &fhuhjac[ifhuhjac], &inc,
//                          &zero, &Ftmp[iFtmp], &inc);
//                }
            sz[2] = ngf*nfe;
            sz[3] = nch*nch;
            sz[4] = npf*npf*nfe;
            DGEMM(&chn, &chn, &na, &sz[3], &ngf, &one, shapfgdotshapfc, &na, &fhuhjac[is*ngf],
                  &sz[2], &zero, &Ftmp[is*npf*npf], &sz[4]);

            if (app.flag_q==1) {
                double *shapMiCf = &temp.shapMiCf[0];
                double *fhushapMiCf = &temp.fhushapMiCf[0];
                double *Df = &temp.Df[0];
                double *shapMiEf = &temp.shapMiEf[0];
                double *fhushapMiEf = &temp.fhushapMiEf[0];
                double *Ff = &temp.Ff[0];

                sz[0] = nfe*ngf;
                sz[1] = sz[0]*ncu;
                sz[2] = sz[1]*npv;
                sz[3] = ncu*ncu*nfe*ngf;
                sz[4] = npv*nfe*ngf;
                sz[5] = ndf*nfe*ngf;
                if (nd==1) {
                    for (k=0; k<ncu; k++)
                        for (n=0; n<npv; n++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[n*sz[0]+is*ngf+i];

                    for (n=0; n<ndf; n++)
                        for (k=0; k<ncu; k++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[n*sz[0]+is*ngf+i];
                }
                if (nd==2) {
                    for (k=0; k<ncu; k++)
                        for (n=0; n<npv; n++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[n*sz[0]+is*ngf+i]+
                                            fhujac[2*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[sz[4]+n*sz[0]+is*ngf+i];

                    for (n=0; n<ndf; n++)
                        for (k=0; k<ncu; k++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[n*sz[0]+is*ngf+i]+
                                            fhujac[2*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[sz[5]+n*sz[0]+is*ngf+i];
                }
                if (nd==3) {
                    for (k=0; k<ncu; k++)
                        for (n=0; n<npv; n++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiCf[k*sz[2]+n*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[n*sz[0]+is*ngf+i]+
                                            fhujac[2*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[sz[4]+n*sz[0]+is*ngf+i]+
                                            fhujac[3*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiCf[2*sz[4]+n*sz[0]+is*ngf+i];

                    for (n=0; n<ndf; n++)
                        for (k=0; k<ncu; k++)
                            for (j=0; j<ncu; j++)
                                for (i=0; i<ngf; i++)
                                    fhushapMiEf[n*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i] =
                                            fhujac[sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[n*sz[0]+is*ngf+i]+
                                            fhujac[2*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[sz[5]+n*sz[0]+is*ngf+i]+
                                            fhujac[3*sz[3]+k*sz[1]+j*sz[0]+is*ngf+i]*shapMiEf[2*sz[5]+n*sz[0]+is*ngf+i];
                }

                sz[1] = npf*nfe;
                sz[2] = ngf*nfe;
                sz[3] = ncu*npv*ncu;
                DGEMM(&chn, &chn, &npf, &sz[3], &ngf, &one, shapfg, &npf, &fhushapMiCf[is*ngf],
                      &sz[2], &zero, &Df[is*npf], &sz[1]);

                sz[1] = npf*nfe;
                sz[2] = ngf*nfe;
                sz[3] = ncu*ncu*ndf;
                DGEMM(&chn, &chn, &npf, &sz[3], &ngf, &one, shapfg, &npf, &fhushapMiEf[is*ngf],
                      &sz[2], &zero, &Ff[is*npf], &sz[1]);

                // shapfg: npf // ngf
                // fhushapMiCf: ngf // nfe / ncu / npv / ncu
                // fhushapMiEf: ngf // nfe / ncu / ncu / ndf
                // Df: npf // nfe / ncu / npv / ncu
                // Ff: npf // nfe / ncu / ncu / ndf
            }
        }
    elemtimes[26] += clock() - t;

    for (j=0; j<nch; j++)
        for (is=0; is<nfe; is++)
            for (k=0; k<npf; k++)
                Rh[is*npf*nch+k*nch+j] = -Rutmp[j*ndf+is*npf+k];

    t = clock();
    sz[2] = nch*ndf*npv*ncu;    
    for (i=0; i<sz[2]; i++)
        GK[i] = 0.0;
    //memset(&GK[0], 0.0, sizeof(double)*sz[2]);
    sz[2] = npf*nch;
    sz[3] = ndf*nch;
    sz[4] = npv*ndf*nch;
    sz[5] = npf*npf;
    sz[6] = ndf*npf;
    sz[7] = nch*ndf*npf;
    t2 = 0;
    for (m=0; m<ncu; m++)
        for (j=0; j<nch; j++)
            for (is=0; is<nfe; is++) {
                t00 = m*sz[4] + is*sz[2] + j;
                for (k=0; k<npf; k++) {
                    t0 = t00 + master.perm[is][k]*sz[3];
                    for (ks=0; ks<npf; ks++) {
                        t1 = t0 + ks*nch;
                        GK[t1] = BDtmp[t2];
                        t2++;
                    }
                }
            }
    elemtimes[27] += clock() - t;

    t = clock();
    sz[0] = nch*ndf*nch*ndf;
    for (i=0; i<sz[0]; i++)
        H[i] = 0.0;
    //memset(&H[0], 0.0, sizeof(double)*sz[0]);
    sz[0] = ndf*nch;
    sz[1] = nch*sz[0];
    sz[2] = npf*npf;
    sz[3] = ndf*npf;
    sz[4] = nch*sz[3];
    for (m=0; m<nch; m++)
        for (j=0; j<nch; j++)
            for (is=0; is<nfe; is++)
                for (k=0; k<npf; k++)
                    for (ks=0; ks<npf; ks++)
                        H[(is*npf+k)*sz[1]+m*sz[0]+(is*npf+ks)*nch+j] =
                                Ftmp[m*sz[4]+j*sz[3]+is*sz[2]+k*npf+ks];
    elemtimes[28] += clock() - t;

    if (app.flag_q==1) {
        double *Df = &temp.Df[0];
        double *Ff = &temp.Ff[0];

        t = clock();
        sz[0] = ndf*ncu;
        sz[1] = npv*ndf*ncu;
        for (m=0; m<ncu; m++)
            for (j=0; j<npv; j++)
                for (i=0; i<ncu; i++)
                    for (k=0; k<ndf; k++)
                        elem.GK[m*sz[1]+j*sz[0]+k*ncu+i] += Df[m*sz[1]+j*sz[0]+i*ndf+k];

        sz[0] = ndf*ncu;
        sz[1] = ncu*ndf*ncu;
        for (j=0; j<ndf; j++)
            for (m=0; m<ncu; m++)
                for (i=0; i<ncu; i++)
                    for (k=0; k<ndf; k++)
                        elem.H[j*sz[1]+m*sz[0]+k*ncu+i] -= Ff[j*sz[1]+m*sz[0]+i*ndf+k];
        elemtimes[31] += clock() - t;
    }
}

void qint(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app, tempstruct &temp, Int ie)
{
    Int inc = 1, i, j, is, m, n, k, t2, iA, iB, iC, info;
    char chn = 'N';
    double one = 1.0, zero = 0.0, minusone = -1.0;

    /* Get dimensions */
    Int npv, ngv, ncd, nc, ncu, nco, ncq, nch, npf, nfe, ngf, nd, ndf, nd1, sza, szb;
    nd = app.nd;
    ncd = app.ncd;
    nfe = master.nfe;
    npv = master.npv;
    npf = master.npf[0];
    ngf = master.ngf[0];
    ngv = master.ngv;
    nc  = app.nc;
    ncu = app.ncu;
    nch = app.nch;
    nco = app.nco;
    ndf = npf*nfe;
    nd1 = nd+1;
    sza = ngf*nfe;
    szb = npf*npf;
            
    elementgeom(mesh, master, app, temp, ie);
    facegeom(mesh, master, app, temp, ie);

    double *M  = &elem.M[0];
    double *C  = &elem.C[0];
    double *E  = &elem.E[0];
    double fc_q = app.fc_q;

    double *shapvgdotshapvl, *shapfgdotshapfc;
    shapvgdotshapvl = &master.shapvgdotshapvl[0];
    shapfgdotshapfc = &master.shapfgdotshapfc[0][0];

    double *jacv, *Xx, *jacf, *nlg, *nlgjac, *Etmp;
    jacv = &temp.jacg[0];
    Xx = &temp.Xxg[0];
    jacf = &temp.jacgf[0];
    nlg = &temp.nlgf[0];
    Etmp = &temp.Etmp[0];
    nlgjac = &temp.nlgjac[0];
    //Int *perm = &mesh.perm[0];

    Int sz[10];
    sz[1] = ngv*nd;
    sz[2] = npv*npv;
    sz[6] = npv*npv*ngv;

    /* Compute the Mass matrix */
    DGEMV(&chn, &sz[2], &ngv, &fc_q, &shapvgdotshapvl[0], &sz[2], jacv,
                    &inc, &zero, &M[0], &inc);

    /* Compute the Convection matrix. TODO: Compare performance with implementation without loop */
    for (i=0; i<nd; i++) {
        iA = sz[6];
        iB = i*ngv;
        iC = i*sz[2];

        DGEMV(&chn, &sz[2], &ngv, &one, &shapvgdotshapvl[iA], &sz[2], &Xx[iB],
              &inc, &zero, &C[iC], &inc);

        for (j=1; j<nd; j++) {
            iA = (j+1)*sz[6];
            iB = j*sz[1] + i*ngv;
            DGEMV(&chn, &sz[2], &ngv, &one, &shapvgdotshapvl[iA], &sz[2], &Xx[iB],
              &inc, &one, &C[iC], &inc);
        }
    }

    for (i=0; i<npv*ndf*nd; i++)
        E[i] = 0.0;

    for (i=0; i<nd; i++) {
        for (j=0; j<sza; j++)
            nlgjac[j] = nlg[i*sza+j]*jacf[j];

        DGEMM(&chn, &chn, &szb, &nfe, &ngf, &one, shapfgdotshapfc, &szb, nlgjac,
               &ngf, &zero, Etmp, &szb);

        iC = i*npv*ndf;
        for (is=0; is<nfe; is++)
            for (j=0; j<npf; j++)
                for (k=0; k<npf; k++) {
                    n = (is*npf+j)*npv + master.perm[is][k];
                    E[iC+n] = Etmp[is*szb+j*npf+k];
                }
    }
}

void getQ(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app,
          solstruct &sol, tempstruct &temp, double* u, double* q, double* UH, 
          double* sh, Int ie)
{
    /* TODO: Include the option to store MiC, MiE instead of recomputing them (e.g. through an efficiencyFlag). This is important for quasi-Newton */
    /* TODO: No source term in q equation? */
    Int inc = 1, i, j, k, n, iB, iC, info;
    char chn = 'N', charu = 'U';
    double one = 1.0, zero = 0.0, minusone = -1.0;
    
    /* Get dimensions */
    Int nd = master.nd;
    Int nfe = master.nfe;
    Int npv = master.npv;
    Int npf = master.npf[0];
    Int nc  = app.nc;
    Int ncu = app.ncu;
    Int nch = app.nch;
    Int ndf = master.ndf;
    Int ncq = app.ncq;
    
    double *M  = &elem.M[0];
    double *C  = &elem.C[0];
    double *E  = &elem.E[0];
    
    qint(elem, mesh, master, app, temp, ie);
    
//     double *u, *q, *sh, *uh;
//     u = &sol.UDG[ie*npv*nc];
//     q = &sol.UDG[ie*npv*nc+npv*ncu];
//     sh = &sol.SH[ie*npv*nc+npv*ncu];
//     uh = &temp.uh[0];    
//     for (j=0; j<ndf; j++) {
//         n = mesh.elcon[ie*ndf+j];
//         for (k=0; k<nch; k++)
//             uh[k*ndf+j] = sol.UH[n*nch+k];
//     }
    
    double *uh = &temp.uh[0];
    for (j=0; j<ndf; j++) {
        n = mesh.elcon[ie*ndf+j];
        for (k=0; k<nch; k++)
            uh[k*ndf+j] = UH[n*nch+k];
    }        
    
    DGEMM(&chn, &chn, &npv, &ncq, &npv, &one, M, &npv, sh, &npv,
              &zero, q, &npv);
    for (i=0; i<nd; i++) {
        iB = i*npv*ncu;
        iC = i*npv*npv;
        DGEMM(&chn, &chn, &npv, &ncu, &npv, &one, &C[iC], &npv, u, &npv,
              &one, &q[iB], &npv);
        iC = i*npv*ndf;
        DGEMM(&chn, &chn, &npv, &nch, &ndf, &minusone, &E[iC], &npv, uh, &ndf,
              &one, &q[iB], &npv);
    }
    
    /* Factor M using Cholesky decomposition */
    DPOTRF(&charu,&npv,M,&npv,&info);

    /* Get q */
    n = ncu*nd;
    DPOTRS(&charu,&npv,&n,M,&npv,q,&npv,&info);
    
    /* TODO: update q for wave problems ?? */
}

void assembleElementMatrixVector(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app,
                solstruct &sol, tempstruct &temp, Int ie, double* elemtimes)
{
    clock_t t;
    int computeJacobian = 1;
    
    t = clock();
    if (app.flag_q==1) 
        qint(elem, mesh, master, app, temp, ie);     // About 10 times faster than elementint for 2D triangles (this result applies for different p's)      
    elemtimes[0] += clock() - t;
    
    t = clock();
    elementint(elem, mesh, master, app, sol, temp, &elemtimes[0], ie, computeJacobian);
    elemtimes[1] += clock() - t;        
    
    t = clock();
    faceint(elem, mesh, master, app, sol, temp, &elemtimes[0], ie, computeJacobian);
    elemtimes[2] += clock() - t;        
}

void assembleElementVector(elemstruct &elem, meshstruct &mesh, masterstruct &master, appstruct &app,
                solstruct &sol, tempstruct &temp, Int ie)
{
    int computeJacobian = 0;
    double elemtimes[32];
    
    if (app.flag_q==1) 
        qint(elem, mesh, master, app, temp, ie);     // About 10 times faster than elementint for 2D triangles (this result applies for different p's)      
    
    elementint(elem, mesh, master, app, sol, temp, &elemtimes[0], ie, computeJacobian);
    faceint(elem, mesh, master, app, sol, temp, &elemtimes[0], ie, computeJacobian);
}

#endif
