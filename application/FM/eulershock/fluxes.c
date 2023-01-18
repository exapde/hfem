
void source_euler(double *sr, double *sr_u, double *pg, double *udg, 
        double *param, double time, int ng, int ncu, int nc, int nd, int ne)
{
    
    int i, j, k, l, m, n;
    
    for (j=0; j<ne; j++) {     // loop over each element
        k = j*ng*nc;           // pointer to element j
        l = j*ng*nc*nc;
        for (i=0; i<ng; i++) { // loop over each gauss points       
            m   = i + k;
            n   = i + l;
            
            sr[0*ng+m]  = 0.0;
            sr[1*ng+m]  = 0.0;    
            sr[2*ng+m]  = 0.0;
            sr[3*ng+m]  = 0.0;           
            
            sr_u[0*ng+n] = 0.0;
            sr_u[1*ng+n] = 0.0;
            sr_u[2*ng+n] = 0.0;
            sr_u[3*ng+n] = 0.0;        
            sr_u[4*ng+n] = 0.0;
            sr_u[5*ng+n] = 0.0;
            sr_u[6*ng+n] = 0.0;
            sr_u[7*ng+n] = 0.0;        
            sr_u[8*ng+n] = 0.0;
            sr_u[9*ng+n] = 0.0;
            sr_u[10*ng+n] = 0.0;
            sr_u[11*ng+n] = 0.0;        
            sr_u[12*ng+n] = 0.0;
            sr_u[13*ng+n] = 0.0;
            sr_u[14*ng+n] = 0.0;
            sr_u[15*ng+n] = 0.0;
                
        }
    }
}

void flux_euler(double *fx, double *fy, double *fx_u, double *fy_u, double *pg, 
        double *udg, double *param, double time, int ng, int ncu, int nc, int nd, int ne)
{
    
    double gam, gam1, c23, c12;
    double r, ru, rv, rE;    
    double r1, u, v, E, uu, vv, uv, q, p, h;    
    int    i, j, k, l, m, n;
            
    // compute parameters
    gam  = param[0];  
    //gam  = 1.4;
    //Re   = param[2];
    //Pr   = param[3];
    //Minf = param[4];
    gam1 = gam - 1.0;
    //Re1  = 1/Re;
    //M2   = Minf*Minf; 
    c23  = 2.0/3.0;
    c12  = 1.0/2.0;
    //fc   = 1.0/(gam1*M2*Re*Pr);
    
    for (j=0; j<ne; j++) {     // loop over each element
        k = j*ng*nc;           // pointer to element j
        l = j*ng*nc*nc;
        for (i=0; i<ng; i++) { // loop over each gauss points       
            m   = i + k;
            n   = i + l;
            r   = udg[0*ng+m];
            ru  = udg[1*ng+m];
            rv  = udg[2*ng+m];
            rE  = udg[3*ng+m];

            r1   = 1.0/r;
            u    = ru*r1;
            v    = rv*r1;
            E    = rE*r1;
            uu   = u*u;
            vv   = v*v;
            uv   = u*v;
            q    = 0.5*(uu+vv);
            p    = gam1*(rE-r*q);
            h    = E+p*r1;

            fx[0*ng+m]  = ru;
            fx[1*ng+m]  = ru*u + p;    
            fx[2*ng+m]  = rv*u;
            fx[3*ng+m]  = ru*h;
            
            fy[0*ng+m]  = rv;
            fy[1*ng+m]  = ru*v;
            fy[2*ng+m]  = rv*v + p;    
            fy[3*ng+m]  = rv*h;          

            fx_u[0*ng+n] = 0.0;
            fx_u[1*ng+n] = 0.5*((gam-3)*uu+gam1*vv);
            fx_u[2*ng+n] = -uv;
            fx_u[3*ng+n] = 2*gam1*u*q-gam*E*u;        
            fx_u[4*ng+n] = 1.0;
            fx_u[5*ng+n] = (3-gam)*u;
            fx_u[6*ng+n] = v;
            fx_u[7*ng+n] = gam*E-0.5*gam1*(3*uu+vv);        
            fx_u[8*ng+n] = 0.0;
            fx_u[9*ng+n] = -gam1*v;
            fx_u[10*ng+n] = u;
            fx_u[11*ng+n] = -gam1*u*v;        
            fx_u[12*ng+n] = 0.0;
            fx_u[13*ng+n] = gam1;
            fx_u[14*ng+n] = 0.0;
            fx_u[15*ng+n] = gam*u;
    
            fy_u[0*ng+n] = 0.0;
            fy_u[1*ng+n] = -uv;
            fy_u[2*ng+n] =  0.5*((gam-3)*vv+gam1*uu);
            fy_u[3*ng+n] = 2*gam1*v*q-gam*E*v;        
            fy_u[4*ng+n] = 0.0;
            fy_u[5*ng+n] = v;
            fy_u[6*ng+n] = -gam1*u;
            fy_u[7*ng+n] = -gam1*uv;        
            fy_u[8*ng+n] = 1.0;
            fy_u[9*ng+n] = u;
            fy_u[10*ng+n] = (3-gam)*v;
            fy_u[11*ng+n] = gam*E-0.5*gam1*(3*vv+uu);        
            fy_u[12*ng+n] = 0.0;
            fy_u[13*ng+n] = 0.0;
            fy_u[14*ng+n] = gam1;
            fy_u[15*ng+n] = gam*v;
            
        }
    }
}


void fhat_euler(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *uhg, double *nlg, 
          double *param, double time, int ng, int nch, int ncu, int nc, int nd, int ne)
{
    
    double gam, epslm, gam1, signun;
    double r, ru, rv, rE, ur, uru, urv, urE;    
    double r1, u, v, E, uu, vv, uv, q, p, h;     
    double nx, ny;
    double fx[4], fy[4], fx_uh[16], fy_uh[16];
    double um1, um2, um3, um4, vm1, vm2, vm3, vm4;
    double r1m1, r1m2, r1m3, r1m4, qm1, qm2, qm3, qm4;
    double pm1, pm2, pm3, pm4, c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4, un, unm1, unm2, unm3, unm4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;            
    
    int    i, j, k, l, m, n, d, e, dim;
    
    dim  = 2;
    gam  = param[0];
    epslm= param[1];       
    gam1 = gam - 1.0;
    
    for (j=0; j<ne; j++) {     // loop over each element
        k = j*ng*nc;           // pointer to element j
        d = j*ng*dim;
        l = j*ng*nc*nc;                      
        for (i=0; i<ng; i++) { // loop over each gauss point    
            m   = i + k;
            e   = i + d;
            n   = i + l;
            
            r   = uhg[0*ng+m];
            ru  = uhg[1*ng+m];
            rv  = uhg[2*ng+m];
            rE  = uhg[3*ng+m];
            
            ur  = udg[0*ng+m];
            uru = udg[1*ng+m];
            urv = udg[2*ng+m];
            urE = udg[3*ng+m];
            
            nx  = nlg[0*ng+e];              
            ny  = nlg[1*ng+e];   
            
            r1   = 1.0/r;
            u    = ru*r1;
            v    = rv*r1;
            E    = rE*r1;
            uu   = u*u;
            vv   = v*v;
            uv   = u*v;
            q    = 0.5*(uu+vv);
            p    = gam1*(rE-r*q);
            h    = E+p*r1;                                                                                    
        
            fx[0]  = ru;
            fx[1]  = ru*u + p;    
            fx[2]  = rv*u;
            fx[3]  = ru*h;                        
            
            fy[0]  = rv;
            fy[1]  = ru*v;
            fy[2]  = rv*v + p;    
            fy[3]  = rv*h;          

            fx_uh[0] = 0.0;
            fx_uh[1] = 0.5*((gam-3)*uu+gam1*vv);
            fx_uh[2] = -uv;
            fx_uh[3] = 2*gam1*u*q-gam*E*u;        
            fx_uh[4] = 1.0;
            fx_uh[5] = (3-gam)*u;
            fx_uh[6] = v;
            fx_uh[7] = gam*E-0.5*gam1*(3*uu+vv);        
            fx_uh[8] = 0.0;
            fx_uh[9] = -gam1*v;
            fx_uh[10] = u;
            fx_uh[11] = -gam1*u*v;        
            fx_uh[12] = 0.0;
            fx_uh[13] = gam1;
            fx_uh[14] = 0.0;
            fx_uh[15] = gam*u;

            fy_uh[0] = 0.0;
            fy_uh[1] = -uv;
            fy_uh[2] =  0.5*((gam-3)*vv+gam1*uu);
            fy_uh[3] = 2*gam1*v*q-gam*E*v;        
            fy_uh[4] = 0.0;
            fy_uh[5] = v;
            fy_uh[6] = -gam1*u;
            fy_uh[7] = -gam1*uv;        
            fy_uh[8] = 1.0;
            fy_uh[9] = u;
            fy_uh[10] = (3-gam)*v;
            fy_uh[11] = gam*E-0.5*gam1*(3*vv+uu);        
            fy_uh[12] = 0.0;
            fy_uh[13] = 0.0;
            fy_uh[14] = gam1;
            fy_uh[15] = gam*v;
            
                        
            r1m1 = -1.0/(r*r);
            r1m2 = 0;
            r1m3 = 0;
            r1m4 = 0;
                        
            um1 = ru*r1m1;
            um2 = r1;
            um3 = 0.0;
            um4 = 0.0;
            
            vm1 = rv*r1m1;
            vm2 = 0.0;
            vm3 = r1;
            vm4 = 0.0;
                        
            qm1 = u*um1 + v*vm1;
            qm2 = u*um2 + v*vm2;
            qm3 = u*um3 + v*vm3;
            qm4 = u*um4 + v*vm4;
            
            pm1  = gam1*(   -  q  - r*qm1);
            pm2  = gam1*(         - r*qm2);
            pm3  = gam1*(         - r*qm3);
            pm4  = gam1*(1.0      - r*qm4);

            c2   = gam*p*r1;
            c2m1 = gam*(pm1*r1 + p*r1m1);
            c2m2 = gam*(pm2*r1 + p*r1m2);
            c2m3 = gam*(pm3*r1 + p*r1m3);
            c2m4 = gam*(pm4*r1 + p*r1m4);

            c    = sqrt(c2);
            cm1  = 0.5*c2m1/c;
            cm2  = 0.5*c2m2/c;
            cm3  = 0.5*c2m3/c;
            cm4  = 0.5*c2m4/c;

            un   = u*nx   + v*ny;
            unm1 = um1*nx + vm1*ny;
            unm2 = um2*nx + vm2*ny;
            unm3 = um3*nx + vm3*ny;
            unm4 = um4*nx + vm4*ny;

            rlam   = fabs(un)+c;
            signun = (un < 0) ? -1 : 1;   
            if (un==0) signun = 0;
            rlamm1 = signun*unm1+cm1;
            rlamm2 = signun*unm2+cm2;
            rlamm3 = signun*unm3+cm3;
            rlamm4 = signun*unm4+cm4;

            if (epslm>0) {
                rlam = 0.5*(rlam*rlam/(epslm*c)+epslm*c);
                rlamm1 = rlam*rlamm1/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm1);
                rlamm2 = rlam*rlamm2/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm2);
                rlamm3 = rlam*rlamm3/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm3);
                rlamm4 = rlam*rlamm4/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm4);
            }   
                                    
            fh[0*ng+m]  = fx[0]*nx + fy[0]*ny + rlam*(ur - r);
            fh[1*ng+m]  = fx[1]*nx + fy[1]*ny + rlam*(uru - ru);    
            fh[2*ng+m]  = fx[2]*nx + fy[2]*ny + rlam*(urv - rv);
            fh[3*ng+m]  = fx[3]*nx + fy[3]*ny + rlam*(urE - rE);            
            
            fh_u[0*ng+n] = rlam;
            fh_u[1*ng+n] = 0.0;
            fh_u[2*ng+n] = 0.0;
            fh_u[3*ng+n] = 0.0;        
            fh_u[4*ng+n] = 0.0;
            fh_u[5*ng+n] = rlam;
            fh_u[6*ng+n] = 0.0;
            fh_u[7*ng+n] = 0.0;        
            fh_u[8*ng+n] = 0.0;
            fh_u[9*ng+n] = 0.0;
            fh_u[10*ng+n] = rlam;
            fh_u[11*ng+n] = 0.0;        
            fh_u[12*ng+n] = 0.0;
            fh_u[13*ng+n] = 0.0;
            fh_u[14*ng+n] = 0.0;
            fh_u[15*ng+n] = rlam;
                        
            fh_uh[0*ng+n] = fx_uh[0]*nx + fy_uh[0]*ny + rlamm1*(ur - r) - rlam;
            fh_uh[1*ng+n] = fx_uh[1]*nx + fy_uh[1]*ny + rlamm1*(uru - ru);
            fh_uh[2*ng+n] = fx_uh[2]*nx + fy_uh[2]*ny + rlamm1*(urv - rv);
            fh_uh[3*ng+n] = fx_uh[3]*nx + fy_uh[3]*ny + rlamm1*(urE - rE);        
            fh_uh[4*ng+n] = fx_uh[4]*nx + fy_uh[4]*ny + rlamm2*(ur - r);
            fh_uh[5*ng+n] = fx_uh[5]*nx + fy_uh[5]*ny + rlamm2*(uru - ru) - rlam;
            fh_uh[6*ng+n] = fx_uh[6]*nx + fy_uh[6]*ny + rlamm2*(urv - rv);
            fh_uh[7*ng+n] = fx_uh[7]*nx + fy_uh[7]*ny + rlamm2*(urE - rE);        
            fh_uh[8*ng+n] = fx_uh[8]*nx + fy_uh[8]*ny + rlamm3*(ur - r);
            fh_uh[9*ng+n] = fx_uh[9]*nx + fy_uh[9]*ny + rlamm3*(uru - ru);
            fh_uh[10*ng+n] = fx_uh[10]*nx + fy_uh[10]*ny + rlamm3*(urv - rv) - rlam;
            fh_uh[11*ng+n] = fx_uh[11]*nx + fy_uh[11]*ny + rlamm3*(urE - rE);              
            fh_uh[12*ng+n] = fx_uh[12]*nx + fy_uh[12]*ny + rlamm4*(ur - r);
            fh_uh[13*ng+n] = fx_uh[13]*nx + fy_uh[13]*ny + rlamm4*(uru - ru);
            fh_uh[14*ng+n] = fx_uh[14]*nx + fy_uh[14]*ny + rlamm4*(urv - rv);
            fh_uh[15*ng+n] = fx_uh[15]*nx + fy_uh[15]*ny + rlamm4*(urE - rE) - rlam;
                                    
        }
    }
    
}


void getan_euler(double *An, double *Anm, double *uhg, double *param, double *nl, int absolute, int ng, int nc, int ne)
{
    
    double gam, epslm, gam1, signun;
    double r, ru, rv, rE, ur, uru, urv, urE;    
    double r1, u, v, E, uu, vv, uv, q, p, h;     
    double nx, ny;    
    double um1, um2, um3, um4, vm1, vm2, vm3, vm4;
    double r1m1, r1m2, r1m3, r1m4, qm1, qm2, qm3, qm4;
    double pm1, pm2, pm3, pm4, c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4, un, unm1, unm2, unm3, unm4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;            
    
    int    i, j, k, l, m, n, d, e, dim;
    
    dim  = 2;
    gam  = param[0];    
    epslm= param[1];
    gam1 = gam - 1.0;
    
    for (j=0; j<ne; j++) {     // loop over each element
        k = j*ng*nc;           // pointer to element j
        d = j*ng*dim;
        l = j*ng*nc*nc;                      
        for (i=0; i<ng; i++) { // loop over each gauss point    
            m   = i + k;
            e   = i + d;
            n   = i + l;
            
            r   = uhg[0*ng+m];
            ru  = uhg[1*ng+m];
            rv  = uhg[2*ng+m];
            rE  = uhg[3*ng+m];
                        
            nx  = nl[0*ng+e];              
            ny  = nl[1*ng+e];   
            
            r1   = 1.0/r;
            u    = ru*r1;
            v    = rv*r1;
            E    = rE*r1;
            uu   = u*u;
            vv   = v*v;
            uv   = u*v;
            q    = 0.5*(uu+vv);
            p    = gam1*(rE-r*q);
            h    = E+p*r1;                                                                                    
                    
                        
            r1m1 = -1.0/(r*r);
            r1m2 = 0;
            r1m3 = 0;
            r1m4 = 0;
                        
            um1 = ru*r1m1;
            um2 = r1;
            um3 = 0.0;
            um4 = 0.0;
            
            vm1 = rv*r1m1;
            vm2 = 0.0;
            vm3 = r1;
            vm4 = 0.0;
                        
            qm1 = u*um1 + v*vm1;
            qm2 = u*um2 + v*vm2;
            qm3 = u*um3 + v*vm3;
            qm4 = u*um4 + v*vm4;
            
            pm1  = gam1*(   -  q  - r*qm1);
            pm2  = gam1*(         - r*qm2);
            pm3  = gam1*(         - r*qm3);
            pm4  = gam1*(1.0      - r*qm4);

            c2   = gam*p*r1;
            c2m1 = gam*(pm1*r1 + p*r1m1);
            c2m2 = gam*(pm2*r1 + p*r1m2);
            c2m3 = gam*(pm3*r1 + p*r1m3);
            c2m4 = gam*(pm4*r1 + p*r1m4);

            c    = sqrt(c2);
            cm1  = 0.5*c2m1/c;
            cm2  = 0.5*c2m2/c;
            cm3  = 0.5*c2m3/c;
            cm4  = 0.5*c2m4/c;

            un   = u*nx   + v*ny;
            unm1 = um1*nx + vm1*ny;
            unm2 = um2*nx + vm2*ny;
            unm3 = um3*nx + vm3*ny;
            unm4 = um4*nx + vm4*ny;
            
            if (absolute==1) {                
                rlam   = fabs(un)+c;
                signun = (un < 0) ? -1 : 1;   
                if (un==0) signun = 0;
                rlamm1 = signun*unm1+cm1;
                rlamm2 = signun*unm2+cm2;
                rlamm3 = signun*unm3+cm3;
                rlamm4 = signun*unm4+cm4;

                if (epslm>0) {
                    rlam = 0.5*(rlam*rlam/(epslm*c)+epslm*c);
                    rlamm1 = rlam*rlamm1/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam*rlamm2/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam*rlamm3/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam*rlamm4/(epslm*c) + 0.5*(1 - rlam*rlam/(epslm*c*epslm*c))*(epslm*cm4);
                }
            }
            else {
                rlam   = un-c;                
                rlamm1 = unm1-cm1;
                rlamm2 = unm2-cm2;
                rlamm3 = unm3-cm3;
                rlamm4 = unm4-cm4;
            }
            
            An[j*ng+i]  = rlam;
            Anm[0*ng+m] = rlamm1;
            Anm[1*ng+m] = rlamm2;    
            Anm[2*ng+m] = rlamm3;
            Anm[3*ng+m] = rlamm4;
                                                
        }
    }
    
}


void fbou_euler(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int ne)
{        
    
    double gam, epslm;
    double r, ru, rv, rE, ur, uru, urv, urE;           
    double nx, ny;        
    int    i, j, k, l, m, n, d, e, dim;
    
    double rlamm, rlamm1, rlamm2, rlamm3, rlamm4;            
    double rlamp, rlamp1, rlamp2, rlamp3, rlamp4;            
    double unu, uinf1, uinf2, uinf3, uinf4;
    
    double *an  = new double[ng*ne];
    double *An  = new double[ng*ne];
    double *anm = new double[ng*nc*ne];
    double *Anm = new double[ng*nc*ne];
        
    dim  = 2;
    gam  = param[0];
    epslm= param[1];        
                    
    switch (ib) {
        case 1 : // far-field boundary conditions
                                   
            getan_euler(an,anm,uhg,param,nl,1,ng,nc,ne);
            getan_euler(An,Anm,uhg,param,nl,0,ng,nc,ne);
            
            for (j=0; j<ne; j++) {     // loop over each element
                k = j*ng*nc;           // pointer to element j
                d = j*ng*dim;
                l = j*ng*nc*nc;                      
                for (i=0; i<ng; i++) { // loop over each gauss point    
                    m   = i + k;
                    e   = i + d;
                    n   = i + l;

                    r   = uhg[0*ng+m];
                    ru  = uhg[1*ng+m];
                    rv  = uhg[2*ng+m];
                    rE  = uhg[3*ng+m];

                    ur  = udg[0*ng+m];
                    uru = udg[1*ng+m];
                    urv = udg[2*ng+m];
                    urE = udg[3*ng+m];

                    nx  = nl[0*ng+e];              
                    ny  = nl[1*ng+e];   
                                        
                    rlamp = an[j*ng+i] + An[j*ng+i];
                    rlamm = an[j*ng+i] - An[j*ng+i];
                    
                    //fh = multiprod(an+An,(u-uh),[2,3],2) - multiprod(an-An,(uinf-uh),[2,3],2);
                    fh[0*ng+m] = rlamp*(ur - r)   - rlamm*(ui[0] - r);
                    fh[1*ng+m] = rlamp*(uru - ru) - rlamm*(ui[1] - ru);    
                    fh[2*ng+m] = rlamp*(urv - rv) - rlamm*(ui[2] - rv);
                    fh[3*ng+m] = rlamp*(urE - rE) - rlamm*(ui[3] - rE);
                                        
                    //fh_udg = an+An;
                    fh_u[0*ng+n] = rlamp;
                    fh_u[1*ng+n] = 0.0;
                    fh_u[2*ng+n] = 0.0;
                    fh_u[3*ng+n] = 0.0;        
                    fh_u[4*ng+n] = 0.0;
                    fh_u[5*ng+n] = rlamp;
                    fh_u[6*ng+n] = 0.0;
                    fh_u[7*ng+n] = 0.0;        
                    fh_u[8*ng+n] = 0.0;
                    fh_u[9*ng+n] = 0.0;
                    fh_u[10*ng+n] = rlamp;
                    fh_u[11*ng+n] = 0.0;        
                    fh_u[12*ng+n] = 0.0;
                    fh_u[13*ng+n] = 0.0;
                    fh_u[14*ng+n] = 0.0;
                    fh_u[15*ng+n] = rlamp;
                    
                    rlamp1 = anm[0*ng+m] + Anm[0*ng+m];
                    rlamp2 = anm[1*ng+m] + Anm[1*ng+m];
                    rlamp3 = anm[2*ng+m] + Anm[2*ng+m];
                    rlamp4 = anm[3*ng+m] + Anm[3*ng+m];
                    
                    rlamm1 = anm[0*ng+m] - Anm[0*ng+m];
                    rlamm2 = anm[1*ng+m] - Anm[1*ng+m];
                    rlamm3 = anm[2*ng+m] - Anm[2*ng+m];
                    rlamm4 = anm[3*ng+m] - Anm[3*ng+m];                    
                    
                    //fh_uh = multiprod(anm+Anm,(u-uh),[2,3],2) - multiprod(anm-Anm,(uinf-uh),[2,3],2) - 2*An;
                    fh_uh[0*ng+n]  = rlamp1*(ur - r)   - rlamm1*(ui[0] - r)  - 2*An[0*ng+m];
                    fh_uh[1*ng+n]  = rlamp1*(uru - ru) - rlamm1*(ui[1] - ru);
                    fh_uh[2*ng+n]  = rlamp1*(urv - rv) - rlamm1*(ui[2] - rv);
                    fh_uh[3*ng+n]  = rlamp1*(urE - rE) - rlamm1*(ui[3] - rE);        
                    fh_uh[4*ng+n]  = rlamp2*(ur - r)   - rlamm2*(ui[0] - r);
                    fh_uh[5*ng+n]  = rlamp2*(uru - ru) - rlamm2*(ui[1] - ru) - 2*An[0*ng+m];
                    fh_uh[6*ng+n]  = rlamp2*(urv - rv) - rlamm2*(ui[2] - rv);
                    fh_uh[7*ng+n]  = rlamp2*(urE - rE) - rlamm2*(ui[3] - rE);        
                    fh_uh[8*ng+n]  = rlamp3*(ur - r)   - rlamm3*(ui[0] - r);
                    fh_uh[9*ng+n]  = rlamp3*(uru - ru) - rlamm3*(ui[1] - ru);
                    fh_uh[10*ng+n] = rlamp3*(urv - rv) - rlamm3*(ui[2] - rv) - 2*An[0*ng+m];
                    fh_uh[11*ng+n] = rlamp3*(urE - rE) - rlamm3*(ui[3] - rE);              
                    fh_uh[12*ng+n] = rlamp4*(ur - r)   - rlamm4*(ui[0] - r);
                    fh_uh[13*ng+n] = rlamp4*(uru - ru) - rlamm4*(ui[1] - ru);
                    fh_uh[14*ng+n] = rlamp4*(urv - rv) - rlamm4*(ui[2] - rv);
                    fh_uh[15*ng+n] = rlamp4*(urE - rE) - rlamm4*(ui[3] - rE) - 2*An[0*ng+m];
                    
                }
            }
                                    
            break;
            
        case 2: // wall boundary conditions                        
            
            for (j=0; j<ne; j++) {     // loop over each element
                k = j*ng*nc;           // pointer to element j
                d = j*ng*dim;
                l = j*ng*nc*nc;                      
                for (i=0; i<ng; i++) { // loop over each gauss point    
                    m   = i + k;
                    e   = i + d;
                    n   = i + l;

                    r   = uhg[0*ng+m];
                    ru  = uhg[1*ng+m];
                    rv  = uhg[2*ng+m];
                    rE  = uhg[3*ng+m];

                    ur  = udg[0*ng+m];
                    uru = udg[1*ng+m];
                    urv = udg[2*ng+m];
                    urE = udg[3*ng+m];

                    nx  = nl[0*ng+e];              
                    ny  = nl[1*ng+e];   
                    
                    unu   = uru*nx + urv*ny;
                    uinf1 = ur;
                    uinf2 = uru - unu*nx;
                    uinf3 = urv - unu*ny;
                    uinf4 = urE;
                    
                    // fh = uinf - uh;
                    fh[0*ng+m] = uinf1 - r;
                    fh[1*ng+m] = uinf2 - ru;    
                    fh[2*ng+m] = uinf3 - rv;
                    fh[3*ng+m] = uinf4 - rE;
                    
                    // fh_udg = uinfu;
                    fh_u[0*ng+n] = 1.0;
                    fh_u[1*ng+n] = 0.0;
                    fh_u[2*ng+n] = 0.0;
                    fh_u[3*ng+n] = 0.0;        
                    fh_u[4*ng+n] = 0.0;
                    fh_u[5*ng+n] = 1.0 - nx*nx;
                    fh_u[6*ng+n] = 0.0 - nx*ny;
                    fh_u[7*ng+n] = 0.0;        
                    fh_u[8*ng+n] = 0.0;
                    fh_u[9*ng+n] = 0.0 - nx*ny;
                    fh_u[10*ng+n] = 1.0 - ny*ny;
                    fh_u[11*ng+n] = 0.0;        
                    fh_u[12*ng+n] = 0.0;
                    fh_u[13*ng+n] = 0.0;
                    fh_u[14*ng+n] = 0.0;
                    fh_u[15*ng+n] = 1.0;
                    
                    //fh_uh 
                    fh_uh[0*ng+n]  = -1.0;
                    fh_uh[1*ng+n]  = 0.0;
                    fh_uh[2*ng+n]  = 0.0;
                    fh_uh[3*ng+n]  = 0.0;        
                    fh_uh[4*ng+n]  = 0.0;
                    fh_uh[5*ng+n]  = -1.0;
                    fh_uh[6*ng+n]  = 0.0;
                    fh_uh[7*ng+n]  = 0.0;        
                    fh_uh[8*ng+n]  = 0.0;
                    fh_uh[9*ng+n]  = 0.0;
                    fh_uh[10*ng+n] = -1.0;
                    fh_uh[11*ng+n] = 0.0;              
                    fh_uh[12*ng+n] = 0.0;
                    fh_uh[13*ng+n] = 0.0;
                    fh_uh[14*ng+n] = 0.0;
                    fh_uh[15*ng+n] = -1.0;
                    
                }
            }
            
            break;
            
        case 3: // Prescribed pressure
                                            
            getan_euler(an,anm,uhg,param,nl,1,ng,nc,ne);
            getan_euler(An,Anm,uhg,param,nl,0,ng,nc,ne);
            
            for (j=0; j<ne; j++) {     // loop over each element
                k = j*ng*nc;           // pointer to element j
                d = j*ng*dim;
                l = j*ng*nc*nc;                      
                for (i=0; i<ng; i++) { // loop over each gauss point    
                    m   = i + k;
                    e   = i + d;
                    n   = i + l;

                    r   = uhg[0*ng+m];
                    ru  = uhg[1*ng+m];
                    rv  = uhg[2*ng+m];
                    rE  = uhg[3*ng+m];

                    ur  = udg[0*ng+m];
                    uru = udg[1*ng+m];
                    urv = udg[2*ng+m];
                    urE = udg[3*ng+m];

                    nx  = nl[0*ng+e];              
                    ny  = nl[1*ng+e];   
                                        
                    rlamp = an[j*ng+i] + An[j*ng+i];
                    rlamm = an[j*ng+i] - An[j*ng+i];
                    
                    uinf1 = ur;
                    uinf2 = uru;
                    uinf3 = urv;
                    uinf4 = ui[3]/(gam-1) + 0.5*uru*uru/ur + 0.5*urv*urv/ur;
                               
                    // fh = multiprod(an+An,(u-uh),[2,3],2) - multiprod(an-An,(uinf-uh),[2,3],2);
                    fh[0*ng+m] = rlamp*(ur - r)   - rlamm*(uinf1 - r);
                    fh[1*ng+m] = rlamp*(uru - ru) - rlamm*(uinf2 - ru);    
                    fh[2*ng+m] = rlamp*(urv - rv) - rlamm*(uinf3 - rv);
                    fh[3*ng+m] = rlamp*(urE - rE) - rlamm*(uinf4 - rE);
                          
                    //uinfu = multiprod(one,eye(nch),2,[1,2]);
                    //uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zer];
        
                    // fh_udg = an+An - multiprod(an-An,uinfu,[2,3]);
                    fh_u[0*ng+n] = rlamp - rlamm;
                    fh_u[1*ng+n] = 0.0;
                    fh_u[2*ng+n] = 0.0;
                    fh_u[3*ng+n] = rlamm*(0.5*(uru*uru+urv*urv)/(ur*ur));        
                    fh_u[4*ng+n] = 0.0;
                    fh_u[5*ng+n] = rlamp - rlamm;
                    fh_u[6*ng+n] = 0.0;
                    fh_u[7*ng+n] = -rlamm*uru/ur;        
                    fh_u[8*ng+n] = 0.0;
                    fh_u[9*ng+n] = 0.0;
                    fh_u[10*ng+n] = rlamp - rlamm;
                    fh_u[11*ng+n] = -rlamm*urv/ur;        
                    fh_u[12*ng+n] = 0.0;
                    fh_u[13*ng+n] = 0.0;
                    fh_u[14*ng+n] = 0.0;
                    fh_u[15*ng+n] = rlamp - rlamm;
                    
                    rlamp1 = anm[0*ng+m] + Anm[0*ng+m];
                    rlamp2 = anm[1*ng+m] + Anm[1*ng+m];
                    rlamp3 = anm[2*ng+m] + Anm[2*ng+m];
                    rlamp4 = anm[3*ng+m] + Anm[3*ng+m];
                    
                    rlamm1 = anm[0*ng+m] - Anm[0*ng+m];
                    rlamm2 = anm[1*ng+m] - Anm[1*ng+m];
                    rlamm3 = anm[2*ng+m] - Anm[2*ng+m];
                    rlamm4 = anm[3*ng+m] - Anm[3*ng+m];                    
                    
                    //fh_uh = multiprod(anm+Anm,(u-uh),[2,3],2) - multiprod(anm-Anm,(uinf-uh),[2,3],2) - 2*An;
                    fh_uh[0*ng+n]  = rlamp1*(ur - r)   - rlamm1*(uinf1 - r)  - 2*An[0*ng+m];
                    fh_uh[1*ng+n]  = rlamp1*(uru - ru) - rlamm1*(uinf2 - ru);
                    fh_uh[2*ng+n]  = rlamp1*(urv - rv) - rlamm1*(uinf3 - rv);
                    fh_uh[3*ng+n]  = rlamp1*(urE - rE) - rlamm1*(uinf4 - rE);        
                    fh_uh[4*ng+n]  = rlamp2*(ur - r)   - rlamm2*(uinf1 - r);
                    fh_uh[5*ng+n]  = rlamp2*(uru - ru) - rlamm2*(uinf2 - ru) - 2*An[0*ng+m];
                    fh_uh[6*ng+n]  = rlamp2*(urv - rv) - rlamm2*(uinf3 - rv);
                    fh_uh[7*ng+n]  = rlamp2*(urE - rE) - rlamm2*(uinf4 - rE);        
                    fh_uh[8*ng+n]  = rlamp3*(ur - r)   - rlamm3*(uinf1 - r);
                    fh_uh[9*ng+n]  = rlamp3*(uru - ru) - rlamm3*(uinf2 - ru);
                    fh_uh[10*ng+n] = rlamp3*(urv - rv) - rlamm3*(uinf3 - rv) - 2*An[0*ng+m];
                    fh_uh[11*ng+n] = rlamp3*(urE - rE) - rlamm3*(uinf4 - rE);              
                    fh_uh[12*ng+n] = rlamp4*(ur - r)   - rlamm4*(uinf1 - r);
                    fh_uh[13*ng+n] = rlamp4*(uru - ru) - rlamm4*(uinf2 - ru);
                    fh_uh[14*ng+n] = rlamp4*(urv - rv) - rlamm4*(uinf3 - rv);
                    fh_uh[15*ng+n] = rlamp4*(urE - rE) - rlamm4*(uinf4 - rE) - 2*An[0*ng+m];
                    
                }
            }                                    
            
            break;
            
        default :
                        
            printf("This error is in %s on line %d\n",__FILE__, __LINE__);
            printf("Boundary condition %d is not implemented yet.", ib);            
            exit(-1);
            
            break;
                        
    }
    
    delete[] an; delete[] An; delete[] anm; delete[] Anm;
    
}
