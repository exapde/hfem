void flux_ns(double *f, double *f_u, double *pg, double *udg, double *param,
        double time, int ng, int nch, int nc, int nd)
{
    
    double gam, gam1, c23, c12, Re, Pr, Minf, Re1, M2, fc, zero;
    double r, ru, rv, rE;    
    double r1, u, v, E, uu, vv, uv, q, p, h;    
    int    i, j, k, l, m, n, a, b;
    double rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u_r, u_ru, v_r, v_rv;        
    double ux, vx, Ex, qx, px, Tx, uy, vy, Ey, qy, py, Ty;
    double txx, txy, tyy, txx_r, txx_ru, txx_rv, txx_rE;
    double txx_rx, txx_rux, txx_rvx, txx_rEx, txx_ry, txx_ruy, txx_rvy, txx_rEy;
    double txy_r, txy_ru, txy_rv, txy_rE;
    double txy_rx, txy_rux, txy_rvx, txy_rEx, txy_ry, txy_ruy, txy_rvy, txy_rEy;
    double tyy_r, tyy_ru, tyy_rv, tyy_rE;
    double tyy_rx, tyy_rux, tyy_rvx, tyy_rEx, tyy_ry, tyy_ruy, tyy_rvy, tyy_rEy;
    double Tx_r, Tx_ru, Tx_rv, Tx_rE, Tx_rx, Tx_rux, Tx_rvx, Tx_rEx;
    double Ty_r, Ty_ru, Ty_rv, Ty_rE, Ty_ry, Ty_ruy, Ty_rvy, Ty_rEy;    
    
    zero = 0.0;    
    gam  = param[0];      
    Re   = param[2];
    Pr   = param[3];
    Minf = param[4];
    gam1 = gam - 1.0;
    Re1  = 1/Re;
    M2   = Minf*Minf; 
    c23  = 2.0/3.0;
    c12  = 1.0/2.0;
    fc   = 1.0/(gam1*M2*Re*Pr);    

    for (i=0; i<ng; i++) { 
        r   = udg[0*ng+i];
        ru  = udg[1*ng+i];
        rv  = udg[2*ng+i];
        rE  = udg[3*ng+i];

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

        /* inviscid flux vectors */
        f[0*ng+i]  = ru;
        f[1*ng+i]  = ru*u + p;    
        f[2*ng+i]  = rv*u;
        f[3*ng+i]  = ru*h;        
        f[4*ng+i]  = rv;
        f[5*ng+i]  = ru*v;
        f[6*ng+i]  = rv*v + p;    
        f[7*ng+i]  = rv*h;          

        /* derivatives with respect to r */
        f_u[0*ng+i] = 0.0;
        f_u[1*ng+i] = 0.5*((gam-3)*uu+gam1*vv);
        f_u[2*ng+i] = -uv;
        f_u[3*ng+i] = 2*gam1*u*q-gam*E*u;
        f_u[4*ng+i] = 0.0;
        f_u[5*ng+i] = -uv;
        f_u[6*ng+i] =  0.5*((gam-3)*vv+gam1*uu);
        f_u[7*ng+i] = 2*gam1*v*q-gam*E*v;        
        
        /* derivatives with respect to ru */
        f_u[8*ng+i] = 1.0;
        f_u[9*ng+i] = (3-gam)*u;
        f_u[10*ng+i] = v;
        f_u[11*ng+i] = gam*E-0.5*gam1*(3*uu+vv);        
        f_u[12*ng+i] = 0.0;
        f_u[13*ng+i] = v;
        f_u[14*ng+i] = -gam1*u;
        f_u[15*ng+i] = -gam1*uv;        
        
        /* derivatives with respect to rv */
        f_u[16*ng+i] = 0.0;
        f_u[17*ng+i] = -gam1*v;
        f_u[18*ng+i] = u;
        f_u[19*ng+i] = -gam1*uv;        
        f_u[20*ng+i] = 1.0;
        f_u[21*ng+i] = u;
        f_u[22*ng+i] = (3-gam)*v;
        f_u[23*ng+i] = gam*E-0.5*gam1*(3*vv+uu);        
        
        /* derivatives with respect to rE */
        f_u[24*ng+i] = 0.0;
        f_u[25*ng+i] = gam1;
        f_u[26*ng+i] = 0.0;
        f_u[27*ng+i] = gam*u;                        
        f_u[28*ng+i] = 0.0;
        f_u[29*ng+i] = 0.0;
        f_u[30*ng+i] = gam1;
        f_u[31*ng+i] = gam*v;

        rx   = udg[4*ng+i];
        rux  = udg[5*ng+i];
        rvx  = udg[6*ng+i];
        rEx  = udg[7*ng+i];

        ry   = udg[8*ng+i];
        ruy  = udg[9*ng+i];
        rvy  = udg[10*ng+i];
        rEy  = udg[11*ng+i];

        u_r  = -u*r1;
        u_ru =  r1;
        v_r  = -v*r1;
        v_rv =  r1;

        ux  = (rux - rx*u)*r1;
        vx  = (rvx - rx*v)*r1;
        Ex  = (rEx - rx*E)*r1;
        qx  = u*ux + v*vx;
        px  = gam1*(rEx - rx*q - r*qx);
        Tx  = gam*M2*(px*r - p*rx)*(r1*r1);

        uy  = (ruy - ry*u)*r1;
        vy  = (rvy - ry*v)*r1;
        Ey  = (rEy - ry*E)*r1;
        qy  = u*uy + v*vy;
        py  = gam1*(rEy - ry*q - r*qy);
        Ty  = gam*M2*(py*r - p*ry)*(r1*r1);

        txx = Re1*c23*(2*ux - vy);
        txy = Re1*(uy + vx);
        tyy = Re1*c23*(2*vy - ux);

        txx_r  =  Re1*c23*((4*ru*rx-2*rv*ry)-r*(2*rux-rvy))*(r1*r1*r1);
        txx_ru = -Re1*c23*2*rx*r1*r1;
        txx_rv =  Re1*c23*ry*r1*r1;
        txx_rE =  zero;

        txx_rx  = -Re1*c23*2*ru*r1*r1;
        txx_rux =  Re1*c23*2*r1;
        txx_rvx =  zero;
        txx_rEx =  zero;

        txx_ry  =  Re1*c23*rv*r1*r1;
        txx_ruy =  zero;
        txx_rvy = -Re1*c23*r1;
        txx_rEy =  zero;

        txy_r  =  Re1*(2*(ru*ry+rv*rx)-r*(ruy+rvx))*r1*r1*r1;
        txy_ru = -Re1*ry*r1*r1;
        txy_rv = -Re1*rx*r1*r1;
        txy_rE =  zero;

        txy_rx  = -Re1*rv*r1*r1;
        txy_rux =  zero;
        txy_rvx =  Re1*r1;
        txy_rEx =  zero;

        txy_ry  = -Re1*ru*r1*r1;
        txy_ruy =  Re1*r1;
        txy_rvy =  zero;
        txy_rEy =  zero;

        tyy_r  =  Re1*c23*((4*rv*ry-2*ru*rx)-r*(2*rvy-rux))*r1*r1*r1;
        tyy_ru =  Re1*c23*rx*r1*r1;
        tyy_rv = -Re1*c23*2*ry*r1*r1;
        tyy_rE =  zero;

        tyy_rx  =  Re1*c23*ru*r1*r1;
        tyy_rux = -Re1*c23*r1;
        tyy_rvx =  zero;
        tyy_rEx =  zero;

        tyy_ry  = -Re1*c23*2*rv*r1*r1;
        tyy_ruy =  zero;
        tyy_rvy =  Re1*c23*2*r1;
        tyy_rEy =  zero;

        Tx_r  = -M2*gam*gam1*(rEx*r*r-2*rux*r*ru-2*rvx*r*rv-2*rE*rx*r+3*rx*(ru*ru+rv*rv))*(r1*r1*r1*r1);
        Tx_ru = -M2*gam*gam1*(r*rux-2*ru*rx)*(r1*r1*r1);
        Tx_rv = -M2*gam*gam1*(r*rvx-2*rv*rx)*(r1*r1*r1);
        Tx_rE = -M2*gam*gam1*rx*(r1*r1);

        Tx_rx  =  M2*gam*gam1*(ru*ru+rv*rv-r*rE)*(r1*r1*r1);
        Tx_rux = -M2*gam*gam1*ru*(r1*r1);
        Tx_rvx = -M2*gam*gam1*rv*(r1*r1);
        Tx_rEx =  M2*gam*gam1*r1;

        Ty_r  = -M2*gam*gam1*(rEy*r*r-2*ruy*r*ru-2*rvy*r*rv-2*rE*ry*r+3*ry*(ru*ru+rv*rv))*(r1*r1*r1*r1);
        Ty_ru = -M2*gam*gam1*(r*ruy-2*ru*ry)*(r1*r1*r1);
        Ty_rv = -M2*gam*gam1*(r*rvy-2*rv*ry)*(r1*r1*r1);
        Ty_rE = -M2*gam*gam1*ry*(r1*r1);

        Ty_ry  = Tx_rx;
        Ty_ruy = Tx_rux;
        Ty_rvy = Tx_rvx;
        Ty_rEy = Tx_rEx;

        /* inviscid+viscous flux vectors */
        f[0*ng+i]  += 0.0;
        f[1*ng+i]  += txx;    
        f[2*ng+i]  += txy;
        f[3*ng+i]  += u*txx + v*txy + fc*Tx;
        f[4*ng+i]  += 0.0;
        f[5*ng+i]  += txy;
        f[6*ng+i]  += tyy;    
        f[7*ng+i]  += u*txy + v*tyy + fc*Ty;          

        /* derivatives with respect to r */
        f_u[0*ng+i] += 0.0;
        f_u[1*ng+i] += txx_r;
        f_u[2*ng+i] += txy_r;
        f_u[3*ng+i] += u_r*txx  + u*txx_r  + v_r*txy  + v*txy_r  + fc*Tx_r;                              
        f_u[4*ng+i] += 0.0;
        f_u[5*ng+i] += txy_r;
        f_u[6*ng+i] += tyy_r;
        f_u[7*ng+i] += u_r*txy  + u*txy_r  + v_r*tyy  + v*tyy_r  + fc*Ty_r;        
        
        /* derivatives with respect to ru */
        f_u[8*ng+i] += 0.0;
        f_u[9*ng+i] += txx_ru;
        f_u[10*ng+i] += txy_ru;
        f_u[11*ng+i] += u_ru*txx + u*txx_ru + v*txy_ru + fc*Tx_ru;  
        f_u[12*ng+i] += 0.0;
        f_u[13*ng+i] += txy_ru;
        f_u[14*ng+i] += tyy_ru; 
        f_u[15*ng+i] += u_ru*txy + u*txy_ru + v*tyy_ru + fc*Ty_ru;        
        
        /* derivatives with respect to rv */
        f_u[16*ng+i] += 0.0;
        f_u[17*ng+i] += txx_rv;
        f_u[18*ng+i] += txy_rv;
        f_u[19*ng+i] += u*txx_rv + v_rv*txy + v*txy_rv + fc*Tx_rv;  
        f_u[20*ng+i] += 0.0;
        f_u[21*ng+i] += txy_rv;
        f_u[22*ng+i] += tyy_rv;            
        f_u[23*ng+i] += u*txy_rv + v_rv*tyy + v*tyy_rv + fc*Ty_rv;       
        
        /* derivatives with respect to rE */
        f_u[24*ng+i] += 0.0;
        f_u[25*ng+i] += txx_rE;
        f_u[26*ng+i] += txy_rE;
        f_u[27*ng+i] += u*txx_rE + v*txy_rE + fc*Tx_rE;                                        
        f_u[28*ng+i] += 0.0;
        f_u[29*ng+i] += txy_rE;
        f_u[30*ng+i] += tyy_rE;
        f_u[31*ng+i] += u*txy_rE + v*tyy_rE + fc*Ty_rE;
        
        /* derivatives with respect to rx */
        f_u[32*ng+i] = 0.0;
        f_u[33*ng+i] = txx_rx;
        f_u[34*ng+i] = txy_rx;
        f_u[35*ng+i] = u*txx_rx  + v*txy_rx  + fc*Tx_rx;                              
        f_u[36*ng+i] = 0.0;
        f_u[37*ng+i] = txy_rx; 
        f_u[38*ng+i] = tyy_rx; 
        f_u[39*ng+i] = u*txy_rx  + v*tyy_rx;
        
        /* derivatives with respect to rux */
        f_u[40*ng+i] = 0.0;
        f_u[41*ng+i] = txx_rux;
        f_u[42*ng+i] = txy_rux;
        f_u[43*ng+i] = u*txx_rux + v*txy_rux + fc*Tx_rux;        
        f_u[44*ng+i] = 0.0;
        f_u[45*ng+i] = txy_rux; 
        f_u[46*ng+i] = tyy_rux; 
        f_u[47*ng+i] = u*txy_rux + v*tyy_rux;
        
        /* derivatives with respect to rvx */
        f_u[48*ng+i] = 0.0;
        f_u[49*ng+i] = txx_rvx;
        f_u[50*ng+i] = txy_rvx;
        f_u[51*ng+i] = u*txx_rvx + v*txy_rvx + fc*Tx_rvx;
        f_u[52*ng+i] = 0.0;
        f_u[53*ng+i] = txy_rvx; 
        f_u[54*ng+i] = tyy_rvx; 
        f_u[55*ng+i] = u*txy_rvx + v*tyy_rvx;
        
        /* derivatives with respect to rEx */
        f_u[56*ng+i] = 0.0;
        f_u[57*ng+i] = txx_rEx;
        f_u[58*ng+i] = txy_rEx;
        f_u[59*ng+i] = u*txx_rEx + v*txy_rEx + fc*Tx_rEx;
        f_u[60*ng+i] = 0.0;
        f_u[61*ng+i] = txy_rEx;
        f_u[62*ng+i] = tyy_rEx;
        f_u[63*ng+i] = u*txy_rEx + v*tyy_rEx; 
        
        /* derivatives with respect to ry */
        f_u[64*ng+i] = 0.0;
        f_u[65*ng+i] = txx_ry; 
        f_u[66*ng+i] = txy_ry;
        f_u[67*ng+i] = u*txx_ry  + v*txy_ry;
        f_u[68*ng+i] = 0.0;
        f_u[69*ng+i] = txy_ry; 
        f_u[70*ng+i] = tyy_ry; 
        f_u[71*ng+i] = u*txy_ry  + v*tyy_ry  + fc*Ty_ry;
        
        /* derivatives with respect to ruy */
        f_u[72*ng+i] = 0.0;
        f_u[73*ng+i] = txx_ruy;
        f_u[74*ng+i] = txy_ruy;
        f_u[75*ng+i] = u*txx_ruy + v*txy_ruy;
        f_u[76*ng+i] = 0.0;
        f_u[77*ng+i] = txy_ruy; 
        f_u[78*ng+i] = tyy_ruy; 
        f_u[79*ng+i] = u*txy_ruy + v*tyy_ruy + fc*Ty_ruy;
        
        /* derivatives with respect to rvy */
        f_u[80*ng+i] = 0.0;
        f_u[81*ng+i] = txx_rvy; 
        f_u[82*ng+i] = txy_rvy; 
        f_u[83*ng+i] = u*txx_rvy + v*txy_rvy;
        f_u[84*ng+i] = 0.0;
        f_u[85*ng+i] = txy_rvy;
        f_u[86*ng+i] = tyy_rvy;
        f_u[87*ng+i] = u*txy_rvy + v*tyy_rvy + fc*Ty_rvy;
        
        /* derivatives with respect to rEy */
        f_u[88*ng+i] = 0.0;
        f_u[89*ng+i] = txx_rEy; 
        f_u[90*ng+i] = txy_rEy; 
        f_u[91*ng+i] = u*txx_rEy + v*txy_rEy;                                                
        f_u[92*ng+i] = 0.0;
        f_u[93*ng+i] = txy_rEy;
        f_u[94*ng+i] = tyy_rEy;
        f_u[95*ng+i] = u*txy_rEy + v*tyy_rEy + fc*Ty_rEy;            
    }

}

void fhat_ns(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *uhg, double *nlg, 
          double *param, double time, int ng, int nch, int nc, int nd)
{
    
    double gam, epslm, gam1, signun, zero;
    double r, ru, rv, rE, ur, uru, urv, urE;    
    double r1, u, v, E, uu, vv, uv, q, p, h;     
    double nx, ny;
    double fx[4], fy[4], fx_uh[16], fy_uh[16], fx_q[32], fy_q[32];
    double um1, um2, um3, um4, vm1, vm2, vm3, vm4;
    double r1m1, r1m2, r1m3, r1m4, qm1, qm2, qm3, qm4;
    double pm1, pm2, pm3, pm4, c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4, un, unm1, unm2, unm3, unm4;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;            
    
    int    i;
    
    double rx, rux, rvx, rEx, ry, ruy, rvy, rEy, u_r, u_ru, v_r, v_rv;        
    double ux, vx, Ex, qx, px, Tx, uy, vy, Ey, qy, py, Ty;
    double txx, txy, tyy, txx_r, txx_ru, txx_rv, txx_rE;
    double txx_rx, txx_rux, txx_rvx, txx_rEx, txx_ry, txx_ruy, txx_rvy, txx_rEy;
    double txy_r, txy_ru, txy_rv, txy_rE;
    double txy_rx, txy_rux, txy_rvx, txy_rEx, txy_ry, txy_ruy, txy_rvy, txy_rEy;
    double tyy_r, tyy_ru, tyy_rv, tyy_rE;
    double tyy_rx, tyy_rux, tyy_rvx, tyy_rEx, tyy_ry, tyy_ruy, tyy_rvy, tyy_rEy;
    double Tx_r, Tx_ru, Tx_rv, Tx_rE, Tx_rx, Tx_rux, Tx_rvx, Tx_rEx;
    double Ty_r, Ty_ru, Ty_rv, Ty_rE, Ty_ry, Ty_ruy, Ty_rvy, Ty_rEy;    
    
    double Re, Pr, Minf, Re1, M2, c23, c12, fc;
            
    zero = 0.0;        
    gam  = param[0];      
    epslm= param[1];       
    Re   = param[2];
    Pr   = param[3];
    Minf = param[4];    
    gam1 = gam - 1.0;
    Re1  = 1/Re;
    M2   = Minf*Minf; 
    c23  = 2.0/3.0;
    c12  = 1.0/2.0;
    fc   = 1.0/(gam1*M2*Re*Pr);    
    
    for (i=0; i<ng; i++) { 
        r   = uhg[0*ng+i];
        ru  = uhg[1*ng+i];
        rv  = uhg[2*ng+i];
        rE  = uhg[3*ng+i];

        ur  = udg[0*ng+i];
        uru = udg[1*ng+i];
        urv = udg[2*ng+i];
        urE = udg[3*ng+i];

        nx  = nlg[0*ng+i];              
        ny  = nlg[1*ng+i];   

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
        
        rlam = param[5];
        rlamm1 = 0.0;
        rlamm2 = 0.0;
        rlamm3 = 0.0;
        rlamm4 = 0.0;
        
        rx   = udg[4*ng+i];
        rux  = udg[5*ng+i];
        rvx  = udg[6*ng+i];
        rEx  = udg[7*ng+i];

        ry   = udg[8*ng+i];
        ruy  = udg[9*ng+i];
        rvy  = udg[10*ng+i];
        rEy  = udg[11*ng+i];

        u_r  = -u*r1;
        u_ru =  r1;
        v_r  = -v*r1;
        v_rv =  r1;

        ux  = (rux - rx*u)*r1;
        vx  = (rvx - rx*v)*r1;
        Ex  = (rEx - rx*E)*r1;
        qx  = u*ux + v*vx;
        px  = gam1*(rEx - rx*q - r*qx);
        Tx  = gam*M2*(px*r - p*rx)*(r1*r1);

        uy  = (ruy - ry*u)*r1;
        vy  = (rvy - ry*v)*r1;
        Ey  = (rEy - ry*E)*r1;
        qy  = u*uy + v*vy;
        py  = gam1*(rEy - ry*q - r*qy);
        Ty  = gam*M2*(py*r - p*ry)*(r1*r1);

        txx = Re1*c23*(2*ux - vy);
        txy = Re1*(uy + vx);
        tyy = Re1*c23*(2*vy - ux);

        txx_r  =  Re1*c23*((4*ru*rx-2*rv*ry)-r*(2*rux-rvy))*(r1*r1*r1);
        txx_ru = -Re1*c23*2*rx*r1*r1;
        txx_rv =  Re1*c23*ry*r1*r1;
        txx_rE =  zero;

        txx_rx  = -Re1*c23*2*ru*r1*r1;
        txx_rux =  Re1*c23*2*r1;
        txx_rvx =  zero;
        txx_rEx =  zero;

        txx_ry  =  Re1*c23*rv*r1*r1;
        txx_ruy =  zero;
        txx_rvy = -Re1*c23*r1;
        txx_rEy =  zero;

        txy_r  =  Re1*(2*(ru*ry+rv*rx)-r*(ruy+rvx))*r1*r1*r1;
        txy_ru = -Re1*ry*r1*r1;
        txy_rv = -Re1*rx*r1*r1;
        txy_rE =  zero;

        txy_rx  = -Re1*rv*r1*r1;
        txy_rux =  zero;
        txy_rvx =  Re1*r1;
        txy_rEx =  zero;

        txy_ry  = -Re1*ru*r1*r1;
        txy_ruy =  Re1*r1;
        txy_rvy =  zero;
        txy_rEy =  zero;

        tyy_r  =  Re1*c23*((4*rv*ry-2*ru*rx)-r*(2*rvy-rux))*r1*r1*r1;
        tyy_ru =  Re1*c23*rx*r1*r1;
        tyy_rv = -Re1*c23*2*ry*r1*r1;
        tyy_rE =  zero;

        tyy_rx  =  Re1*c23*ru*r1*r1;
        tyy_rux = -Re1*c23*r1;
        tyy_rvx =  zero;
        tyy_rEx =  zero;

        tyy_ry  = -Re1*c23*2*rv*r1*r1;
        tyy_ruy =  zero;
        tyy_rvy =  Re1*c23*2*r1;
        tyy_rEy =  zero;

        Tx_r  = -M2*gam*gam1*(rEx*r*r-2*rux*r*ru-2*rvx*r*rv-2*rE*rx*r+3*rx*(ru*ru+rv*rv))*(r1*r1*r1*r1);
        Tx_ru = -M2*gam*gam1*(r*rux-2*ru*rx)*(r1*r1*r1);
        Tx_rv = -M2*gam*gam1*(r*rvx-2*rv*rx)*(r1*r1*r1);
        Tx_rE = -M2*gam*gam1*rx*(r1*r1);

        Tx_rx  =  M2*gam*gam1*(ru*ru+rv*rv-r*rE)*(r1*r1*r1);
        Tx_rux = -M2*gam*gam1*ru*(r1*r1);
        Tx_rvx = -M2*gam*gam1*rv*(r1*r1);
        Tx_rEx =  M2*gam*gam1*r1;

        Ty_r  = -M2*gam*gam1*(rEy*r*r-2*ruy*r*ru-2*rvy*r*rv-2*rE*ry*r+3*ry*(ru*ru+rv*rv))*(r1*r1*r1*r1);
        Ty_ru = -M2*gam*gam1*(r*ruy-2*ru*ry)*(r1*r1*r1);
        Ty_rv = -M2*gam*gam1*(r*rvy-2*rv*ry)*(r1*r1*r1);
        Ty_rE = -M2*gam*gam1*ry*(r1*r1);

        Ty_ry  = Tx_rx;
        Ty_ruy = Tx_rux;
        Ty_rvy = Tx_rvx;
        Ty_rEy = Tx_rEx;

        fx[0]  += 0.0;
        fx[1]  += txx;    
        fx[2]  += txy;
        fx[3]  += u*txx + v*txy + fc*Tx;

        fy[0]  += 0.0;
        fy[1]  += txy;
        fy[2]  += tyy;    
        fy[3]  += u*txy + v*tyy + fc*Ty;          

        fh[0*ng+i]  = fx[0]*nx + fy[0]*ny + rlam*(ur - r);
        fh[1*ng+i]  = fx[1]*nx + fy[1]*ny + rlam*(uru - ru);    
        fh[2*ng+i]  = fx[2]*nx + fy[2]*ny + rlam*(urv - rv);
        fh[3*ng+i]  = fx[3]*nx + fy[3]*ny + rlam*(urE - rE);            

        fh_u[0*ng+i] = rlam;
        fh_u[1*ng+i] = 0.0;
        fh_u[2*ng+i] = 0.0;
        fh_u[3*ng+i] = 0.0;        
        fh_u[4*ng+i] = 0.0;
        fh_u[5*ng+i] = rlam;
        fh_u[6*ng+i] = 0.0;
        fh_u[7*ng+i] = 0.0;        
        fh_u[8*ng+i] = 0.0;
        fh_u[9*ng+i] = 0.0;
        fh_u[10*ng+i] = rlam;
        fh_u[11*ng+i] = 0.0;        
        fh_u[12*ng+i] = 0.0;
        fh_u[13*ng+i] = 0.0;
        fh_u[14*ng+i] = 0.0;
        fh_u[15*ng+i] = rlam;

        fx_uh[0] += 0.0;
        fx_uh[1] += txx_r;
        fx_uh[2] += txy_r;
        fx_uh[3] += u_r*txx  + u*txx_r  + v_r*txy  + v*txy_r  + fc*Tx_r;                              
        fx_uh[4] += 0.0;
        fx_uh[5] += txx_ru;
        fx_uh[6] += txy_ru;
        fx_uh[7] += u_ru*txx + u*txx_ru + v*txy_ru + fc*Tx_ru;        
        fx_uh[8] += 0.0;
        fx_uh[9] += txx_rv;
        fx_uh[10] += txy_rv;
        fx_uh[11] += u*txx_rv + v_rv*txy + v*txy_rv + fc*Tx_rv;        
        fx_uh[12] += 0.0;
        fx_uh[13] += txx_rE;
        fx_uh[14] += txy_rE;
        fx_uh[15] += u*txx_rE + v*txy_rE + fc*Tx_rE;

        fy_uh[0] += 0.0;
        fy_uh[1] += txy_r;
        fy_uh[2] += tyy_r;
        fy_uh[3] += u_r*txy  + u*txy_r  + v_r*tyy  + v*tyy_r  + fc*Ty_r;        
        fy_uh[4] += 0.0;
        fy_uh[5] += txy_ru;
        fy_uh[6] += tyy_ru; 
        fy_uh[7] += u_ru*txy + u*txy_ru + v*tyy_ru + fc*Ty_ru;        
        fy_uh[8] += 0.0;
        fy_uh[9] += txy_rv;
        fy_uh[10] += tyy_rv;            
        fy_uh[11] += u*txy_rv + v_rv*tyy + v*tyy_rv + fc*Ty_rv;       
        fy_uh[12] += 0.0;
        fy_uh[13] += txy_rE;
        fy_uh[14] += tyy_rE;
        fy_uh[15] += u*txy_rE + v*tyy_rE + fc*Ty_rE;

        fh_uh[0*ng+i] = fx_uh[0]*nx + fy_uh[0]*ny + rlamm1*(ur - r) - rlam;
        fh_uh[1*ng+i] = fx_uh[1]*nx + fy_uh[1]*ny + rlamm1*(uru - ru);
        fh_uh[2*ng+i] = fx_uh[2]*nx + fy_uh[2]*ny + rlamm1*(urv - rv);
        fh_uh[3*ng+i] = fx_uh[3]*nx + fy_uh[3]*ny + rlamm1*(urE - rE);        
        fh_uh[4*ng+i] = fx_uh[4]*nx + fy_uh[4]*ny + rlamm2*(ur - r);
        fh_uh[5*ng+i] = fx_uh[5]*nx + fy_uh[5]*ny + rlamm2*(uru - ru) - rlam;
        fh_uh[6*ng+i] = fx_uh[6]*nx + fy_uh[6]*ny + rlamm2*(urv - rv);
        fh_uh[7*ng+i] = fx_uh[7]*nx + fy_uh[7]*ny + rlamm2*(urE - rE);        
        fh_uh[8*ng+i] = fx_uh[8]*nx + fy_uh[8]*ny + rlamm3*(ur - r);
        fh_uh[9*ng+i] = fx_uh[9]*nx + fy_uh[9]*ny + rlamm3*(uru - ru);
        fh_uh[10*ng+i] = fx_uh[10]*nx + fy_uh[10]*ny + rlamm3*(urv - rv) - rlam;
        fh_uh[11*ng+i] = fx_uh[11]*nx + fy_uh[11]*ny + rlamm3*(urE - rE);              
        fh_uh[12*ng+i] = fx_uh[12]*nx + fy_uh[12]*ny + rlamm4*(ur - r);
        fh_uh[13*ng+i] = fx_uh[13]*nx + fy_uh[13]*ny + rlamm4*(uru - ru);
        fh_uh[14*ng+i] = fx_uh[14]*nx + fy_uh[14]*ny + rlamm4*(urv - rv);
        fh_uh[15*ng+i] = fx_uh[15]*nx + fy_uh[15]*ny + rlamm4*(urE - rE) - rlam;                     

        fx_q[0] = 0.0;
        fx_q[1] = txx_rx;
        fx_q[2] = txy_rx;
        fx_q[3] = u*txx_rx  + v*txy_rx  + fc*Tx_rx;                              
        fx_q[4] = 0.0;
        fx_q[5] = txx_rux;
        fx_q[6] = txy_rux;
        fx_q[7] = u*txx_rux + v*txy_rux + fc*Tx_rux;        
        fx_q[8] = 0.0;
        fx_q[9] = txx_rvx;
        fx_q[10] = txy_rvx;
        fx_q[11] = u*txx_rvx + v*txy_rvx + fc*Tx_rvx;
        fx_q[12] = 0.0;
        fx_q[13] = txx_rEx;
        fx_q[14] = txy_rEx;
        fx_q[15] = u*txx_rEx + v*txy_rEx + fc*Tx_rEx;

        fx_q[16] = 0.0;
        fx_q[17] = txx_ry; 
        fx_q[18] = txy_ry;
        fx_q[19] = u*txx_ry  + v*txy_ry;
        fx_q[20] = 0.0;
        fx_q[21] = txx_ruy;
        fx_q[22] = txy_ruy;
        fx_q[23] = u*txx_ruy + v*txy_ruy;
        fx_q[24] = 0.0;
        fx_q[25] = txx_rvy; 
        fx_q[26] = txy_rvy; 
        fx_q[27] = u*txx_rvy + v*txy_rvy;
        fx_q[28] = 0.0;
        fx_q[29] = txx_rEy; 
        fx_q[30] = txy_rEy; 
        fx_q[31] = u*txx_rEy + v*txy_rEy;

        fy_q[0] = 0.0;
        fy_q[1] = txy_rx; 
        fy_q[2] = tyy_rx; 
        fy_q[3] = u*txy_rx  + v*tyy_rx;
        fy_q[4] = 0.0;
        fy_q[5] = txy_rux; 
        fy_q[6] = tyy_rux; 
        fy_q[7] = u*txy_rux + v*tyy_rux;
        fy_q[8] = 0.0;
        fy_q[9] = txy_rvx; 
        fy_q[10] = tyy_rvx; 
        fy_q[11] = u*txy_rvx + v*tyy_rvx;
        fy_q[12] = 0.0;
        fy_q[13] = txy_rEx;
        fy_q[14] = tyy_rEx;
        fy_q[15] = u*txy_rEx + v*tyy_rEx; 

        fy_q[16] = 0.0;
        fy_q[17] = txy_ry; 
        fy_q[18] = tyy_ry; 
        fy_q[19] = u*txy_ry  + v*tyy_ry  + fc*Ty_ry;
        fy_q[20] = 0.0;
        fy_q[21] = txy_ruy; 
        fy_q[22] = tyy_ruy; 
        fy_q[23] = u*txy_ruy + v*tyy_ruy + fc*Ty_ruy;
        fy_q[24] = 0.0;
        fy_q[25] = txy_rvy;
        fy_q[26] = tyy_rvy;
        fy_q[27] = u*txy_rvy + v*tyy_rvy + fc*Ty_rvy;
        fy_q[28] = 0.0;
        fy_q[29] = txy_rEy;
        fy_q[30] = tyy_rEy;
        fy_q[31] = u*txy_rEy + v*tyy_rEy + fc*Ty_rEy;                                    

        fh_u[16*ng+i] = fx_q[0]*nx + fy_q[0]*ny;
        fh_u[17*ng+i] = fx_q[1]*nx + fy_q[1]*ny;
        fh_u[18*ng+i] = fx_q[2]*nx + fy_q[2]*ny;
        fh_u[19*ng+i] = fx_q[3]*nx + fy_q[3]*ny;        
        fh_u[20*ng+i] = fx_q[4]*nx + fy_q[4]*ny;
        fh_u[21*ng+i] = fx_q[5]*nx + fy_q[5]*ny;
        fh_u[22*ng+i] = fx_q[6]*nx + fy_q[6]*ny;
        fh_u[23*ng+i] = fx_q[7]*nx + fy_q[7]*ny;        
        fh_u[24*ng+i] = fx_q[8]*nx + fy_q[8]*ny;
        fh_u[25*ng+i] = fx_q[9]*nx + fy_q[9]*ny;
        fh_u[26*ng+i] = fx_q[10]*nx + fy_q[10]*ny;
        fh_u[27*ng+i] = fx_q[11]*nx + fy_q[11]*ny;              
        fh_u[28*ng+i] = fx_q[12]*nx + fy_q[12]*ny;
        fh_u[29*ng+i] = fx_q[13]*nx + fy_q[13]*ny;
        fh_u[30*ng+i] = fx_q[14]*nx + fy_q[14]*ny;
        fh_u[31*ng+i] = fx_q[15]*nx + fy_q[15]*ny;

        fh_u[32*ng+i] = fx_q[16]*nx + fy_q[16]*ny;
        fh_u[33*ng+i] = fx_q[17]*nx + fy_q[17]*ny;
        fh_u[34*ng+i] = fx_q[18]*nx + fy_q[18]*ny;
        fh_u[35*ng+i] = fx_q[19]*nx + fy_q[19]*ny;        
        fh_u[36*ng+i] = fx_q[20]*nx + fy_q[20]*ny;
        fh_u[37*ng+i] = fx_q[21]*nx + fy_q[21]*ny;
        fh_u[38*ng+i] = fx_q[22]*nx + fy_q[22]*ny;
        fh_u[39*ng+i] = fx_q[23]*nx + fy_q[23]*ny;        
        fh_u[40*ng+i] = fx_q[24]*nx + fy_q[24]*ny;
        fh_u[41*ng+i] = fx_q[25]*nx + fy_q[25]*ny;
        fh_u[42*ng+i] = fx_q[26]*nx + fy_q[26]*ny;
        fh_u[43*ng+i] = fx_q[27]*nx + fy_q[27]*ny;              
        fh_u[44*ng+i] = fx_q[28]*nx + fy_q[28]*ny;
        fh_u[45*ng+i] = fx_q[29]*nx + fy_q[29]*ny;
        fh_u[46*ng+i] = fx_q[30]*nx + fy_q[30]*ny;
        fh_u[47*ng+i] = fx_q[31]*nx + fy_q[31]*ny;            
    }    
    
}

void getan_ns(double *An, double *Anm, double *uhg, double *param, double *nl, int absolute, int ng, int nch)
{
    
    double gam, epslm, gam1, signun, ic;
    double r, ru, rv, rE;    
    double r1, u, v, E, uu, vv, q, p, h;     
    double nx, ny;    
    double um1, um2, um3, um4, vm1, vm2, vm3, vm4;
    double Em1, Em2, Em3, Em4, hm1, hm2, hm3, hm4;    
    double r1m1, r1m2, r1m3, r1m4, qm1, qm2, qm3, qm4;
    double s1, s1m1, s1m2, s1m3, s1m4, s2, s2m1, s2m2, s2m3, s2m4;
    double pm1, pm2, pm3, pm4, c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4, un, unm1, unm2, unm3, unm4;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4;    
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4;            
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4;            
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4;            
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;            
    
    int    i;    
        
    gam  = param[0];    
    epslm= param[1];
    gam1 = gam - 1.0;
    
    for (i=0; i<ng; i++) { 
        r   = uhg[0*ng+i];
        ru  = uhg[1*ng+i];
        rv  = uhg[2*ng+i];
        rE  = uhg[3*ng+i];

        nx  = nl[0*ng+i];              
        ny  = nl[1*ng+i];   

        r1   = 1.0/r;
        u    = ru*r1;
        v    = rv*r1;
        E    = rE*r1;
        uu   = u*u;
        vv   = v*v;        
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

        Em1  =          rE*r1m1;
        Em2  =          rE*r1m2;
        Em3  =          rE*r1m3;
        Em4  =     r1 + rE*r1m4;

        qm1 = u*um1 + v*vm1;
        qm2 = u*um2 + v*vm2;
        qm3 = u*um3 + v*vm3;
        qm4 = u*um4 + v*vm4;

        pm1  = gam1*(   -  q  - r*qm1);
        pm2  = gam1*(         - r*qm2);
        pm3  = gam1*(         - r*qm3);
        pm4  = gam1*(1.0      - r*qm4);
        
        hm1  = Em1 + pm1*r1 + p*r1m1;
        hm2  = Em2 + pm2*r1 + p*r1m2;
        hm3  = Em3 + pm3*r1 + p*r1m3;
        hm4  = Em4 + pm4*r1 + p*r1m4;

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
            rlam1   = fabs(un+c);
            signun = (un+c < 0) ? -1.0 : 1.0;
            if (un+c==0) signun = 0;
            rlam1m1 = signun*(unm1+cm1);
            rlam1m2 = signun*(unm2+cm2);
            rlam1m3 = signun*(unm3+cm3);
            rlam1m4 = signun*(unm4+cm4);

            if (epslm>0) {
                rlam = 0.5*(rlam1*rlam1/(epslm*c)+epslm*c);
                rlamm1 = rlam1*rlam1m1/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm1);
                rlamm2 = rlam1*rlam1m2/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm2);
                rlamm3 = rlam1*rlam1m3/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm3);
                rlamm4 = rlam1*rlam1m4/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm4);
                ic = (rlam1 < epslm*c) ? 1.0 : 0.0;
                rlam1 = ic*rlam + (1-ic)*rlam1;
                rlam1m1 = ic*rlamm1 + (1-ic)*rlam1m1;
                rlam1m2 = ic*rlamm2 + (1-ic)*rlam1m2;
                rlam1m3 = ic*rlamm3 + (1-ic)*rlam1m3;
                rlam1m4 = ic*rlamm4 + (1-ic)*rlam1m4;
            }
        }
        else {
            rlam1   = un+c;
            rlam1m1 = unm1+cm1;
            rlam1m2 = unm2+cm2;
            rlam1m3 = unm3+cm3;
            rlam1m4 = unm4+cm4;
        }
        
        if (absolute==1) {            
            rlam2   = fabs(un-c);
            signun = (un-c < 0) ? -1.0 : 1.0;
            if (un-c==0) signun = 0;
            rlam2m1 = signun*(unm1-cm1);
            rlam2m2 = signun*(unm2-cm2);
            rlam2m3 = signun*(unm3-cm3);
            rlam2m4 = signun*(unm4-cm4);

            if (epslm>0) {
                rlam = 0.5*(rlam2*rlam2/(epslm*c)+epslm*c);
                rlamm1 = rlam2*rlam2m1/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm1);
                rlamm2 = rlam2*rlam2m2/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm2);
                rlamm3 = rlam2*rlam2m3/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm3);
                rlamm4 = rlam2*rlam2m4/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm4);
                ic = (rlam2 < epslm*c) ? 1.0 : 0.0;
                rlam2 = ic*rlam + (1-ic)*rlam2;
                rlam2m1 = ic*rlamm1 + (1-ic)*rlam2m1;
                rlam2m2 = ic*rlamm2 + (1-ic)*rlam2m2;
                rlam2m3 = ic*rlamm3 + (1-ic)*rlam2m3;
                rlam2m4 = ic*rlamm4 + (1-ic)*rlam2m4;
            }
        }
        else {
            rlam2   = un-c;
            rlam2m1 = unm1-cm1;
            rlam2m2 = unm2-cm2;
            rlam2m3 = unm3-cm3;
            rlam2m4 = unm4-cm4;
        }

        if (absolute==1) {
            rlam3   = fabs(un);
            signun = (un < 0) ? -1.0 : 1.0;
            if (un==0) signun = 0;
            rlam3m1 = signun*unm1;
            rlam3m2 = signun*unm2;
            rlam3m3 = signun*unm3;
            rlam3m4 = signun*unm4;

            if (epslm>0) {
                rlam = 0.5*(rlam3*rlam3/(epslm*c)+epslm*c);
                rlamm1 = rlam3*rlam3m1/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm1);
                rlamm2 = rlam3*rlam3m2/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm2);
                rlamm3 = rlam3*rlam3m3/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm3);
                rlamm4 = rlam3*rlam3m4/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm4);
                ic = (rlam3 < epslm*c) ? 1.0 : 0.0;
                rlam3 = ic*rlam + (1-ic)*rlam3;
                rlam3m1 = ic*rlamm1 + (1-ic)*rlam3m1;
                rlam3m2 = ic*rlamm2 + (1-ic)*rlam3m2;
                rlam3m3 = ic*rlamm3 + (1-ic)*rlam3m3;
                rlam3m4 = ic*rlamm4 + (1-ic)*rlam3m4;
            }
        }
        else {
            rlam3   = un;
            rlam3m1 = unm1;
            rlam3m2 = unm2;
            rlam3m3 = unm3;
            rlam3m4 = unm4;
        }

        s1      = 0.5*(rlam1   + rlam2);
        s1m1    = 0.5*(rlam1m1 + rlam2m1);
        s1m2    = 0.5*(rlam1m2 + rlam2m2);
        s1m3    = 0.5*(rlam1m3 + rlam2m3);
        s1m4    = 0.5*(rlam1m4 + rlam2m4);

        s2      = 0.5*(rlam1   - rlam2);
        s2m1    = 0.5*(rlam1m1 - rlam2m1);
        s2m2    = 0.5*(rlam1m2 - rlam2m2);
        s2m3    = 0.5*(rlam1m3 - rlam2m3);
        s2m4    = 0.5*(rlam1m4 - rlam2m4);

        cc1   = gam1*(s1-rlam3)*q/c2-(s2*un/c);
        cc1m1 = gam1*((s1m1-rlam3m1)*q/c2 + (s1-rlam3)*qm1/c2 - (s1-rlam3)*q*c2m1/(c2*c2)) - s2m1*un/c - s2*unm1/c + s2*un*cm1/(c*c);
        cc1m2 = gam1*((s1m2-rlam3m2)*q/c2 + (s1-rlam3)*qm2/c2 - (s1-rlam3)*q*c2m2/(c2*c2)) - s2m2*un/c - s2*unm2/c + s2*un*cm2/(c*c);
        cc1m3 = gam1*((s1m3-rlam3m3)*q/c2 + (s1-rlam3)*qm3/c2 - (s1-rlam3)*q*c2m3/(c2*c2)) - s2m3*un/c - s2*unm3/c + s2*un*cm3/(c*c);
        cc1m4 = gam1*((s1m4-rlam3m4)*q/c2 + (s1-rlam3)*qm4/c2 - (s1-rlam3)*q*c2m4/(c2*c2)) - s2m4*un/c - s2*unm4/c + s2*un*cm4/(c*c);

        cc2   = gam1*s2*q/c-(s1-rlam3)*un;
        cc2m1 = gam1*(s2m1*q/c + s2*qm1/c - s2*q*cm1/(c*c)) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
        cc2m2 = gam1*(s2m2*q/c + s2*qm2/c - s2*q*cm2/(c*c)) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
        cc2m3 = gam1*(s2m3*q/c + s2*qm3/c - s2*q*cm3/(c*c)) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
        cc2m4 = gam1*(s2m4*q/c + s2*qm4/c - s2*q*cm4/(c*c)) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;
        
        An[0*ng+i] = rlam3+cc1;
        An[1*ng+i] = cc1*u+cc2*nx;
        An[2*ng+i] = cc1*v+cc2*ny;
        An[3*ng+i] = cc1*h+cc2*un;        
        
        Anm[0*ng+i] = rlam3m1+cc1m1;
        Anm[1*ng+i] = cc1m1*u+cc1*um1+cc2m1*nx;
        Anm[2*ng+i] = cc1m1*v+cc1*vm1+cc2m1*ny;
        Anm[3*ng+i] = cc1m1*h+cc1*hm1+cc2m1*un+cc2*unm1;        
        
        Anm[16*ng+i] = rlam3m2+cc1m2;
        Anm[17*ng+i] = cc1m2*u+cc1*um2+cc2m2*nx;
        Anm[18*ng+i] = cc1m2*v+cc1*vm2+cc2m2*ny;
        Anm[19*ng+i] = cc1m2*h+cc1*hm2+cc2m2*un+cc2*unm2;        
        
        Anm[32*ng+i] = rlam3m3+cc1m3;
        Anm[33*ng+i] = cc1m3*u+cc1*um3+cc2m3*nx;
        Anm[34*ng+i] = cc1m3*v+cc1*vm3+cc2m3*ny;
        Anm[35*ng+i] = cc1m3*h+cc1*hm3+cc2m3*un+cc2*unm3;        
        
        Anm[48*ng+i] = rlam3m4+cc1m4;
        Anm[49*ng+i] = cc1m4*u+cc1*um4+cc2m4*nx;
        Anm[50*ng+i] = cc1m4*v+cc1*vm4+cc2m4*ny;
        Anm[51*ng+i] = cc1m4*h+cc1*hm4+cc2m4*un+cc2*unm4;                                                

        cc1   = -gam1*(s1-rlam3)*u/c2+(s2*nx/c);
        cc1m1 = -gam1*((s1m1-rlam3m1)*u/c2 + (s1-rlam3)*um1/c2 - (s1-rlam3)*u*c2m1/(c2*c2)) + s2m1*nx/c - s2*nx*cm1/(c*c);
        cc1m2 = -gam1*((s1m2-rlam3m2)*u/c2 + (s1-rlam3)*um2/c2 - (s1-rlam3)*u*c2m2/(c2*c2)) + s2m2*nx/c - s2*nx*cm2/(c*c);
        cc1m3 = -gam1*((s1m3-rlam3m3)*u/c2 + (s1-rlam3)*um3/c2 - (s1-rlam3)*u*c2m3/(c2*c2)) + s2m3*nx/c - s2*nx*cm3/(c*c);
        cc1m4 = -gam1*((s1m4-rlam3m4)*u/c2 + (s1-rlam3)*um4/c2 - (s1-rlam3)*u*c2m4/(c2*c2)) + s2m4*nx/c - s2*nx*cm4/(c*c);

        cc2   = -gam1*s2*u/c + (s1-rlam3)*nx;
        cc2m1 = -gam1*(s2m1*u/c + s2*um1/c - s2*u*cm1/(c*c)) + (s1m1-rlam3m1)*nx;
        cc2m2 = -gam1*(s2m2*u/c + s2*um2/c - s2*u*cm2/(c*c)) + (s1m2-rlam3m2)*nx;
        cc2m3 = -gam1*(s2m3*u/c + s2*um3/c - s2*u*cm3/(c*c)) + (s1m3-rlam3m3)*nx;
        cc2m4 = -gam1*(s2m4*u/c + s2*um4/c - s2*u*cm4/(c*c)) + (s1m4-rlam3m4)*nx;

        An[4*ng+i] = cc1;
        An[5*ng+i] = rlam3+cc1*u+cc2*nx;
        An[6*ng+i] = cc1*v+cc2*ny;
        An[7*ng+i] = cc1*h+cc2*un;                
        
        Anm[4*ng+i] = cc1m1;
        Anm[5*ng+i] = rlam3m1+cc1m1*u+cc1*um1+cc2m1*nx;
        Anm[6*ng+i] = cc1m1*v+cc1*vm1+cc2m1*ny;
        Anm[7*ng+i] = cc1m1*h+cc1*hm1+cc2m1*un+cc2*unm1;                
        
        Anm[20*ng+i] = cc1m2;
        Anm[21*ng+i] = rlam3m2+cc1m2*u+cc1*um2+cc2m2*nx;
        Anm[22*ng+i] = cc1m2*v+cc1*vm2+cc2m2*ny;
        Anm[23*ng+i] = cc1m2*h+cc1*hm2+cc2m2*un+cc2*unm2;                        
        
        Anm[36*ng+i] = cc1m3;
        Anm[37*ng+i] = rlam3m3+cc1m3*u+cc1*um3+cc2m3*nx;
        Anm[38*ng+i] = cc1m3*v+cc1*vm3+cc2m3*ny;
        Anm[39*ng+i] = cc1m3*h+cc1*hm3+cc2m3*un+cc2*unm3;                
        
        Anm[52*ng+i] = cc1m4;
        Anm[53*ng+i] = rlam3m4+cc1m4*u+cc1*um4+cc2m4*nx;
        Anm[54*ng+i] = cc1m4*v+cc1*vm4+cc2m4*ny;
        Anm[55*ng+i] = cc1m4*h+cc1*hm4+cc2m4*un+cc2*unm4;                                                        
        
        cc1   = -gam1*(s1-rlam3)*v/c2+(s2*ny/c);
        cc1m1 = -gam1*((s1m1-rlam3m1)*v/c2 + (s1-rlam3)*vm1/c2 - (s1-rlam3)*v*c2m1/(c2*c2)) + s2m1*ny/c - s2*ny*cm1/(c*c);
        cc1m2 = -gam1*((s1m2-rlam3m2)*v/c2 + (s1-rlam3)*vm2/c2 - (s1-rlam3)*v*c2m2/(c2*c2)) + s2m2*ny/c - s2*ny*cm2/(c*c);
        cc1m3 = -gam1*((s1m3-rlam3m3)*v/c2 + (s1-rlam3)*vm3/c2 - (s1-rlam3)*v*c2m3/(c2*c2)) + s2m3*ny/c - s2*ny*cm3/(c*c);
        cc1m4 = -gam1*((s1m4-rlam3m4)*v/c2 + (s1-rlam3)*vm4/c2 - (s1-rlam3)*v*c2m4/(c2*c2)) + s2m4*ny/c - s2*ny*cm4/(c*c);

        cc2   = -gam1*s2*v/c+(s1-rlam3)*ny;
        cc2m1 = -gam1*(s2m1*v/c + s2*vm1/c - s2*v*cm1/(c*c)) + (s1m1-rlam3m1)*ny;
        cc2m2 = -gam1*(s2m2*v/c + s2*vm2/c - s2*v*cm2/(c*c)) + (s1m2-rlam3m2)*ny;
        cc2m3 = -gam1*(s2m3*v/c + s2*vm3/c - s2*v*cm3/(c*c)) + (s1m3-rlam3m3)*ny;
        cc2m4 = -gam1*(s2m4*v/c + s2*vm4/c - s2*v*cm4/(c*c)) + (s1m4-rlam3m4)*ny;

        An[8*ng+i] = cc1;
        An[9*ng+i] = cc1*u+cc2*nx;
        An[10*ng+i] = rlam3+cc1*v+cc2*ny;
        An[11*ng+i] = cc1*h+cc2*un;                
        
        Anm[8*ng+i] = cc1m1;
        Anm[9*ng+i] = cc1m1*u+cc1*um1+cc2m1*nx;
        Anm[10*ng+i] = rlam3m1+cc1m1*v+cc1*vm1+cc2m1*ny;
        Anm[11*ng+i] = cc1m1*h+cc1*hm1+cc2m1*un+cc2*unm1;                
        
        Anm[24*ng+i] = cc1m2;
        Anm[25*ng+i] = cc1m2*u+cc1*um2+cc2m2*nx;
        Anm[26*ng+i] = rlam3m2+cc1m2*v+cc1*vm2+cc2m2*ny;
        Anm[27*ng+i] = cc1m2*h+cc1*hm2+cc2m2*un+cc2*unm2;                        
        
        Anm[40*ng+i] = cc1m3;
        Anm[41*ng+i] = cc1m3*u+cc1*um3+cc2m3*nx;
        Anm[42*ng+i] = rlam3m3+cc1m3*v+cc1*vm3+cc2m3*ny;
        Anm[43*ng+i] = cc1m3*h+cc1*hm3+cc2m3*un+cc2*unm3;                
        
        Anm[56*ng+i] = cc1m4;
        Anm[57*ng+i] = cc1m4*u+cc1*um4+cc2m4*nx;
        Anm[58*ng+i] = rlam3m4+cc1m4*v+cc1*vm4+cc2m4*ny;
        Anm[59*ng+i] = cc1m4*h+cc1*hm4+cc2m4*un+cc2*unm4;                                                        
        
        cc1   = gam1*(s1-rlam3)/c2;
        cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/(c2*c2));
        cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/(c2*c2));
        cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/(c2*c2));
        cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/(c2*c2));

        cc2   = gam1*s2/c;
        cc2m1 = gam1*(s2m1/c - s2*cm1/(c*c));
        cc2m2 = gam1*(s2m2/c - s2*cm2/(c*c));
        cc2m3 = gam1*(s2m3/c - s2*cm3/(c*c));
        cc2m4 = gam1*(s2m4/c - s2*cm4/(c*c));

        An[12*ng+i] = cc1;
        An[13*ng+i] = cc1*u+cc2*nx;
        An[14*ng+i] = cc1*v+cc2*ny;
        An[15*ng+i] = rlam3+cc1*h+cc2*un;                
        
        Anm[12*ng+i] = cc1m1;
        Anm[13*ng+i] = cc1m1*u+cc1*um1+cc2m1*nx;
        Anm[14*ng+i] = cc1m1*v+cc1*vm1+cc2m1*ny;
        Anm[15*ng+i] = rlam3m1+cc1m1*h+cc1*hm1+cc2m1*un+cc2*unm1;                
        
        Anm[28*ng+i] = cc1m2;
        Anm[29*ng+i] = cc1m2*u+cc1*um2+cc2m2*nx;
        Anm[30*ng+i] = cc1m2*v+cc1*vm2+cc2m2*ny;
        Anm[31*ng+i] = rlam3m2+cc1m2*h+cc1*hm2+cc2m2*un+cc2*unm2;                        
        
        Anm[44*ng+i] = cc1m3;
        Anm[45*ng+i] = cc1m3*u+cc1*um3+cc2m3*nx;
        Anm[46*ng+i] = cc1m3*v+cc1*vm3+cc2m3*ny;
        Anm[47*ng+i] = rlam3m3+cc1m3*h+cc1*hm3+cc2m3*un+cc2*unm3;                
        
        Anm[60*ng+i] = cc1m4;
        Anm[61*ng+i] = cc1m4*u+cc1*um4+cc2m4*nx;
        Anm[62*ng+i] = cc1m4*v+cc1*vm4+cc2m4*ny;
        Anm[63*ng+i] = rlam3m4+cc1m4*h+cc1*hm4+cc2m4*un+cc2*unm4;                                                                
    }    
    
}

void fbou_ns(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd)
{        
    
    double gam, epslm, gam1, Minf, M2;    
    double nx, ny, tm;        
    int    i, j, k, m, n, nm, nk;
                   
    gam  = param[0];
    epslm= param[1];                           
    Minf = param[4];    
    gam1 = gam-1.0;    
    M2   = Minf*Minf;            

    
    if (ib==1) {

        double *an  = (double *) malloc( ng*nch*nch * sizeof(double) ); 
        double *An  = (double *) malloc( ng*nch*nch * sizeof(double) );     
        double *anm = (double *) malloc( ng*nch*nch*nch * sizeof(double) );     
        double *Anm = (double *) malloc( ng*nch*nch*nch * sizeof(double) );                        
        double *u   = (double *) malloc( nch * sizeof(double) ); 
        double *uh  = (double *) malloc( nch * sizeof(double) ); 
        double *rlamm  = (double *) malloc( nch*nch * sizeof(double) ); 
        double *rlamp  = (double *) malloc( nch*nch * sizeof(double) ); 
        double *rlammm  = (double *) malloc( nch*nch*nch * sizeof(double) ); 
        double *rlampm  = (double *) malloc( nch*nch*nch * sizeof(double) ); 
            
        getan_ns(an,anm,uhg,param,nl,0,ng,nch);
        getan_ns(An,Anm,uhg,param,nl,1,ng,nch);                        
        
        for (j=0; j<nc; j++)
            for (k=0; k<nch; k++) 
                for (i=0; i<ng; i++) {
                    nm = j*nch*ng+k*ng+i;
                    fh_u[nm] =0;
                }
        
        for (i=0; i<ng; i++) { 
            u[0] = udg[0*ng+i];
            u[1] = udg[1*ng+i];
            u[2] = udg[2*ng+i];
            u[3] = udg[3*ng+i];

            uh[0]  = uhg[0*ng+i];
            uh[1]  = uhg[1*ng+i];
            uh[2]  = uhg[2*ng+i];
            uh[3]  = uhg[3*ng+i];

            for (j=0; j<nch; j++)
                for (k=0; k<nch; k++) {                        
                    nm = j*nch*ng+k*ng+i;
                    tm = an[nm] + An[nm];
                    fh_u[nm]  = tm;

                    nk = j*nch+k;
                    rlamp[nk] = tm;
                    rlamm[nk] = an[nm] - An[nm];      

                    for (n=0; n<nch; n++) {
                        nm = n*nch*nch*ng+j*nch*ng+k*ng+i;
                        nk = n*nch*nch+j*nch+k;
                        rlampm[nk] =  anm[nm] + Anm[nm];
                        rlammm[nk] =  anm[nm] - Anm[nm];
                    }
                }

            for (j=0; j<nch; j++)
                fh[j*ng+i] = rlamp[0+j]*(u[0] - uh[0]) + rlamp[4+j]*(u[1] - uh[1]) +     
                             rlamp[8+j]*(u[2] - uh[2]) + rlamp[12+j]*(u[3] - uh[3]) -     
                             rlamm[0+j]*(ui[0] - uh[0]) - rlamm[4+j]*(ui[1] - uh[1]) -     
                             rlamm[8+j]*(ui[2] - uh[2]) - rlamm[12+j]*(ui[3] - uh[3]);                                                                                  

            //fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
            
            for (j=0; j<nch; j++)
                for (k=0; k<nch; k++) {
                    nm = j*nch*ng+k*ng+i;
                    nk = j*nch*nch+k;
                    fh_uh[nm] = rlampm[0+nk]*(u[0] - uh[0]) + rlampm[4+nk]*(u[1] - uh[1]) +     
                                rlampm[8+nk]*(u[2] - uh[2]) + rlampm[12+nk]*(u[3] - uh[3]) -     
                                rlammm[0+nk]*(ui[0] - uh[0]) - rlammm[4+nk]*(ui[1] - uh[1]) -     
                                rlammm[8+nk]*(ui[2] - uh[2]) - rlammm[12+nk]*(ui[3] - uh[3]) -                                                                  
                                2*An[nm];                                                                                  
                }                
        }            
        
        free(an); free(An); free(anm); free(Anm);
        free(uh); free(u); free(rlamm); free(rlamp);
        free(rlammm); free(rlampm);   
    }                                     
    else if(ib==2) { 

        fhat_ns(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd);                          

        for (i=0; i<ng; i++) {

            fh[0*ng+i] =  udg[0*ng+i]-uhg[0*ng+i];
            fh[1*ng+i] = -uhg[1*ng+i];
            fh[2*ng+i] = -uhg[2*ng+i];

            fh_uh[0*ng+i] = -1.0;
            fh_uh[1*ng+i] = 0.0;
            fh_uh[2*ng+i] = 0.0; 
            fh_uh[4*ng+i] = 0.0; 
            fh_uh[5*ng+i] = -1.0;
            fh_uh[6*ng+i] = 0.0;
            fh_uh[8*ng+i] = 0.0;
            fh_uh[9*ng+i] = 0.0;
            fh_uh[10*ng+i] = -1.0;
            fh_uh[12*ng+i] = 0.0;
            fh_uh[13*ng+i] = 0.0;
            fh_uh[14*ng+i] = 0.0;

            fh_u[0*ng+i] = 1.0;
            fh_u[1*ng+i] = 0.0;
            fh_u[2*ng+i] = 0.0;
            fh_u[4*ng+i] = 0.0;
            fh_u[5*ng+i] = 0.0;
            fh_u[6*ng+i] = 0.0;                    
            fh_u[8*ng+i] = 0.0;
            fh_u[9*ng+i] = 0.0;
            fh_u[10*ng+i] = 0.0;
            fh_u[12*ng+i] = 0.0;
            fh_u[13*ng+i] = 0.0;
            fh_u[14*ng+i] = 0.0;                                                           

            fh_u[16*ng+i] = 0.0;
            fh_u[17*ng+i] = 0.0;
            fh_u[18*ng+i] = 0.0;                    
            fh_u[20*ng+i] = 0.0;
            fh_u[21*ng+i] = 0.0;
            fh_u[22*ng+i] = 0.0;                    
            fh_u[24*ng+i] = 0.0;
            fh_u[25*ng+i] = 0.0;
            fh_u[26*ng+i] = 0.0;                    
            fh_u[28*ng+i] = 0.0;
            fh_u[29*ng+i] = 0.0;
            fh_u[30*ng+i] = 0.0;                    

            fh_u[32*ng+i] = 0.0;
            fh_u[33*ng+i] = 0.0;
            fh_u[34*ng+i] = 0.0;               
            fh_u[36*ng+i] = 0.0;
            fh_u[37*ng+i] = 0.0;
            fh_u[38*ng+i] = 0.0;                    
            fh_u[40*ng+i] = 0.0;
            fh_u[41*ng+i] = 0.0;
            fh_u[42*ng+i] = 0.0;                    
            fh_u[44*ng+i] = 0.0;
            fh_u[45*ng+i] = 0.0;
            fh_u[46*ng+i] = 0.0;
        }            
    }
    else if (ib==3) { 

        fhat_ns(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd);

        for (i=0; i<ng; i++) { 
            fh[1*ng+i] = uhg[1*ng+i];
            fh[2*ng+i] = uhg[2*ng+i];
            fh[3*ng+i] = gam*gam1*M2*uhg[3*ng+i]/uhg[0*ng+i] - ui[3];

            fh_uh[1*ng+i] = 0.0;
            fh_uh[2*ng+i] = 0.0;
            fh_uh[3*ng+i] = -gam*gam1*M2*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
            fh_uh[5*ng+i] = 1.0;
            fh_uh[6*ng+i] = 0.0;
            fh_uh[7*ng+i] = 0.0;
            fh_uh[9*ng+i] = 0.0;
            fh_uh[10*ng+i] = 1.0;
            fh_uh[11*ng+i] = 0.0;
            fh_uh[13*ng+i] = 0.0;
            fh_uh[14*ng+i] = 0.0;
            fh_uh[15*ng+i] = gam*gam1*M2/uhg[0*ng+i];

            fh_u[1*ng+i] = 0.0;
            fh_u[2*ng+i] = 0.0;
            fh_u[3*ng+i] = 0.0;
            fh_u[5*ng+i] = 0.0;
            fh_u[6*ng+i] = 0.0;                    
            fh_u[7*ng+i] = 0.0;                    
            fh_u[9*ng+i] = 0.0;
            fh_u[10*ng+i] = 0.0;
            fh_u[11*ng+i] = 0.0;                    
            fh_u[13*ng+i] = 0.0;
            fh_u[14*ng+i] = 0.0;                                        
            fh_u[15*ng+i] = 0.0;                                        

            fh_u[17*ng+i] = 0.0;
            fh_u[18*ng+i] = 0.0;                    
            fh_u[19*ng+i] = 0.0;                    
            fh_u[21*ng+i] = 0.0;
            fh_u[22*ng+i] = 0.0;                    
            fh_u[23*ng+i] = 0.0;                    
            fh_u[25*ng+i] = 0.0;
            fh_u[26*ng+i] = 0.0;                    
            fh_u[27*ng+i] = 0.0;                    
            fh_u[29*ng+i] = 0.0;
            fh_u[30*ng+i] = 0.0;                    
            fh_u[31*ng+i] = 0.0;                    

            fh_u[33*ng+i] = 0.0;
            fh_u[34*ng+i] = 0.0;                    
            fh_u[35*ng+i] = 0.0;                    
            fh_u[37*ng+i] = 0.0;
            fh_u[38*ng+i] = 0.0;                    
            fh_u[39*ng+i] = 0.0;                    
            fh_u[41*ng+i] = 0.0;
            fh_u[42*ng+i] = 0.0;                    
            fh_u[43*ng+i] = 0.0;                    
            fh_u[45*ng+i] = 0.0;
            fh_u[46*ng+i] = 0.0;
            fh_u[47*ng+i] = 0.0;

        }                         
    }
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}

