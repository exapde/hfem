void getan_ns2d(double *An, double *Anm, double *uh, double *nl, double *param, int fabsflag, int ng)
{    
    double gam, epslm, gam1, signun, ic;
    double nx, ny;    
    double r, ru, rv, rE;          
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double E, Em1, Em2, Em3, Em4;
    double af, afm1, afm2, afm3, afm4;
    double p, pm1, pm2, pm3, pm4;
    double h, hm1, hm2, hm3, hm4;            
    double s1, s1m1, s1m2, s1m3, s1m4;
    double s2, s2m1, s2m2, s2m3, s2m4;    
    double c2, c2m1, c2m2, c2m3, c2m4;
    double c, cm1, cm2, cm3, cm4;
    double un, unm1, unm2, unm3, unm4;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4;    
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4;           
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4;            
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4;            
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4;                 
    
    gam   = param[0];
    epslm = param[1];
    gam1 = gam - 1.0;
    
    int    i;                          
    for (i=0; i<ng; i++) { 
        r   = uhg[0*ng+i];
        ru  = uhg[1*ng+i];
        rv  = uhg[2*ng+i];
        rE  = uhg[3*ng+i];

        nx  = nl[0*ng+i];              
        ny  = nl[1*ng+i];   
        
        r1   = 1.0/r;
        r1m1 = -1.0/(r**2);
        r1m2 = 0.0;
        r1m3 = 0.0;
        r1m4 = 0.0;

        uv   = ru*r1;
        uvm1 =          ru*r1m1;
        uvm2 =     r1 + ru*r1m2;
        uvm3 =          ru*r1m3;
        uvm4 =          ru*r1m4;

        vv   = rv*r1;
        vvm1 =          rv*r1m1;
        vvm2 =          rv*r1m2;
        vvm3 =     r1 + rv*r1m3;
        vvm4 =          rv*r1m4;

        E    = rE*r1;
        Em1  =          rE*r1m1;
        Em2  =          rE*r1m2;
        Em3  =          rE*r1m3;
        Em4  =     r1 + rE*r1m4;

        af   = 0.5*(uv*uv+vv*vv);
        afm1 = uv*uvm1 + vv*vvm1;
        afm2 = uv*uvm2 + vv*vvm2;
        afm3 = uv*uvm3 + vv*vvm3;
        afm4 = uv*uvm4 + vv*vvm4;

        p    = gam1*(rE -r*af);
        pm1  = gam1*(   -   af - r*afm1);
        pm2  = gam1*(          - r*afm2);
        pm3  = gam1*(          - r*afm3);
        pm4  = gam1*(1.0       - r*afm4);

        h    = E   + p*r1;
        hm1  = Em1 + pm1*r1 + p*r1m1;
        hm2  = Em2 + pm2*r1 + p*r1m2;
        hm3  = Em3 + pm3*r1 + p*r1m3;
        hm4  = Em4 + pm4*r1 + p*r1m4;

        c2   = gam* p*r1;
        c2m1 = gam*(pm1*r1 + p*r1m1);
        c2m2 = gam*(pm2*r1 + p*r1m2);
        c2m3 = gam*(pm3*r1 + p*r1m3);
        c2m4 = gam*(pm4*r1 + p*r1m4);

        c    = sqrt(c2);
        cm1  = 0.5*c2m1/c;
        cm2  = 0.5*c2m2/c;
        cm3  = 0.5*c2m3/c;
        cm4  = 0.5*c2m4/c;

        un   = uv*nx   + vv*ny;
        unm1 = uvm1*nx + vvm1*ny;
        unm2 = uvm2*nx + vvm2*ny;
        unm3 = uvm3*nx + vvm3*ny;
        unm4 = uvm4*nx + vvm4*ny;

        if fabsflag {
            rlam1   = fabs(un+c);
            signun = (un+c < 0) ? -1.0 : 1.0;
            rlam1m1 = signun*(unm1+cm1);
            rlam1m2 = signun*(unm2+cm2);
            rlam1m3 = signun*(unm3+cm3);
            rlam1m4 = signun*(unm4+cm4);

            if epslm>0 {
                rlam = 0.5*(rlam1*rlam1/(epslm*c)+epslm*c);
                rlamm1 = rlam1*rlam1m1/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c)**2)*(epslm*cm1);
                rlamm2 = rlam1*rlam1m2/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c)**2)*(epslm*cm2);
                rlamm3 = rlam1*rlam1m3/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c)**2)*(epslm*cm3);
                rlamm4 = rlam1*rlam1m4/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c)**2)*(epslm*cm4);
                ic = rlam1 < epslm*c;
                rlam1 = ic*rlam + (1-ic)*rlam1;
                rlam1m1 = ic*rlamm1 + (1-ic)*rlam1m1;
                rlam1m2 = ic*rlamm2 + (1-ic)*rlam1m2;
                rlam1m3 = ic*rlamm3 + (1-ic)*rlam1m3;
                rlam1m4 = ic*rlamm4 + (1-ic)*rlam1m4;
            }

            rlam2   = fabs(un-c);
            signun = (un-c < 0) ? -1.0 : 1.0;
            rlam2m1 = signun*(unm1-cm1);
            rlam2m2 = signun*(unm2-cm2);
            rlam2m3 = signun*(unm3-cm3);
            rlam2m4 = signun*(unm4-cm4);

            if epslm>0 {
                rlam = 0.5*(rlam2*rlam2/(epslm*c)+epslm*c);
                rlamm1 = rlam2*rlam2m1/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c)**2)*(epslm*cm1);
                rlamm2 = rlam2*rlam2m2/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c)**2)*(epslm*cm2);
                rlamm3 = rlam2*rlam2m3/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c)**2)*(epslm*cm3);
                rlamm4 = rlam2*rlam2m4/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c)**2)*(epslm*cm4);
                ic = rlam2 < epslm*c;
                rlam2 = ic*rlam + (1-ic)*rlam2;
                rlam2m1 = ic*rlamm1 + (1-ic)*rlam2m1;
                rlam2m2 = ic*rlamm2 + (1-ic)*rlam2m2;
                rlam2m3 = ic*rlamm3 + (1-ic)*rlam2m3;
                rlam2m4 = ic*rlamm4 + (1-ic)*rlam2m4;
            }

            rlam3   = fabs(un);
            signun = (un < 0) ? -1.0 : 1.0;
            rlam3m1 = signun*unm1;
            rlam3m2 = signun*unm2;
            rlam3m3 = signun*unm3;
            rlam3m4 = signun*unm4;

            if epslm>0 {
                rlam = 0.5*(rlam3*rlam3/(epslm*c)+epslm*c);
                rlamm1 = rlam3*rlam3m1/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c)**2)*(epslm*cm1);
                rlamm2 = rlam3*rlam3m2/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c)**2)*(epslm*cm2);
                rlamm3 = rlam3*rlam3m3/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c)**2)*(epslm*cm3);
                rlamm4 = rlam3*rlam3m4/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c)**2)*(epslm*cm4);
                ic = rlam3 < epslm*c;
                rlam3 = ic*rlam + (1-ic)*rlam3;
                rlam3m1 = ic*rlamm1 + (1-ic)*rlam3m1;
                rlam3m2 = ic*rlamm2 + (1-ic)*rlam3m2;
                rlam3m3 = ic*rlamm3 + (1-ic)*rlam3m3;
                rlam3m4 = ic*rlamm4 + (1-ic)*rlam3m4;
            }
        }
        else {
            rlam1   = un+c;
            rlam1m1 = unm1+cm1;
            rlam1m2 = unm2+cm2;
            rlam1m3 = unm3+cm3;
            rlam1m4 = unm4+cm4;

            rlam2   = un-c;
            rlam2m1 = unm1-cm1;
            rlam2m2 = unm2-cm2;
            rlam2m3 = unm3-cm3;
            rlam2m4 = unm4-cm4;

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

        cc1   = gam1*(s1-rlam3)*af/c2-(s2*un/c);
        cc1m1 = gam1*((s1m1-rlam3m1)*af/c2 + (s1-rlam3)*afm1/c2 - (s1-rlam3)*af*c2m1/c2**2) - s2m1*un/c - s2*unm1/c + s2*un*cm1/c**2;
        cc1m2 = gam1*((s1m2-rlam3m2)*af/c2 + (s1-rlam3)*afm2/c2 - (s1-rlam3)*af*c2m2/c2**2) - s2m2*un/c - s2*unm2/c + s2*un*cm2/c**2;
        cc1m3 = gam1*((s1m3-rlam3m3)*af/c2 + (s1-rlam3)*afm3/c2 - (s1-rlam3)*af*c2m3/c2**2) - s2m3*un/c - s2*unm3/c + s2*un*cm3/c**2;
        cc1m4 = gam1*((s1m4-rlam3m4)*af/c2 + (s1-rlam3)*afm4/c2 - (s1-rlam3)*af*c2m4/c2**2) - s2m4*un/c - s2*unm4/c + s2*un*cm4/c**2;

        cc2   = gam1*s2*af/c-(s1-rlam3)*un;
        cc2m1 = gam1*(s2m1*af/c + s2*afm1/c - s2*af*cm1/c**2) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
        cc2m2 = gam1*(s2m2*af/c + s2*afm2/c - s2*af*cm2/c**2) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
        cc2m3 = gam1*(s2m3*af/c + s2*afm3/c - s2*af*cm3/c**2) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
        cc2m4 = gam1*(s2m4*af/c + s2*afm4/c - s2*af*cm4/c**2) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;

        An[0*ng+i]  = rlam3+cc1;
        An[1*ng+i]  = cc1*uv+cc2*nx;
        An[2*ng+i]  = cc1*vv+cc2*ny;
        An[3*ng+i]  = cc1*h+cc2*un;
        Anm[0*ng+i] = rlam3m1+cc1m1;
        Anm[1*ng+i] = cc1m1*uv+cc1*uvm1+cc2m1*nx;
        Anm[2*ng+i] = cc1m1*vv+cc1*vvm1+cc2m1*ny;
        Anm[3*ng+i] = cc1m1*h+cc1*hm1+cc2m1*un+cc2*unm1;
        Anm[16*ng+i] = rlam3m2+cc1m2;
        Anm[17*ng+i] = cc1m2*uv+cc1*uvm2+cc2m2*nx;
        Anm[18*ng+i] = cc1m2*vv+cc1*vvm2+cc2m2*ny;
        Anm[19*ng+i] = cc1m2*h+cc1*hm2+cc2m2*un+cc2*unm2;
        Anm[32*ng+i] = rlam3m3+cc1m3;
        Anm[33*ng+i] = cc1m3*uv+cc1*uvm3+cc2m3*nx;
        Anm[34*ng+i] = cc1m3*vv+cc1*vvm3+cc2m3*ny;
        Anm[35*ng+i] = cc1m3*h+cc1*hm3+cc2m3*un+cc2*unm3;
        Anm[48*ng+i] = rlam3m4+cc1m4;
        Anm[49*ng+i] = cc1m4*uv+cc1*uvm4+cc2m4*nx;
        Anm[50*ng+i] = cc1m4*vv+cc1*vvm4+cc2m4*ny;
        Anm[51*ng+i] = cc1m4*h+cc1*hm4+cc2m4*un+cc2*unm4;

        cc1   = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
        cc1m1 = -gam1*((s1m1-rlam3m1)*uv/c2 + (s1-rlam3)*uvm1/c2 - (s1-rlam3)*uv*c2m1/c2**2) + s2m1*nx/c - s2*nx*cm1/c**2;
        cc1m2 = -gam1*((s1m2-rlam3m2)*uv/c2 + (s1-rlam3)*uvm2/c2 - (s1-rlam3)*uv*c2m2/c2**2) + s2m2*nx/c - s2*nx*cm2/c**2;
        cc1m3 = -gam1*((s1m3-rlam3m3)*uv/c2 + (s1-rlam3)*uvm3/c2 - (s1-rlam3)*uv*c2m3/c2**2) + s2m3*nx/c - s2*nx*cm3/c**2;
        cc1m4 = -gam1*((s1m4-rlam3m4)*uv/c2 + (s1-rlam3)*uvm4/c2 - (s1-rlam3)*uv*c2m4/c2**2) + s2m4*nx/c - s2*nx*cm4/c**2;

        cc2   = -gam1*s2*uv/c + (s1-rlam3)*nx;
        cc2m1 = -gam1*(s2m1*uv/c + s2*uvm1/c - s2*uv*cm1/c**2) + (s1m1-rlam3m1)*nx;
        cc2m2 = -gam1*(s2m2*uv/c + s2*uvm2/c - s2*uv*cm2/c**2) + (s1m2-rlam3m2)*nx;
        cc2m3 = -gam1*(s2m3*uv/c + s2*uvm3/c - s2*uv*cm3/c**2) + (s1m3-rlam3m3)*nx;
        cc2m4 = -gam1*(s2m4*uv/c + s2*uvm4/c - s2*uv*cm4/c**2) + (s1m4-rlam3m4)*nx;

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
        
        cc1   = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
        cc1m1 = -gam1*((s1m1-rlam3m1)*vv/c2 + (s1-rlam3)*vvm1/c2 - (s1-rlam3)*vv*c2m1/c2**2) + s2m1*ny/c - s2*ny*cm1/c**2;
        cc1m2 = -gam1*((s1m2-rlam3m2)*vv/c2 + (s1-rlam3)*vvm2/c2 - (s1-rlam3)*vv*c2m2/c2**2) + s2m2*ny/c - s2*ny*cm2/c**2;
        cc1m3 = -gam1*((s1m3-rlam3m3)*vv/c2 + (s1-rlam3)*vvm3/c2 - (s1-rlam3)*vv*c2m3/c2**2) + s2m3*ny/c - s2*ny*cm3/c**2;
        cc1m4 = -gam1*((s1m4-rlam3m4)*vv/c2 + (s1-rlam3)*vvm4/c2 - (s1-rlam3)*vv*c2m4/c2**2) + s2m4*ny/c - s2*ny*cm4/c**2;

        cc2   = -gam1*s2*vv/c+(s1-rlam3)*ny;
        cc2m1 = -gam1*(s2m1*vv/c + s2*vvm1/c - s2*vv*cm1/c**2) + (s1m1-rlam3m1)*ny;
        cc2m2 = -gam1*(s2m2*vv/c + s2*vvm2/c - s2*vv*cm2/c**2) + (s1m2-rlam3m2)*ny;
        cc2m3 = -gam1*(s2m3*vv/c + s2*vvm3/c - s2*vv*cm3/c**2) + (s1m3-rlam3m3)*ny;
        cc2m4 = -gam1*(s2m4*vv/c + s2*vvm4/c - s2*vv*cm4/c**2) + (s1m4-rlam3m4)*ny;

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
        cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/c2**2);
        cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/c2**2);
        cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/c2**2);
        cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/c2**2);

        cc2   = gam1*s2/c;
        cc2m1 = gam1*(s2m1/c - s2*cm1/c**2);
        cc2m2 = gam1*(s2m2/c - s2*cm2/c**2);
        cc2m3 = gam1*(s2m3/c - s2*cm3/c**2);
        cc2m4 = gam1*(s2m4/c - s2*cm4/c**2);

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
