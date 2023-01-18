
// Written by: C. Nguyen & P. Fernandez

void getan2dNEW(double *An, double *Anm, double *uh, double* pg, double *nl, appstruct &app, double *param, int absflag, int numPoints, int nd, int computeJacobian)
{
    double gam, epslm, gam1, signun, ic, Dc2Reg_Dc2;
    double nx, ny;
    double r, ru, rv, rE;
    double r1, r1m1, r1m2, r1m3, r1m4;
    double uv, uvm1, uvm2, uvm3, uvm4;
    double vv, vvm1, vvm2, vvm3, vvm4;
    double Vgx, Vgy, Vgn;
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
    double pi = 3.141592653589793;
    double b = 10.0;

    Int ALEflag = app.ALEflag;

    int    i;
    for (i=0; i<numPoints; i++) {
        r   = uh[0*numPoints+i];
        ru  = uh[1*numPoints+i];
        rv  = uh[2*numPoints+i];
        rE  = uh[3*numPoints+i];

        nx  = nl[0*numPoints+i];
        ny  = nl[1*numPoints+i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
        }

        r1   = 1.0/r;
        uv   = ru*r1;
        vv   = rv*r1;
        E    = rE*r1;
        af   = 0.5*(uv*uv+vv*vv);
        p    = gam1*(rE -r*af);
        h    = E   + p*r1;
        c2   = gam* p*r1;
                Dc2Reg_Dc2 = atan(b*c2)/pi + 1.0/2.0 + (b*c2) / ((1+b*c2*b*c2) * pi);
                c2 = max ( c2*(atan(b*c2)/pi + 1.0/2.0) - atan(b)/pi + 1.0/2.0 , 1.0e-6);
        if (c2 < 0)
            printf("c2 = %g was detected. u = %g, v = %g, p = %g, r = %g, E = %g, rE = %g.\n", c2, uv, vv, p, r, E, rE);
        c    = sqrt(c2);
        un   = uv*nx   + vv*ny;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = r1 + rE * r1m4;

            afm1 = uv * uvm1 + vv * vvm1;
            afm2 = uv * uvm2 + vv * vvm2;
            afm3 = uv * uvm3 + vv * vvm3;
            afm4 = uv * uvm4 + vv * vvm4;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (1.0 - r * afm4);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            
                    c2m1 *= Dc2Reg_Dc2;
                    c2m2 *= Dc2Reg_Dc2;
                    c2m3 *= Dc2Reg_Dc2;
                    c2m4 *= Dc2Reg_Dc2;

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;

            unm1 = uvm1 * nx + vvm1 * ny;
            unm2 = uvm2 * nx + vvm2 * ny;
            unm3 = uvm3 * nx + vvm3 * ny;
            unm4 = uvm4 * nx + vvm4 * ny;
        }

        Vgn = Vgx * nx + Vgy * ny;

        if (absflag==1)
        {
            rlam1   = fabs(un - Vgn + c);
            if (computeJacobian == 1) {
                signun = (un - Vgn + c < 0) ? -1.0 : 1.0;
                rlam1m1 = signun * (unm1 + cm1);
                rlam1m2 = signun * (unm2 + cm2);
                rlam1m3 = signun * (unm3 + cm3);
                rlam1m4 = signun * (unm4 + cm4);
            }

            if (epslm>0) {
                rlam = 0.5*(rlam1*rlam1/(epslm*c)+epslm*c);

                if (computeJacobian == 1) {
                    rlamm1 = rlam1*rlam1m1/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam1*rlam1m2/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam1*rlam1m3/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam1*rlam1m4/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm4);
                }

                ic = rlam1 < epslm*c;
                rlam1 = ic*rlam + (1-ic)*rlam1;

                if (computeJacobian == 1) {
                    rlam1m1 = ic*rlamm1 + (1-ic)*rlam1m1;
                    rlam1m2 = ic*rlamm2 + (1-ic)*rlam1m2;
                    rlam1m3 = ic*rlamm3 + (1-ic)*rlam1m3;
                    rlam1m4 = ic*rlamm4 + (1-ic)*rlam1m4;
                }
            }

            rlam2   = fabs(un - Vgn - c);
            if (computeJacobian == 1) {
                signun = (un - Vgn - c < 0) ? -1.0 : 1.0;
                rlam2m1 = signun * (unm1 - cm1);
                rlam2m2 = signun * (unm2 - cm2);
                rlam2m3 = signun * (unm3 - cm3);
                rlam2m4 = signun * (unm4 - cm4);
            }

            if (epslm>0) {
                rlam = 0.5*(rlam2*rlam2/(epslm*c)+epslm*c);

                if (computeJacobian == 1) {
                    rlamm1 = rlam2*rlam2m1/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam2*rlam2m2/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam2*rlam2m3/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam2*rlam2m4/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm4);
                }

                ic = rlam2 < epslm*c;
                rlam2 = ic*rlam + (1-ic)*rlam2;

                if (computeJacobian == 1) {
                    rlam2m1 = ic*rlamm1 + (1-ic)*rlam2m1;
                    rlam2m2 = ic*rlamm2 + (1-ic)*rlam2m2;
                    rlam2m3 = ic*rlamm3 + (1-ic)*rlam2m3;
                    rlam2m4 = ic*rlamm4 + (1-ic)*rlam2m4;
                }
            }

            rlam3   = fabs(un - Vgn);
            if (computeJacobian == 1) {
                signun = (un - Vgn < 0) ? -1.0 : 1.0;
                rlam3m1 = signun * unm1;
                rlam3m2 = signun * unm2;
                rlam3m3 = signun * unm3;
                rlam3m4 = signun * unm4;
            }

            if (epslm>0) {
                rlam = 0.5*(rlam3*rlam3/(epslm*c)+epslm*c);

                if (computeJacobian == 1) {
                    rlamm1 = rlam3*rlam3m1/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam3*rlam3m2/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam3*rlam3m3/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam3*rlam3m4/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm4);
                }

                ic = rlam3 < epslm*c;
                rlam3 = ic*rlam + (1-ic)*rlam3;

                if (computeJacobian == 1) {
                    rlam3m1 = ic*rlamm1 + (1-ic)*rlam3m1;
                    rlam3m2 = ic*rlamm2 + (1-ic)*rlam3m2;
                    rlam3m3 = ic*rlamm3 + (1-ic)*rlam3m3;
                    rlam3m4 = ic*rlamm4 + (1-ic)*rlam3m4;
                }
            }
        }
        else {
            rlam1   = un - Vgn + c;
            rlam2   = un - Vgn - c;
            rlam3   = un - Vgn;

            if (computeJacobian == 1) {
                rlam1m1 = unm1 + cm1;
                rlam1m2 = unm2 + cm2;
                rlam1m3 = unm3 + cm3;
                rlam1m4 = unm4 + cm4;

                rlam2m1 = unm1 - cm1;
                rlam2m2 = unm2 - cm2;
                rlam2m3 = unm3 - cm3;
                rlam2m4 = unm4 - cm4;

                rlam3m1 = unm1;
                rlam3m2 = unm2;
                rlam3m3 = unm3;
                rlam3m4 = unm4;
            }
        }

        s1 = 0.5*(rlam1   + rlam2);
        s2 = 0.5*(rlam1   - rlam2);

        if (computeJacobian == 1) {
            s1m1    = 0.5*(rlam1m1 + rlam2m1);
            s1m2    = 0.5*(rlam1m2 + rlam2m2);
            s1m3    = 0.5*(rlam1m3 + rlam2m3);
            s1m4    = 0.5*(rlam1m4 + rlam2m4);

            s2m1    = 0.5*(rlam1m1 - rlam2m1);
            s2m2    = 0.5*(rlam1m2 - rlam2m2);
            s2m3    = 0.5*(rlam1m3 - rlam2m3);
            s2m4    = 0.5*(rlam1m4 - rlam2m4);
        }

        cc1   = gam1*(s1-rlam3)*af/c2-(s2*un/c);
        cc2   = gam1*s2*af/c-(s1-rlam3)*un;

        An[0*numPoints+i]  = rlam3+cc1;
        An[1*numPoints+i]  = cc1*uv+cc2*nx;
        An[2*numPoints+i]  = cc1*vv+cc2*ny;
        An[3*numPoints+i]  = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)*af/c2 + (s1-rlam3)*afm1/c2 - (s1-rlam3)*af*c2m1/(c2*c2)) - s2m1*un/c - s2*unm1/c + s2*un*cm1/(c*c);
            cc1m2 = gam1*((s1m2-rlam3m2)*af/c2 + (s1-rlam3)*afm2/c2 - (s1-rlam3)*af*c2m2/(c2*c2)) - s2m2*un/c - s2*unm2/c + s2*un*cm2/(c*c);
            cc1m3 = gam1*((s1m3-rlam3m3)*af/c2 + (s1-rlam3)*afm3/c2 - (s1-rlam3)*af*c2m3/(c2*c2)) - s2m3*un/c - s2*unm3/c + s2*un*cm3/(c*c);
            cc1m4 = gam1*((s1m4-rlam3m4)*af/c2 + (s1-rlam3)*afm4/c2 - (s1-rlam3)*af*c2m4/(c2*c2)) - s2m4*un/c - s2*unm4/c + s2*un*cm4/(c*c);

            cc2m1 = gam1*(s2m1*af/c + s2*afm1/c - s2*af*cm1/(c*c)) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
            cc2m2 = gam1*(s2m2*af/c + s2*afm2/c - s2*af*cm2/(c*c)) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
            cc2m3 = gam1*(s2m3*af/c + s2*afm3/c - s2*af*cm3/(c*c)) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
            cc2m4 = gam1*(s2m4*af/c + s2*afm4/c - s2*af*cm4/(c*c)) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;

            Anm[0 * numPoints + i] = rlam3m1 + cc1m1;
            Anm[1 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[2 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[3 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[16 * numPoints + i] = rlam3m2 + cc1m2;
            Anm[17 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[18 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[19 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[32 * numPoints + i] = rlam3m3 + cc1m3;
            Anm[33 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[34 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[35 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[48 * numPoints + i] = rlam3m4 + cc1m4;
            Anm[49 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[50 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[51 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
        }

        cc1   = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
        cc2   = -gam1*s2*uv/c + (s1-rlam3)*nx;

        An[4*numPoints+i] = cc1;
        An[5*numPoints+i] = rlam3+cc1*uv+cc2*nx;
        An[6*numPoints+i] = cc1*vv+cc2*ny;
        An[7*numPoints+i] = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*uv/c2 + (s1-rlam3)*uvm1/c2 - (s1-rlam3)*uv*c2m1/(c2*c2)) + s2m1*nx/c - s2*nx*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*uv/c2 + (s1-rlam3)*uvm2/c2 - (s1-rlam3)*uv*c2m2/(c2*c2)) + s2m2*nx/c - s2*nx*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*uv/c2 + (s1-rlam3)*uvm3/c2 - (s1-rlam3)*uv*c2m3/(c2*c2)) + s2m3*nx/c - s2*nx*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*uv/c2 + (s1-rlam3)*uvm4/c2 - (s1-rlam3)*uv*c2m4/(c2*c2)) + s2m4*nx/c - s2*nx*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*uv/c + s2*uvm1/c - s2*uv*cm1/(c*c)) + (s1m1-rlam3m1)*nx;
            cc2m2 = -gam1*(s2m2*uv/c + s2*uvm2/c - s2*uv*cm2/(c*c)) + (s1m2-rlam3m2)*nx;
            cc2m3 = -gam1*(s2m3*uv/c + s2*uvm3/c - s2*uv*cm3/(c*c)) + (s1m3-rlam3m3)*nx;
            cc2m4 = -gam1*(s2m4*uv/c + s2*uvm4/c - s2*uv*cm4/(c*c)) + (s1m4-rlam3m4)*nx;

            Anm[4 * numPoints + i] = cc1m1;
            Anm[5 * numPoints + i] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[6 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[7 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[20 * numPoints + i] = cc1m2;
            Anm[21 * numPoints + i] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[22 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[23 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[36 * numPoints + i] = cc1m3;
            Anm[37 * numPoints + i] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[38 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[39 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[52 * numPoints + i] = cc1m4;
            Anm[53 * numPoints + i] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[54 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[55 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
        }

        cc1   = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
        cc2   = -gam1*s2*vv/c+(s1-rlam3)*ny;

        An[8*numPoints+i] = cc1;
        An[9*numPoints+i] = cc1*uv+cc2*nx;
        An[10*numPoints+i] = rlam3+cc1*vv+cc2*ny;
        An[11*numPoints+i] = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*vv/c2 + (s1-rlam3)*vvm1/c2 - (s1-rlam3)*vv*c2m1/(c2*c2)) + s2m1*ny/c - s2*ny*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*vv/c2 + (s1-rlam3)*vvm2/c2 - (s1-rlam3)*vv*c2m2/(c2*c2)) + s2m2*ny/c - s2*ny*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*vv/c2 + (s1-rlam3)*vvm3/c2 - (s1-rlam3)*vv*c2m3/(c2*c2)) + s2m3*ny/c - s2*ny*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*vv/c2 + (s1-rlam3)*vvm4/c2 - (s1-rlam3)*vv*c2m4/(c2*c2)) + s2m4*ny/c - s2*ny*cm4/(c*c);

            cc2m1 = -gam1*(s2m1*vv/c + s2*vvm1/c - s2*vv*cm1/(c*c)) + (s1m1-rlam3m1)*ny;
            cc2m2 = -gam1*(s2m2*vv/c + s2*vvm2/c - s2*vv*cm2/(c*c)) + (s1m2-rlam3m2)*ny;
            cc2m3 = -gam1*(s2m3*vv/c + s2*vvm3/c - s2*vv*cm3/(c*c)) + (s1m3-rlam3m3)*ny;
            cc2m4 = -gam1*(s2m4*vv/c + s2*vvm4/c - s2*vv*cm4/(c*c)) + (s1m4-rlam3m4)*ny;

            Anm[8 * numPoints + i] = cc1m1;
            Anm[9 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[10 * numPoints + i] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[11 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[24 * numPoints + i] = cc1m2;
            Anm[25 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[26 * numPoints + i] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[27 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[40 * numPoints + i] = cc1m3;
            Anm[41 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[42 * numPoints + i] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[43 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[56 * numPoints + i] = cc1m4;
            Anm[57 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[58 * numPoints + i] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[59 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
        }

        cc1   = gam1*(s1-rlam3)/c2;
        cc2   = gam1*s2/c;

        An[12*numPoints+i] = cc1;
        An[13*numPoints+i] = cc1*uv+cc2*nx;
        An[14*numPoints+i] = cc1*vv+cc2*ny;
        An[15*numPoints+i] = rlam3+cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/(c2*c2));
            cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/(c2*c2));
            cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/(c2*c2));
            cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/(c2*c2));

            cc2m1 = gam1*(s2m1/c - s2*cm1/(c*c));
            cc2m2 = gam1*(s2m2/c - s2*cm2/(c*c));
            cc2m3 = gam1*(s2m3/c - s2*cm3/(c*c));
            cc2m4 = gam1*(s2m4/c - s2*cm4/(c*c));

            Anm[12 * numPoints + i] = cc1m1;
            Anm[13 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[14 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[15 * numPoints + i] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[28 * numPoints + i] = cc1m2;
            Anm[29 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[30 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[31 * numPoints + i] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[44 * numPoints + i] = cc1m3;
            Anm[45 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[46 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[47 * numPoints + i] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[60 * numPoints + i] = cc1m4;
            Anm[61 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[62 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[63 * numPoints + i] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
        }
    }
}

void getan3dNEW(double *An, double *Anm, double *uh, double* pg, double *nl, appstruct &app, double *param, int absflag, int numPoints, int nd, int computeJacobian)
{
    double gam, epslm, gam1, signun, ic;
    double nx, ny, nz;
    double r, ru, rv, rw, rE;
    double r1, r1m1, r1m2, r1m3, r1m4, r1m5;
    double uv, uvm1, uvm2, uvm3, uvm4, uvm5;
    double vv, vvm1, vvm2, vvm3, vvm4, vvm5;
    double wv, wvm1, wvm2, wvm3, wvm4, wvm5;
    double Vgx, Vgy, Vgz, Vgn;
    double E, Em1, Em2, Em3, Em4, Em5;
    double af, afm1, afm2, afm3, afm4, afm5;
    double p, pm1, pm2, pm3, pm4, pm5;
    double h, hm1, hm2, hm3, hm4, hm5;
    double s1, s1m1, s1m2, s1m3, s1m4, s1m5;
    double s2, s2m1, s2m2, s2m3, s2m4, s2m5;
    double c2, c2m1, c2m2, c2m3, c2m4, c2m5;
    double c, cm1, cm2, cm3, cm4, cm5;
    double un, unm1, unm2, unm3, unm4, unm5;
    double cc1, cc1m1, cc1m2, cc1m3, cc1m4, cc1m5;
    double cc2, cc2m1, cc2m2, cc2m3, cc2m4, cc2m5;
    double rlam, rlamm1, rlamm2, rlamm3, rlamm4, rlamm5;
    double rlam1, rlam1m1, rlam1m2, rlam1m3, rlam1m4, rlam1m5;
    double rlam2, rlam2m1, rlam2m2, rlam2m3, rlam2m4, rlam2m5;
    double rlam3, rlam3m1, rlam3m2, rlam3m3, rlam3m4, rlam3m5;

    gam   = param[0];
    epslm = param[1];
    gam1 = gam - 1.0;

    Int ALEflag = app.ALEflag;

    int    i;
    for (i=0; i<numPoints; i++) {
        nx   = nl[0*numPoints+i];
        ny   = nl[1*numPoints+i];
        nz   = nl[2*numPoints+i];

        r    = uh[0*numPoints+i];
        ru   = uh[1*numPoints+i];
        rv   = uh[2*numPoints+i];
        rw   = uh[3*numPoints+i];
        rE   = uh[4*numPoints+i];

        if (ALEflag != 0) {
            Vgx = pg[(2 * nd + 0) * numPoints + i];
            Vgy = pg[(2 * nd + 1) * numPoints + i];
            Vgz = pg[(2 * nd + 2) * numPoints + i];
        }
        else {
            Vgx = 0.0;
            Vgy = 0.0;
            Vgz = 0.0;
        }

        r1   = 1.0/r;
        uv   = ru*r1;
        vv   = rv*r1;
        wv   = rw*r1;
        E    = rE*r1;
        af   = 0.5*(uv*uv+vv*vv+wv*wv);
        p    = gam1*(rE -r*af);
        h    = E   + p*r1;
        c2   = gam* p*r1;
        if (c2 < 0)
            printf("c2 = %g was detected. p = %g, r = %g, E = %g, rE = %g.\n", c2, p, r, E, rE);
        c    = sqrt(c2);
        un   = uv*nx   + vv*ny + wv*nz;

        if (computeJacobian == 1) {
            r1m1 = -1.0 / (r * r);
            r1m2 = 0.0;
            r1m3 = 0.0;
            r1m4 = 0.0;
            r1m5 = 0.0;

            uvm1 = ru * r1m1;
            uvm2 = r1 + ru * r1m2;
            uvm3 = ru * r1m3;
            uvm4 = ru * r1m4;
            uvm5 = ru * r1m5;

            vvm1 = rv * r1m1;
            vvm2 = rv * r1m2;
            vvm3 = r1 + rv * r1m3;
            vvm4 = rv * r1m4;
            vvm5 = rv * r1m5;

            wvm1 = rw * r1m1;
            wvm2 = rw * r1m2;
            wvm3 = rw * r1m3;
            wvm4 = r1 + rw * r1m4;
            wvm5 = rw * r1m5;

            Em1 = rE * r1m1;
            Em2 = rE * r1m2;
            Em3 = rE * r1m3;
            Em4 = rE * r1m4;
            Em5 = r1 + rE * r1m5;

            afm1 = uv * uvm1 + vv * vvm1 + wv * wvm1;
            afm2 = uv * uvm2 + vv * vvm2 + wv * wvm2;
            afm3 = uv * uvm3 + vv * vvm3 + wv * wvm3;
            afm4 = uv * uvm4 + vv * vvm4 + wv * wvm4;
            afm5 = uv * uvm5 + vv * vvm5 + wv * wvm5;

            pm1 = gam1 * (-af - r * afm1);
            pm2 = gam1 * (-r * afm2);
            pm3 = gam1 * (-r * afm3);
            pm4 = gam1 * (-r * afm4);
            pm5 = gam1 * (1.0 - r * afm5);

            hm1 = Em1 + pm1 * r1 + p * r1m1;
            hm2 = Em2 + pm2 * r1 + p * r1m2;
            hm3 = Em3 + pm3 * r1 + p * r1m3;
            hm4 = Em4 + pm4 * r1 + p * r1m4;
            hm5 = Em5 + pm5 * r1 + p * r1m5;

            c2m1 = gam * (pm1 * r1 + p * r1m1);
            c2m2 = gam * (pm2 * r1 + p * r1m2);
            c2m3 = gam * (pm3 * r1 + p * r1m3);
            c2m4 = gam * (pm4 * r1 + p * r1m4);
            c2m5 = gam * (pm5 * r1 + p * r1m5);

            cm1 = 0.5 * c2m1 / c;
            cm2 = 0.5 * c2m2 / c;
            cm3 = 0.5 * c2m3 / c;
            cm4 = 0.5 * c2m4 / c;
            cm5 = 0.5 * c2m5 / c;

            unm1 = uvm1 * nx + vvm1 * ny + wvm1 * nz;
            unm2 = uvm2 * nx + vvm2 * ny + wvm2 * nz;
            unm3 = uvm3 * nx + vvm3 * ny + wvm3 * nz;
            unm4 = uvm4 * nx + vvm4 * ny + wvm4 * nz;
            unm5 = uvm5 * nx + vvm5 * ny + wvm5 * nz;
        }

        Vgn = Vgx * nx + Vgy * ny + Vgz * nz;

        if (absflag==1)
        {
            rlam1   = fabs(un - Vgn + c);
            if (computeJacobian == 1) {
                signun = (un - Vgn + c < 0) ? -1.0 : 1.0;
                rlam1m1 = signun * (unm1 + cm1);
                rlam1m2 = signun * (unm2 + cm2);
                rlam1m3 = signun * (unm3 + cm3);
                rlam1m4 = signun * (unm4 + cm4);
                rlam1m5 = signun * (unm5 + cm5);
            }

            if (epslm>0) {
                rlam = 0.5*(rlam1*rlam1/(epslm*c)+epslm*c);
                if (computeJacobian == 1) {
                    rlamm1 = rlam1*rlam1m1/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam1*rlam1m2/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam1*rlam1m3/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam1*rlam1m4/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm4);
                    rlamm5 = rlam1*rlam1m5/(epslm*c) + 0.5*(1 - rlam1*rlam1/(epslm*c*epslm*c))*(epslm*cm5);
                }

                ic = rlam1 < epslm*c;
                rlam1 = ic*rlam + (1-ic)*rlam1;

                if (computeJacobian == 1) {
                    rlam1m1 = ic*rlamm1 + (1-ic)*rlam1m1;
                    rlam1m2 = ic*rlamm2 + (1-ic)*rlam1m2;
                    rlam1m3 = ic*rlamm3 + (1-ic)*rlam1m3;
                    rlam1m4 = ic*rlamm4 + (1-ic)*rlam1m4;
                    rlam1m5 = ic*rlamm5 + (1-ic)*rlam1m5;
                }
            }

            rlam2   = fabs(un - Vgn - c);
            if (computeJacobian == 1) {
                signun = (un - Vgn - c < 0) ? -1.0 : 1.0;
                rlam2m1 = signun * (unm1 - cm1);
                rlam2m2 = signun * (unm2 - cm2);
                rlam2m3 = signun * (unm3 - cm3);
                rlam2m4 = signun * (unm4 - cm4);
                rlam2m5 = signun * (unm5 - cm5);
            }

            if (epslm>0) {
                rlam = 0.5*(rlam2*rlam2/(epslm*c)+epslm*c);
                if (computeJacobian == 1) {
                    rlamm1 = rlam2*rlam2m1/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam2*rlam2m2/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam2*rlam2m3/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam2*rlam2m4/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm4);
                    rlamm5 = rlam2*rlam2m5/(epslm*c) + 0.5*(1 - rlam2*rlam2/(epslm*c*epslm*c))*(epslm*cm5);
                }

                ic = rlam2 < epslm*c;
                rlam2 = ic*rlam + (1-ic)*rlam2;

                if (computeJacobian == 1) {
                    rlam2m1 = ic*rlamm1 + (1-ic)*rlam2m1;
                    rlam2m2 = ic*rlamm2 + (1-ic)*rlam2m2;
                    rlam2m3 = ic*rlamm3 + (1-ic)*rlam2m3;
                    rlam2m4 = ic*rlamm4 + (1-ic)*rlam2m4;
                    rlam2m5 = ic*rlamm5 + (1-ic)*rlam2m5;
                }
            }

            rlam3   = fabs(un - Vgn);
            if (computeJacobian == 1) {
                signun = (un - Vgn < 0) ? -1.0 : 1.0;
                rlam3m1 = signun * unm1;
                rlam3m2 = signun * unm2;
                rlam3m3 = signun * unm3;
                rlam3m4 = signun * unm4;
                rlam3m5 = signun * unm5;
            }

            if (epslm>0) {
                rlam = 0.5*(rlam3*rlam3/(epslm*c)+epslm*c);
                if (computeJacobian == 1) {
                    rlamm1 = rlam3*rlam3m1/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm1);
                    rlamm2 = rlam3*rlam3m2/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm2);
                    rlamm3 = rlam3*rlam3m3/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm3);
                    rlamm4 = rlam3*rlam3m4/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm4);
                    rlamm5 = rlam3*rlam3m5/(epslm*c) + 0.5*(1 - rlam3*rlam3/(epslm*c*epslm*c))*(epslm*cm5);
                }

                ic = rlam3 < epslm*c;
                rlam3 = ic*rlam + (1-ic)*rlam3;

                if (computeJacobian == 1) {
                    rlam3m1 = ic*rlamm1 + (1-ic)*rlam3m1;
                    rlam3m2 = ic*rlamm2 + (1-ic)*rlam3m2;
                    rlam3m3 = ic*rlamm3 + (1-ic)*rlam3m3;
                    rlam3m4 = ic*rlamm4 + (1-ic)*rlam3m4;
                    rlam3m5 = ic*rlamm5 + (1-ic)*rlam3m5;
                }
            }
        }
        else {
            rlam1   = un - Vgn + c;
            rlam2   = un - Vgn - c;
            rlam3   = un - Vgn;

            if (computeJacobian == 1) {
                rlam1m1 = unm1+cm1;
                rlam1m2 = unm2+cm2;
                rlam1m3 = unm3+cm3;
                rlam1m4 = unm4+cm4;
                rlam1m5 = unm5+cm5;

                rlam2m1 = unm1-cm1;
                rlam2m2 = unm2-cm2;
                rlam2m3 = unm3-cm3;
                rlam2m4 = unm4-cm4;
                rlam2m5 = unm5-cm5;

                rlam3m1 = unm1;
                rlam3m2 = unm2;
                rlam3m3 = unm3;
                rlam3m4 = unm4;
                rlam3m5 = unm5;
            }
        }

        s1      = 0.5*(rlam1   + rlam2);
        s2      = 0.5*(rlam1   - rlam2);

        if (computeJacobian == 1) {
            s1m1 = 0.5 * (rlam1m1 + rlam2m1);
            s1m2 = 0.5 * (rlam1m2 + rlam2m2);
            s1m3 = 0.5 * (rlam1m3 + rlam2m3);
            s1m4 = 0.5 * (rlam1m4 + rlam2m4);
            s1m5 = 0.5 * (rlam1m5 + rlam2m5);

            s2m1 = 0.5 * (rlam1m1 - rlam2m1);
            s2m2 = 0.5 * (rlam1m2 - rlam2m2);
            s2m3 = 0.5 * (rlam1m3 - rlam2m3);
            s2m4 = 0.5 * (rlam1m4 - rlam2m4);
            s2m5 = 0.5 * (rlam1m5 - rlam2m5);
        }

        cc1   = gam1*(s1-rlam3)*af/c2-(s2*un/c);
        cc2   = gam1*s2*af/c-(s1-rlam3)*un;

        An[0*numPoints+i]  = rlam3+cc1;
        An[1*numPoints+i]  = cc1*uv+cc2*nx;
        An[2*numPoints+i]  = cc1*vv+cc2*ny;
        An[3*numPoints+i]  = cc1*wv+cc2*nz;
        An[4*numPoints+i]  = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc2m1 = gam1*(s2m1*af/c + s2*afm1/c - s2*af*cm1/(c*c)) - (s1m1-rlam3m1)*un - (s1-rlam3)*unm1;
            cc2m2 = gam1*(s2m2*af/c + s2*afm2/c - s2*af*cm2/(c*c)) - (s1m2-rlam3m2)*un - (s1-rlam3)*unm2;
            cc2m3 = gam1*(s2m3*af/c + s2*afm3/c - s2*af*cm3/(c*c)) - (s1m3-rlam3m3)*un - (s1-rlam3)*unm3;
            cc2m4 = gam1*(s2m4*af/c + s2*afm4/c - s2*af*cm4/(c*c)) - (s1m4-rlam3m4)*un - (s1-rlam3)*unm4;
            cc2m5 = gam1*(s2m5*af/c + s2*afm5/c - s2*af*cm5/(c*c)) - (s1m5-rlam3m5)*un - (s1-rlam3)*unm5;

            cc1m1 = gam1*((s1m1-rlam3m1)*af/c2 + (s1-rlam3)*afm1/c2 - (s1-rlam3)*af*c2m1/(c2*c2)) - s2m1*un/c - s2*unm1/c + s2*un*cm1/(c*c);
            cc1m2 = gam1*((s1m2-rlam3m2)*af/c2 + (s1-rlam3)*afm2/c2 - (s1-rlam3)*af*c2m2/(c2*c2)) - s2m2*un/c - s2*unm2/c + s2*un*cm2/(c*c);
            cc1m3 = gam1*((s1m3-rlam3m3)*af/c2 + (s1-rlam3)*afm3/c2 - (s1-rlam3)*af*c2m3/(c2*c2)) - s2m3*un/c - s2*unm3/c + s2*un*cm3/(c*c);
            cc1m4 = gam1*((s1m4-rlam3m4)*af/c2 + (s1-rlam3)*afm4/c2 - (s1-rlam3)*af*c2m4/(c2*c2)) - s2m4*un/c - s2*unm4/c + s2*un*cm4/(c*c);
            cc1m5 = gam1*((s1m5-rlam3m5)*af/c2 + (s1-rlam3)*afm5/c2 - (s1-rlam3)*af*c2m5/(c2*c2)) - s2m5*un/c - s2*unm5/c + s2*un*cm5/(c*c);

            Anm[0 * numPoints + i] = rlam3m1 + cc1m1;
            Anm[1 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[2 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[3 * numPoints + i] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            Anm[4 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[25 * numPoints + i] = rlam3m2 + cc1m2;
            Anm[26 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[27 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[28 * numPoints + i] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            Anm[29 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[50 * numPoints + i] = rlam3m3 + cc1m3;
            Anm[51 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[52 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[53 * numPoints + i] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            Anm[54 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[75 * numPoints + i] = rlam3m4 + cc1m4;
            Anm[76 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[77 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[78 * numPoints + i] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            Anm[79 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            Anm[100 * numPoints + i] = rlam3m5 + cc1m5;
            Anm[101 * numPoints + i] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            Anm[102 * numPoints + i] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            Anm[103 * numPoints + i] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            Anm[104 * numPoints + i] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;
        }

        cc1   = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
        cc2   = -gam1*s2*uv/c + (s1-rlam3)*nx;

        An[5*numPoints+i]  = cc1;
        An[6*numPoints+i]  = rlam3+cc1*uv+cc2*nx;
        An[7*numPoints+i]  = cc1*vv+cc2*ny;
        An[8*numPoints+i]  = cc1*wv+cc2*nz;
        An[9*numPoints+i]  = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*uv/c2 + (s1-rlam3)*uvm1/c2 - (s1-rlam3)*uv*c2m1/(c2*c2)) + s2m1*nx/c - s2*nx*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*uv/c2 + (s1-rlam3)*uvm2/c2 - (s1-rlam3)*uv*c2m2/(c2*c2)) + s2m2*nx/c - s2*nx*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*uv/c2 + (s1-rlam3)*uvm3/c2 - (s1-rlam3)*uv*c2m3/(c2*c2)) + s2m3*nx/c - s2*nx*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*uv/c2 + (s1-rlam3)*uvm4/c2 - (s1-rlam3)*uv*c2m4/(c2*c2)) + s2m4*nx/c - s2*nx*cm4/(c*c);
            cc1m5 = -gam1*((s1m5-rlam3m5)*uv/c2 + (s1-rlam3)*uvm5/c2 - (s1-rlam3)*uv*c2m5/(c2*c2)) + s2m5*nx/c - s2*nx*cm5/(c*c);

            cc2m1 = -gam1*(s2m1*uv/c + s2*uvm1/c - s2*uv*cm1/(c*c)) + (s1m1-rlam3m1)*nx;
            cc2m2 = -gam1*(s2m2*uv/c + s2*uvm2/c - s2*uv*cm2/(c*c)) + (s1m2-rlam3m2)*nx;
            cc2m3 = -gam1*(s2m3*uv/c + s2*uvm3/c - s2*uv*cm3/(c*c)) + (s1m3-rlam3m3)*nx;
            cc2m4 = -gam1*(s2m4*uv/c + s2*uvm4/c - s2*uv*cm4/(c*c)) + (s1m4-rlam3m4)*nx;
            cc2m5 = -gam1*(s2m5*uv/c + s2*uvm5/c - s2*uv*cm5/(c*c)) + (s1m5-rlam3m5)*nx;

            Anm[5 * numPoints + i] = cc1m1;
            Anm[6 * numPoints + i] = rlam3m1 + cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[7 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[8 * numPoints + i] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            Anm[9 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[30 * numPoints + i] = cc1m2;
            Anm[31 * numPoints + i] = rlam3m2 + cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[32 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[33 * numPoints + i] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            Anm[34 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[55 * numPoints + i] = cc1m3;
            Anm[56 * numPoints + i] = rlam3m3 + cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[57 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[58 * numPoints + i] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            Anm[59 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[80 * numPoints + i] = cc1m4;
            Anm[81 * numPoints + i] = rlam3m4 + cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[82 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[83 * numPoints + i] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            Anm[84 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            Anm[105 * numPoints + i] = cc1m5;
            Anm[106 * numPoints + i] = rlam3m5 + cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            Anm[107 * numPoints + i] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            Anm[108 * numPoints + i] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            Anm[109 * numPoints + i] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;
        }

        cc1   = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
        cc2   = -gam1*s2*vv/c+(s1-rlam3)*ny;

        An[10*numPoints+i]  = cc1;
        An[11*numPoints+i]  = cc1*uv+cc2*nx;
        An[12*numPoints+i]  = rlam3+cc1*vv+cc2*ny;
        An[13*numPoints+i]  = cc1*wv+cc2*nz;
        An[14*numPoints+i]  = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*vv/c2 + (s1-rlam3)*vvm1/c2 - (s1-rlam3)*vv*c2m1/(c2*c2)) + s2m1*ny/c - s2*ny*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*vv/c2 + (s1-rlam3)*vvm2/c2 - (s1-rlam3)*vv*c2m2/(c2*c2)) + s2m2*ny/c - s2*ny*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*vv/c2 + (s1-rlam3)*vvm3/c2 - (s1-rlam3)*vv*c2m3/(c2*c2)) + s2m3*ny/c - s2*ny*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*vv/c2 + (s1-rlam3)*vvm4/c2 - (s1-rlam3)*vv*c2m4/(c2*c2)) + s2m4*ny/c - s2*ny*cm4/(c*c);
            cc1m5 = -gam1*((s1m5-rlam3m5)*vv/c2 + (s1-rlam3)*vvm5/c2 - (s1-rlam3)*vv*c2m5/(c2*c2)) + s2m5*ny/c - s2*ny*cm5/(c*c);

            cc2m1 = -gam1*(s2m1*vv/c + s2*vvm1/c - s2*vv*cm1/(c*c)) + (s1m1-rlam3m1)*ny;
            cc2m2 = -gam1*(s2m2*vv/c + s2*vvm2/c - s2*vv*cm2/(c*c)) + (s1m2-rlam3m2)*ny;
            cc2m3 = -gam1*(s2m3*vv/c + s2*vvm3/c - s2*vv*cm3/(c*c)) + (s1m3-rlam3m3)*ny;
            cc2m4 = -gam1*(s2m4*vv/c + s2*vvm4/c - s2*vv*cm4/(c*c)) + (s1m4-rlam3m4)*ny;
            cc2m5 = -gam1*(s2m5*vv/c + s2*vvm5/c - s2*vv*cm5/(c*c)) + (s1m5-rlam3m5)*ny;

            Anm[10 * numPoints + i] = cc1m1;
            Anm[11 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[12 * numPoints + i] = rlam3m1 + cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[13 * numPoints + i] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            Anm[14 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[35 * numPoints + i] = cc1m2;
            Anm[36 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[37 * numPoints + i] = rlam3m2 + cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[38 * numPoints + i] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            Anm[39 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[60 * numPoints + i] = cc1m3;
            Anm[61 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[62 * numPoints + i] = rlam3m3 + cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[63 * numPoints + i] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            Anm[64 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[85 * numPoints + i] = cc1m4;
            Anm[86 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[87 * numPoints + i] = rlam3m4 + cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[88 * numPoints + i] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            Anm[89 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            Anm[110 * numPoints + i] = cc1m5;
            Anm[111 * numPoints + i] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            Anm[112 * numPoints + i] = rlam3m5 + cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            Anm[113 * numPoints + i] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            Anm[114 * numPoints + i] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;
        }

        cc1   = -gam1*(s1-rlam3)*wv/c2+(s2*nz/c);
        cc2   = -gam1*s2*wv/c+(s1-rlam3)*nz;

        An[15*numPoints+i]  = cc1;
        An[16*numPoints+i]  = cc1*uv+cc2*nx;
        An[17*numPoints+i]  = cc1*vv+cc2*ny;
        An[18*numPoints+i]  = rlam3+cc1*wv+cc2*nz;
        An[19*numPoints+i]  = cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = -gam1*((s1m1-rlam3m1)*wv/c2 + (s1-rlam3)*wvm1/c2 - (s1-rlam3)*wv*c2m1/(c2*c2)) + s2m1*nz/c - s2*nz*cm1/(c*c);
            cc1m2 = -gam1*((s1m2-rlam3m2)*wv/c2 + (s1-rlam3)*wvm2/c2 - (s1-rlam3)*wv*c2m2/(c2*c2)) + s2m2*nz/c - s2*nz*cm2/(c*c);
            cc1m3 = -gam1*((s1m3-rlam3m3)*wv/c2 + (s1-rlam3)*wvm3/c2 - (s1-rlam3)*wv*c2m3/(c2*c2)) + s2m3*nz/c - s2*nz*cm3/(c*c);
            cc1m4 = -gam1*((s1m4-rlam3m4)*wv/c2 + (s1-rlam3)*wvm4/c2 - (s1-rlam3)*wv*c2m4/(c2*c2)) + s2m4*nz/c - s2*nz*cm4/(c*c);
            cc1m5 = -gam1*((s1m5-rlam3m5)*wv/c2 + (s1-rlam3)*wvm5/c2 - (s1-rlam3)*wv*c2m5/(c2*c2)) + s2m5*nz/c - s2*nz*cm5/(c*c);

            cc2m1 = -gam1*(s2m1*wv/c + s2*wvm1/c - s2*wv*cm1/(c*c)) + (s1m1-rlam3m1)*nz;
            cc2m2 = -gam1*(s2m2*wv/c + s2*wvm2/c - s2*wv*cm2/(c*c)) + (s1m2-rlam3m2)*nz;
            cc2m3 = -gam1*(s2m3*wv/c + s2*wvm3/c - s2*wv*cm3/(c*c)) + (s1m3-rlam3m3)*nz;
            cc2m4 = -gam1*(s2m4*wv/c + s2*wvm4/c - s2*wv*cm4/(c*c)) + (s1m4-rlam3m4)*nz;
            cc2m5 = -gam1*(s2m5*wv/c + s2*wvm5/c - s2*wv*cm5/(c*c)) + (s1m5-rlam3m5)*nz;

            Anm[15 * numPoints + i] = cc1m1;
            Anm[16 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[17 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[18 * numPoints + i] = rlam3m1 + cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            Anm[19 * numPoints + i] = cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[40 * numPoints + i] = cc1m2;
            Anm[41 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[42 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[43 * numPoints + i] = rlam3m2 + cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            Anm[44 * numPoints + i] = cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[65 * numPoints + i] = cc1m3;
            Anm[66 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[67 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[68 * numPoints + i] = rlam3m3 + cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            Anm[69 * numPoints + i] = cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[90 * numPoints + i] = cc1m4;
            Anm[91 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[92 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[93 * numPoints + i] = rlam3m4 + cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            Anm[94 * numPoints + i] = cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            Anm[115 * numPoints + i] = cc1m5;
            Anm[116 * numPoints + i] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            Anm[117 * numPoints + i] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            Anm[118 * numPoints + i] = rlam3m5 + cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            Anm[119 * numPoints + i] = cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;
        }

        cc1   = gam1*(s1-rlam3)/c2;
        cc2   = gam1*s2/c;

        An[20*numPoints+i]  = cc1;
        An[21*numPoints+i]  = cc1*uv+cc2*nx;
        An[22*numPoints+i]  = cc1*vv+cc2*ny;
        An[23*numPoints+i]  = cc1*wv+cc2*nz;
        An[24*numPoints+i]  = rlam3+cc1*h+cc2*un;

        if (computeJacobian == 1) {
            cc1m1 = gam1*((s1m1-rlam3m1)/c2 - (s1-rlam3)*c2m1/(c2*c2));
            cc1m2 = gam1*((s1m2-rlam3m2)/c2 - (s1-rlam3)*c2m2/(c2*c2));
            cc1m3 = gam1*((s1m3-rlam3m3)/c2 - (s1-rlam3)*c2m3/(c2*c2));
            cc1m4 = gam1*((s1m4-rlam3m4)/c2 - (s1-rlam3)*c2m4/(c2*c2));
            cc1m5 = gam1*((s1m5-rlam3m5)/c2 - (s1-rlam3)*c2m5/(c2*c2));

            cc2m1 = gam1*(s2m1/c - s2*cm1/(c*c));
            cc2m2 = gam1*(s2m2/c - s2*cm2/(c*c));
            cc2m3 = gam1*(s2m3/c - s2*cm3/(c*c));
            cc2m4 = gam1*(s2m4/c - s2*cm4/(c*c));
            cc2m5 = gam1*(s2m5/c - s2*cm5/(c*c));

            Anm[20 * numPoints + i] = cc1m1;
            Anm[21 * numPoints + i] = cc1m1 * uv + cc1 * uvm1 + cc2m1 * nx;
            Anm[22 * numPoints + i] = cc1m1 * vv + cc1 * vvm1 + cc2m1 * ny;
            Anm[23 * numPoints + i] = cc1m1 * wv + cc1 * wvm1 + cc2m1 * nz;
            Anm[24 * numPoints + i] = rlam3m1 + cc1m1 * h + cc1 * hm1 + cc2m1 * un + cc2 * unm1;
            Anm[45 * numPoints + i] = cc1m2;
            Anm[46 * numPoints + i] = cc1m2 * uv + cc1 * uvm2 + cc2m2 * nx;
            Anm[47 * numPoints + i] = cc1m2 * vv + cc1 * vvm2 + cc2m2 * ny;
            Anm[48 * numPoints + i] = cc1m2 * wv + cc1 * wvm2 + cc2m2 * nz;
            Anm[49 * numPoints + i] = rlam3m2 + cc1m2 * h + cc1 * hm2 + cc2m2 * un + cc2 * unm2;
            Anm[70 * numPoints + i] = cc1m3;
            Anm[71 * numPoints + i] = cc1m3 * uv + cc1 * uvm3 + cc2m3 * nx;
            Anm[72 * numPoints + i] = cc1m3 * vv + cc1 * vvm3 + cc2m3 * ny;
            Anm[73 * numPoints + i] = cc1m3 * wv + cc1 * wvm3 + cc2m3 * nz;
            Anm[74 * numPoints + i] = rlam3m3 + cc1m3 * h + cc1 * hm3 + cc2m3 * un + cc2 * unm3;
            Anm[95 * numPoints + i] = cc1m4;
            Anm[96 * numPoints + i] = cc1m4 * uv + cc1 * uvm4 + cc2m4 * nx;
            Anm[97 * numPoints + i] = cc1m4 * vv + cc1 * vvm4 + cc2m4 * ny;
            Anm[98 * numPoints + i] = cc1m4 * wv + cc1 * wvm4 + cc2m4 * nz;
            Anm[99 * numPoints + i] = rlam3m4 + cc1m4 * h + cc1 * hm4 + cc2m4 * un + cc2 * unm4;
            Anm[120 * numPoints + i] = cc1m5;
            Anm[121 * numPoints + i] = cc1m5 * uv + cc1 * uvm5 + cc2m5 * nx;
            Anm[122 * numPoints + i] = cc1m5 * vv + cc1 * vvm5 + cc2m5 * ny;
            Anm[123 * numPoints + i] = cc1m5 * wv + cc1 * wvm5 + cc2m5 * nz;
            Anm[124 * numPoints + i] = rlam3m5 + cc1m5 * h + cc1 * hm5 + cc2m5 * un + cc2 * unm5;
        }
    }
}

void getanNEW(double *An, double *Anm, double *uh, double* pg, double *nl, appstruct &app, double *param,
              int absflag, int numPoints, int nd, int computeJacobian)
{
    switch (nd) {
        case 2 :
            getan2dNEW(An, Anm, uh, pg, nl, app, param, absflag, numPoints, nd, computeJacobian);
            break;
        case 3 :
            getan3dNEW(An, Anm, uh, pg, nl, app, param, absflag, numPoints, nd, computeJacobian);
            break;
        default :
            break;
    }
}
