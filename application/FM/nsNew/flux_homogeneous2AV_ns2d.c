
// Written by: C. Nguyen & P. Fernandez

void flux_homogeneous2AV_ns2d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
                              double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    // This function assumes the non-dimensionalization is such that v_inf = 1

    double r, ru, rv, rE, rx, ry, rux, rvy, u, v, ux, vy, p;
    double p_reg, D_p_reg_D_p, r_reg, D_r_reg_D_r, trash;
    double c, g, x, l, wallDistance, ff, D_ff_D_x;

    double param1 = param[0];
    double param2 = param[1];
    double param3 = param[2];
    double param4 = param[3];
    double param5 = param[4];
    double param6 = param[5];
    double c_inf = 1/param5;

    double alpha = 1.0e4;
    double beta = 1.0e-2;
    double alpha_ff = 10.0;
    double beta_ff = 0.5;
    double eps_v = 1.0e-8;
    double rRegMin = 1.0e-2;
    double pRegMin = 1.0e-3;
    double rampFactor = app.rampFactor;
    double epsilon0 = 1.0;
    double h0 = 1.0e-3;

    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];

        wallDistance = pg[3*ng+i];
//        error("Need to determine the initial index of wallDistance\n");

        double u1 = udg[0*ng+i];
        double u2 = udg[1*ng+i];
        double u3 = udg[2*ng+i];
        double u4 = udg[3*ng+i];
        double u5 = udg[4*ng+i];
        double u6 = udg[5*ng+i];
        double u7 = udg[6*ng+i];
        double u8 = udg[7*ng+i];
        double u9 = udg[8*ng+i];
        double u10 = udg[9*ng+i];
        double u11 = udg[10*ng+i];
        double u12 = udg[11*ng+i];

        r = u1;
        ru = u2;
        rv = u3;
        rE = u4;
        rx = u5;
        rux = u6;
        ry = u9;
        rvy = u11;

        u = ru/r;
        v = rv/r;
        p = (param1-1)*(rE-r*(u*u+v*v+eps_v*eps_v)/2);
        surrogateMaximum(&p_reg, &trash, &D_p_reg_D_p, pRegMin, p, alpha);
        surrogateMaximum(&r_reg, &trash, &D_r_reg_D_r, rRegMin, r, alpha);

        ux = (rux - u*rx)/r;
        vy = (rvy - v*ry)/r;

        c = sqrt((param1*p_reg)/r_reg);
        l = min(h0,10*wallDistance);
        g = c_inf*sqrt(1.0 + 0.5*log(1.0 + exp(2*((c/c_inf)*(c/c_inf)-1.0))));
        x = (l*(ux+vy))/g;
        ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;

        double t2 = 1.0/(u1*u1);
        double t3 = 1.0/u1;
        double t4 = param1*(1.0/2.0);
        double t5 = t4-1.0/2.0;
        double t6 = u2*u2;
        double t7 = t2*t6;
        double t8 = u3*u3;
        double t9 = t2*t8;
        double t10 = eps_v*eps_v;
        double t11 = t7+t9+t10;

        f[0*ng+i] += epsilon0*ff*rampFactor*u5;
        f[1*ng+i] += epsilon0*ff*rampFactor*u6;
        f[2*ng+i] += epsilon0*ff*rampFactor*u7;
        f[3*ng+i] += -epsilon0*ff*rampFactor*(t5*(t11*u5+t3*u2*(u6-t3*u2*u5)*2.0+t3*u3*(u7-t3*u3*u5)*2.0)-param1*u8);
        f[4*ng+i] += epsilon0*ff*rampFactor*u9;
        f[5*ng+i] += epsilon0*ff*rampFactor*u10;
        f[6*ng+i] += epsilon0*ff*rampFactor*u11;
        f[7*ng+i] += -epsilon0*ff*rampFactor*(t5*(t11*u9+t3*u2*(u10-t3*u2*u9)*2.0+t3*u3*(u11-t3*u3*u9)*2.0)-param1*u12);
    }

    if (computeJacobian == 1) {

        for (int i = 0; i <ng; i++) {
            double x1 = pg[0*ng+i];
            double x2 = pg[1*ng+i];

            wallDistance = pg[3*ng+i];
//            error("Need to determine the initial index of wallDistance\n");

            double u1 = udg[0*ng+i];
            double u2 = udg[1*ng+i];
            double u3 = udg[2*ng+i];
            double u4 = udg[3*ng+i];
            double u5 = udg[4*ng+i];
            double u6 = udg[5*ng+i];
            double u7 = udg[6*ng+i];
            double u8 = udg[7*ng+i];
            double u9 = udg[8*ng+i];
            double u10 = udg[9*ng+i];
            double u11 = udg[10*ng+i];
            double u12 = udg[11*ng+i];

            r = u1;
            ru = u2;
            rv = u3;
            rE = u4;
            rx = u5;
            rux = u6;
            ry = u9;
            rvy = u11;

            u = ru/r;
            v = rv/r;
            p = (param1-1)*(rE-r*(u*u+v*v+eps_v*eps_v)/2);
            surrogateMaximum(&p_reg, &trash, &D_p_reg_D_p, pRegMin, p, alpha);
            surrogateMaximum(&r_reg, &trash, &D_r_reg_D_r, rRegMin, r, alpha);

            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;

            c = sqrt((param1*p_reg)/r_reg);
            l = min(h0,10*wallDistance);
            g = c_inf*sqrt(1.0 + 0.5*log(1.0 + exp(2*((c/c_inf)*(c/c_inf)-1.0))));
            x = (l*(ux+vy))/g;
            ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;
            D_ff_D_x = exp(alpha_ff*(x-beta_ff))/(1 + exp(alpha_ff*(x-beta_ff)));

            double t2 = 1.0/(u1*u1);
            double t3 = 1.0/u1;
            double t4 = 1.0/(u1*u1*u1);
            double t5 = param5*param5;
            double t6 = 1.0/r_reg;
            double t7 = p_reg*param1*t5*t6*2.0;
            double t8 = t7-2.0;
            double t9 = exp(t8);
            double t10 = t9+1.0;
            double t11 = log(t10);
            double t12 = t11*(1.0/2.0);
            double t13 = t12+1.0;
            double t18 = t3*u2*u5;
            double t14 = -t18+u6;
            double t20 = t3*u3*u9;
            double t15 = -t20+u11;
            double t16 = 1.0/t10;
            double t17 = 1.0/pow(t13,3.0/2.0);
            double t19 = t3*t14;
            double t21 = t3*t15;
            double t22 = t19+t21;
            double t23 = u2*u2;
            double t24 = u3*u3;
            double t25 = 1.0/sqrt(t13);
            double t26 = t2*t14;
            double t27 = t2*t15;
            double t42 = t4*u2*u5;
            double t43 = t4*u3*u9;
            double t28 = t26+t27-t42-t43;
            double t29 = 1.0/(r_reg*r_reg);
            double t30 = D_r_reg_D_r*l*p_reg*param1*param5*t5*t9*t16*t17*t22*t29*(1.0/2.0);
            double t31 = param1-1.0;
            double t32 = t2*t23*(1.0/2.0);
            double t33 = t2*t24*(1.0/2.0);
            double t34 = t4*t23*2.0;
            double t35 = t4*t24*2.0;
            double t36 = t34+t35;
            double t37 = eps_v*eps_v;
            double t38 = t37*(1.0/2.0);
            double t45 = t36*u1*(1.0/2.0);
            double t39 = t32+t33+t38-t45;
            double t40 = D_p_reg_D_p*l*param1*param5*t5*t6*t9*t16*t17*t22*t31*t39*(1.0/2.0);
            double t44 = l*param5*t25*t28;
            double t41 = t30+t40-t44;
            double t46 = param1*(1.0/2.0);
            double t47 = t46-1.0/2.0;
            double t58 = t3*u3*u5;
            double t48 = -t58+u7;
            double t49 = t2*t23;
            double t50 = t2*t24;
            double t51 = t37+t49+t50;
            double t64 = t3*u2*u9;
            double t52 = -t64+u10;
            double t53 = l*param5*t2*t25*u5;
            double t55 = D_p_reg_D_p*l*param1*param5*t3*t5*t6*t9*t16*t17*t22*t31*u2*(1.0/2.0);
            double t54 = t53-t55;
            double t56 = t51*u5;
            double t57 = t3*t14*u2*2.0;
            double t59 = t3*t48*u3*2.0;
            double t60 = t56+t57+t59;
            double t61 = param1*u8;
            double t73 = t47*t60;
            double t62 = t61-t73;
            double t63 = t51*u9;
            double t65 = t3*t52*u2*2.0;
            double t66 = t3*t15*u3*2.0;
            double t67 = t63+t65+t66;
            double t68 = param1*u12;
            double t74 = t47*t67;
            double t69 = t68-t74;
            double t70 = l*param5*t2*t25*u9;
            double t72 = D_p_reg_D_p*l*param1*param5*t3*t5*t6*t9*t16*t17*t22*t31*u3*(1.0/2.0);
            double t71 = t70-t72;
            double t75 = epsilon0*ff*rampFactor;
            double t76 = -t37+t49+t50;
            double t77 = epsilon0*ff*rampFactor*t47*t76;
            double t78 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u5;
            double t79 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u6;
            double t80 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u7;
            double t81 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*(t61-t73);
            double t82 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u9;
            double t83 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u10;
            double t84 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*u11;
            double t85 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t25*(t68-t74);
            double t86 = epsilon0*ff*param1*rampFactor;

            f_udg[0*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u5;
            f_udg[1*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u6;
            f_udg[2*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u7;
            f_udg[3*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*t62+epsilon0*ff*rampFactor*t47*(t36*u5+t2*t14*u2*2.0-t4*t23*u5*2.0-t4*t24*u5*2.0+t2*t48*u3*2.0);
            f_udg[4*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u9;
            f_udg[5*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u10;
            f_udg[6*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*u11;
            f_udg[7*ng+i] += D_ff_D_x*epsilon0*rampFactor*t41*t69+epsilon0*ff*rampFactor*t47*(t36*u9+t2*t15*u3*2.0-t4*t23*u9*2.0-t4*t24*u9*2.0+t2*t52*u2*2.0);
            f_udg[8*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u5;
            f_udg[9*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u6;
            f_udg[10*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u7;
            f_udg[11*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*t62-epsilon0*ff*rampFactor*t3*t14*t47*2.0;
            f_udg[12*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u9;
            f_udg[13*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u10;
            f_udg[14*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*u11;
            f_udg[15*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*t69-epsilon0*ff*rampFactor*t3*t47*t52*2.0;
            f_udg[16*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u5;
            f_udg[17*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u6;
            f_udg[18*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u7;
            f_udg[19*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t62*t71-epsilon0*ff*rampFactor*t3*t47*t48*2.0;
            f_udg[20*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u9;
            f_udg[21*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u10;
            f_udg[22*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t71*u11;
            f_udg[23*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t69*t71-epsilon0*ff*rampFactor*t3*t15*t47*2.0;
            f_udg[24*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u5*(-1.0/2.0);
            f_udg[25*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u6*(-1.0/2.0);
            f_udg[26*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u7*(-1.0/2.0);
            f_udg[27*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*t62*(-1.0/2.0);
            f_udg[28*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u9*(-1.0/2.0);
            f_udg[29*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u10*(-1.0/2.0);
            f_udg[30*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*u11*(-1.0/2.0);
            f_udg[31*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t16*t17*t22*t31*t69*(-1.0/2.0);
            f_udg[32*ng+i] += t75-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u5;
            f_udg[33*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u6;
            f_udg[34*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u7;
            f_udg[35*ng+i] += t77-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*t62*u2;
            f_udg[36*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u9;
            f_udg[37*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u10;
            f_udg[38*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u2*u11;
            f_udg[39*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*t69*u2;
            f_udg[40*ng+i] += t78;
            f_udg[41*ng+i] += t75+t79;
            f_udg[42*ng+i] += t80;
            f_udg[43*ng+i] += t81-epsilon0*ff*rampFactor*t3*t47*u2*2.0;
            f_udg[44*ng+i] += t82;
            f_udg[45*ng+i] += t83;
            f_udg[46*ng+i] += t84;
            f_udg[47*ng+i] += t85;
            f_udg[48*ng+i] += 0.0;
            f_udg[49*ng+i] += 0.0;
            f_udg[50*ng+i] += t75;
            f_udg[51*ng+i] += epsilon0*ff*rampFactor*t3*t47*u3*-2.0;
            f_udg[52*ng+i] += 0.0;
            f_udg[53*ng+i] += 0.0;
            f_udg[54*ng+i] += 0.0;
            f_udg[55*ng+i] += 0.0;
            f_udg[56*ng+i] += 0.0;
            f_udg[57*ng+i] += 0.0;
            f_udg[58*ng+i] += 0.0;
            f_udg[59*ng+i] += t86;
            f_udg[60*ng+i] += 0.0;
            f_udg[61*ng+i] += 0.0;
            f_udg[62*ng+i] += 0.0;
            f_udg[63*ng+i] += 0.0;
            f_udg[64*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u5;
            f_udg[65*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u6;
            f_udg[66*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u7;
            f_udg[67*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*t62*u3;
            f_udg[68*ng+i] += t75-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u9;
            f_udg[69*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u10;
            f_udg[70*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*u3*u11;
            f_udg[71*ng+i] += t77-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t25*t69*u3;
            f_udg[72*ng+i] += 0.0;
            f_udg[73*ng+i] += 0.0;
            f_udg[74*ng+i] += 0.0;
            f_udg[75*ng+i] += 0.0;
            f_udg[76*ng+i] += 0.0;
            f_udg[77*ng+i] += t75;
            f_udg[78*ng+i] += 0.0;
            f_udg[79*ng+i] += epsilon0*ff*rampFactor*t3*t47*u2*-2.0;
            f_udg[80*ng+i] += t78;
            f_udg[81*ng+i] += t79;
            f_udg[82*ng+i] += t80;
            f_udg[83*ng+i] += t81;
            f_udg[84*ng+i] += t82;
            f_udg[85*ng+i] += t83;
            f_udg[86*ng+i] += t75+t84;
            f_udg[87*ng+i] += t85-epsilon0*ff*rampFactor*t3*t47*u3*2.0;
            f_udg[88*ng+i] += 0.0;
            f_udg[89*ng+i] += 0.0;
            f_udg[90*ng+i] += 0.0;
            f_udg[91*ng+i] += 0.0;
            f_udg[92*ng+i] += 0.0;
            f_udg[93*ng+i] += 0.0;
            f_udg[94*ng+i] += 0.0;
            f_udg[95*ng+i] += t86;
        }
    }
}
