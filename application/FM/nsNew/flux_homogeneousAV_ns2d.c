
// Written by: C. Nguyen & P. Fernandez

void flux_homogeneousAV_ns2d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
                             double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    double r, ru, rv, rE, rx, ry, rux, rvy, u, v, ux, vy, p;
    double p_reg, D_p_reg_D_p, r_reg, D_r_reg_D_r, trash;
    double c, x, l, wallDistance, ff, D_ff_D_x;

    double param1 = param[0];
    double param2 = param[1];
    double param3 = param[2];
    double param4 = param[3];
    double param5 = param[4];
    double param6 = param[5];

    double alpha = 1.0e4;
    double beta = 1.0e-2;
    double alpha_ff = 10.0;
    double beta_ff = 0.5;
    double eps_v = 1.0e-8;
    double rRegMin = 1.0e-2;
    double pRegMin = 1.0e-3;
    double rampFactor = app.rampFactor;
    double epsilon0 = 1.0;
    //double h0 = 1.0e-3;
    double h0;

    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];

        h0 = pg[2*ng+i];
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
        x = (l*(ux+vy))/c;
        ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;

//        printf("p = %g, p_reg = %g, D_p_reg_D_p = %g, r = %g, r_reg = %g, D_p_reg_D_p = %g, x-beta_ff = %g, ff = %g, D_ff_Dx = %g\n",p, p_reg, D_p_reg_D_p, r, r_reg, D_r_reg_D_r, x-beta_ff, ff, D_ff_D_x);

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

            h0 = pg[2*ng+i];
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
            x = (l*(ux+vy))/c;
            ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;
            D_ff_D_x = exp(alpha_ff*(x-beta_ff))/(1 + exp(alpha_ff*(x-beta_ff)));

            double t2 = 1.0/(u1*u1);
            double t3 = 1.0/u1;
            double t4 = 1.0/(u1*u1*u1);
            double t9 = t3*u2*u5;
            double t5 = -t9+u6;
            double t11 = t3*u3*u9;
            double t6 = -t11+u11;
            double t7 = 1.0/r_reg;
            double t8 = p_reg*param1*t7;
            double t10 = t3*t5;
            double t12 = t3*t6;
            double t13 = t10+t12;
            double t14 = 1.0/pow(t8,3.0/2.0);
            double t15 = u2*u2;
            double t16 = u3*u3;
            double t17 = 1.0/sqrt(t8);
            double t18 = t2*t5;
            double t19 = t2*t6;
            double t34 = t4*u2*u5;
            double t35 = t4*u3*u9;
            double t20 = t18+t19-t34-t35;
            double t21 = 1.0/(r_reg*r_reg);
            double t22 = D_r_reg_D_r*l*p_reg*param1*t13*t14*t21*(1.0/2.0);
            double t23 = param1-1.0;
            double t24 = t2*t15*(1.0/2.0);
            double t25 = t2*t16*(1.0/2.0);
            double t26 = t4*t15*2.0;
            double t27 = t4*t16*2.0;
            double t28 = t26+t27;
            double t29 = eps_v*eps_v;
            double t30 = t29*(1.0/2.0);
            double t37 = t28*u1*(1.0/2.0);
            double t31 = t24+t25+t30-t37;
            double t32 = D_p_reg_D_p*l*param1*t7*t13*t14*t23*t31*(1.0/2.0);
            double t36 = l*t17*t20;
            double t33 = t22+t32-t36;
            double t38 = param1*(1.0/2.0);
            double t39 = t38-1.0/2.0;
            double t50 = t3*u3*u5;
            double t40 = -t50+u7;
            double t41 = t2*t15;
            double t42 = t2*t16;
            double t43 = t29+t41+t42;
            double t56 = t3*u2*u9;
            double t44 = -t56+u10;
            double t45 = l*t2*t17*u5;
            double t47 = D_p_reg_D_p*l*param1*t3*t7*t13*t14*t23*u2*(1.0/2.0);
            double t46 = t45-t47;
            double t48 = t43*u5;
            double t49 = t3*t5*u2*2.0;
            double t51 = t3*t40*u3*2.0;
            double t52 = t48+t49+t51;
            double t53 = param1*u8;
            double t65 = t39*t52;
            double t54 = t53-t65;
            double t55 = t43*u9;
            double t57 = t3*t44*u2*2.0;
            double t58 = t3*t6*u3*2.0;
            double t59 = t55+t57+t58;
            double t60 = param1*u12;
            double t66 = t39*t59;
            double t61 = t60-t66;
            double t62 = l*t2*t17*u9;
            double t64 = D_p_reg_D_p*l*param1*t3*t7*t13*t14*t23*u3*(1.0/2.0);
            double t63 = t62-t64;
            double t67 = epsilon0*ff*rampFactor;
            double t68 = -t29+t41+t42;
            double t69 = epsilon0*ff*rampFactor*t39*t68;
            double t70 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u5;
            double t71 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u6;
            double t72 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u7;
            double t73 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*(t53-t65);
            double t74 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u9;
            double t75 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u10;
            double t76 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*u11;
            double t77 = D_ff_D_x*epsilon0*l*rampFactor*t3*t17*(t60-t66);
            double t78 = epsilon0*ff*param1*rampFactor;

            f_udg[0*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u5;
            f_udg[1*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u6;
            f_udg[2*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u7;
            f_udg[3*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*t54+epsilon0*ff*rampFactor*t39*(t28*u5+t2*t5*u2*2.0-t4*t15*u5*2.0-t4*t16*u5*2.0+t2*t40*u3*2.0);
            f_udg[4*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u9;
            f_udg[5*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u10;
            f_udg[6*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*u11;
            f_udg[7*ng+i] += D_ff_D_x*epsilon0*rampFactor*t33*t61+epsilon0*ff*rampFactor*t39*(t28*u9+t2*t6*u3*2.0-t4*t15*u9*2.0-t4*t16*u9*2.0+t2*t44*u2*2.0);
            f_udg[8*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u5;
            f_udg[9*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u6;
            f_udg[10*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u7;
            f_udg[11*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*t54-epsilon0*ff*rampFactor*t3*t5*t39*2.0;
            f_udg[12*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u9;
            f_udg[13*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u10;
            f_udg[14*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*u11;
            f_udg[15*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t46*t61-epsilon0*ff*rampFactor*t3*t39*t44*2.0;
            f_udg[16*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u5;
            f_udg[17*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u6;
            f_udg[18*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u7;
            f_udg[19*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t54*t63-epsilon0*ff*rampFactor*t3*t39*t40*2.0;
            f_udg[20*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u9;
            f_udg[21*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u10;
            f_udg[22*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t63*u11;
            f_udg[23*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t61*t63-epsilon0*ff*rampFactor*t3*t6*t39*2.0;
            f_udg[24*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u5*(-1.0/2.0);
            f_udg[25*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u6*(-1.0/2.0);
            f_udg[26*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u7*(-1.0/2.0);
            f_udg[27*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*t54*(-1.0/2.0);
            f_udg[28*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u9*(-1.0/2.0);
            f_udg[29*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u10*(-1.0/2.0);
            f_udg[30*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*u11*(-1.0/2.0);
            f_udg[31*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*rampFactor*t7*t13*t14*t23*t61*(-1.0/2.0);
            f_udg[32*ng+i] += t67-D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u5;
            f_udg[33*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u6;
            f_udg[34*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u7;
            f_udg[35*ng+i] += t69-D_ff_D_x*epsilon0*l*rampFactor*t2*t17*t54*u2;
            f_udg[36*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u9;
            f_udg[37*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u10;
            f_udg[38*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u2*u11;
            f_udg[39*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*t61*u2;
            f_udg[40*ng+i] += t70;
            f_udg[41*ng+i] += t67+t71;
            f_udg[42*ng+i] += t72;
            f_udg[43*ng+i] += t73-epsilon0*ff*rampFactor*t3*t39*u2*2.0;
            f_udg[44*ng+i] += t74;
            f_udg[45*ng+i] += t75;
            f_udg[46*ng+i] += t76;
            f_udg[47*ng+i] += t77;
            f_udg[48*ng+i] += 0.0;
            f_udg[49*ng+i] += 0.0;
            f_udg[50*ng+i] += t67;
            f_udg[51*ng+i] += epsilon0*ff*rampFactor*t3*t39*u3*-2.0;
            f_udg[52*ng+i] += 0.0;
            f_udg[53*ng+i] += 0.0;
            f_udg[54*ng+i] += 0.0;
            f_udg[55*ng+i] += 0.0;
            f_udg[56*ng+i] += 0.0;
            f_udg[57*ng+i] += 0.0;
            f_udg[58*ng+i] += 0.0;
            f_udg[59*ng+i] += t78;
            f_udg[60*ng+i] += 0.0;
            f_udg[61*ng+i] += 0.0;
            f_udg[62*ng+i] += 0.0;
            f_udg[63*ng+i] += 0.0;
            f_udg[64*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u5;
            f_udg[65*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u6;
            f_udg[66*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u7;
            f_udg[67*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*t54*u3;
            f_udg[68*ng+i] += t67-D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u9;
            f_udg[69*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u10;
            f_udg[70*ng+i] += -D_ff_D_x*epsilon0*l*rampFactor*t2*t17*u3*u11;
            f_udg[71*ng+i] += t69-D_ff_D_x*epsilon0*l*rampFactor*t2*t17*t61*u3;
            f_udg[72*ng+i] += 0.0;
            f_udg[73*ng+i] += 0.0;
            f_udg[74*ng+i] += 0.0;
            f_udg[75*ng+i] += 0.0;
            f_udg[76*ng+i] += 0.0;
            f_udg[77*ng+i] += t67;
            f_udg[78*ng+i] += 0.0;
            f_udg[79*ng+i] += epsilon0*ff*rampFactor*t3*t39*u2*-2.0;
            f_udg[80*ng+i] += t70;
            f_udg[81*ng+i] += t71;
            f_udg[82*ng+i] += t72;
            f_udg[83*ng+i] += t73;
            f_udg[84*ng+i] += t74;
            f_udg[85*ng+i] += t75;
            f_udg[86*ng+i] += t67+t76;
            f_udg[87*ng+i] += t77-epsilon0*ff*rampFactor*t3*t39*u3*2.0;
            f_udg[88*ng+i] += 0.0;
            f_udg[89*ng+i] += 0.0;
            f_udg[90*ng+i] += 0.0;
            f_udg[91*ng+i] += 0.0;
            f_udg[92*ng+i] += 0.0;
            f_udg[93*ng+i] += 0.0;
            f_udg[94*ng+i] += 0.0;
            f_udg[95*ng+i] += t78;
        }
    }
}
