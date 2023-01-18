
// Written by: C. Nguyen & P. Fernandez

void flux_ns2dNEW(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
			   double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

	/* Allocate dynamic memory */
	double * mu  = new double [ng];
	double * mu_udg;
	double * f_mu;
	if (computeJacobian == 1 && app.viscosityModel != 0) {
		mu_udg = new double [ng * nc];
		f_mu = new double [ng * ncu * nd];
	}

    getViscosity(mu, mu_udg, pg, udg, param, app.viscosityModel, ng, nc, ncu, nd, computeJacobian);
    
    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];

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

        double t2 = u2*u2;
        double t3 = 1.0/(u1*u1);
        double t4 = 1.0/u1;
        double t5 = 1.0/param3;
        double t6 = t2*t3*(1.0/2.0);
        double t7 = u3*u3;
        double t8 = t3*t7*(1.0/2.0);
        double t9 = t6+t8;
        double t23 = t9*u1;
        double t10 = -t23+u4;
        double t11 = param1-1.0;
        double t22 = t4*u3*u5;
        double t12 = -t22+u7;
        double t13 = t4*t12;
        double t25 = t4*u2*u9;
        double t14 = -t25+u10;
        double t15 = t4*t14;
        double t16 = t13+t15;
        double t21 = t4*u2*u5;
        double t17 = -t21+u6;
        double t18 = t4*t17*(4.0/3.0);
        double t29 = t4*u3*u9;
        double t19 = -t29+u11;
        double t20 = t18-t4*t19*(2.0/3.0);
        double t24 = t4*u2*u3;
        double t26 = mu[i]*t5*t16;
        double t27 = t24+t26;
        double t28 = t10*t11;
        double t30 = t4*u4;
        double t31 = t4*t10*t11;
        double t32 = t30+t31;
        double t33 = t4*t17*(2.0/3.0);
        double t34 = t33-t4*t19*(4.0/3.0);
        double t35 = 1.0/param4;
        double t36 = 1.0/t11;

        f[0*ng+i] = u2;
        f[1*ng+i] = t28+t2*t4+mu[i]*t5*t20;
        f[2*ng+i] = t27;
        f[3*ng+i] = t32*u2+mu[i]*t4*t5*t16*u3+mu[i]*t4*t5*t20*u2-mu[i]*param1*t3*t5*t35*t36*(t11*u1*(-u8+t9*u5+u1*(t3*t12*u3+t3*t17*u2))+t10*t11*u5);
        f[4*ng+i] = u3;
        f[5*ng+i] = t27;
        f[6*ng+i] = t28+t4*t7-mu[i]*t5*t34;
        f[7*ng+i] = t32*u3+mu[i]*t4*t5*t16*u2-mu[i]*t4*t5*t34*u3-mu[i]*param1*t3*t5*t35*t36*(t11*u1*(-u12+t9*u9+u1*(t3*t14*u2+t3*t19*u3))+t10*t11*u9);
    }

    if (computeJacobian == 1) {
        for (int i = 0; i <ng; i++) {
            double x1 = pg[0*ng+i];
            double x2 = pg[1*ng+i];

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

            double t2 = 1.0/(u1*u1);
            double t3 = u2*u2;
            double t4 = 1.0/(u1*u1*u1);
            double t5 = u3*u3;
            double t6 = 1.0/u1;
            double t7 = 1.0/param3;
            double t8 = t2*t3*(1.0/2.0);
            double t9 = t2*t5*(1.0/2.0);
            double t10 = param1-1.0;
            double t11 = t3*t4;
            double t12 = t4*t5;
            double t13 = t11+t12;
            double t39 = t13*u1;
            double t14 = t8+t9-t39;
            double t19 = t6*u3*u5;
            double t15 = -t19+u7;
            double t21 = t6*u2*u9;
            double t16 = -t21+u10;
            double t24 = t6*u2*u5;
            double t17 = -t24+u6;
            double t26 = t6*u3*u9;
            double t18 = -t26+u11;
            double t20 = t2*t15;
            double t22 = t2*t16;
            double t41 = t4*u3*u5;
            double t42 = t4*u2*u9;
            double t23 = t20+t22-t41-t42;
            double t25 = t2*t17*(4.0/3.0);
            double t27 = t4*u3*u9*(2.0/3.0);
            double t28 = t25+t27-t2*t18*(2.0/3.0)-t4*u2*u5*(4.0/3.0);
            double t29 = t8+t9;
            double t45 = t29*u1;
            double t30 = -t45+u4;
            double t31 = 1.0/param4;
            double t32 = 1.0/t10;
            double t33 = t29*u5;
            double t34 = t2*t17*u2;
            double t35 = t2*t15*u3;
            double t36 = t34+t35;
            double t37 = t36*u1;
            double t38 = t33+t37-u8;
            double t40 = 1.0/(u1*u1*u1*u1);
            double t43 = -mu[i]*t7*t23-t2*u2*u3;
            double t44 = t2*u4;
            double t46 = t2*t10*t30;
            double t47 = t6*t10*t14;
            double t48 = t44+t46+t47;
            double t49 = t6*t15;
            double t50 = t6*t16;
            double t51 = t49+t50;
            double t52 = t2*t17*(2.0/3.0);
            double t53 = t4*u3*u9*(4.0/3.0);
            double t54 = t52+t53-t2*t18*(4.0/3.0)-t4*u2*u5*(2.0/3.0);
            double t55 = t29*u9;
            double t56 = t2*t16*u2;
            double t57 = t2*t18*u3;
            double t58 = t56+t57;
            double t59 = t58*u1;
            double t60 = t55+t59-u12;
            double t61 = t6*t17*(4.0/3.0);
            double t62 = t61-t6*t18*(2.0/3.0);
            double t63 = t6*u3;
            double t64 = t63-mu[i]*t2*t7*u9;
            double t65 = mu[i]*t6*t7*t51;
            double t66 = t6*u2;
            double t67 = t66-mu[i]*t2*t7*u5;
            double t68 = t6*u4;
            double t69 = t6*t10*t30;
            double t70 = t6*t17*(2.0/3.0);
            double t71 = t70-t6*t18*(4.0/3.0);
            double t72 = t6*t10;
            double t73 = t6+t72;
            double t74 = mu[i]*t6*t7;
            double t75 = mu[i]*t2*t7*u2;
            double t76 = t10*t30;
            double t77 = t10*t14*u1;
            double t78 = t76+t77;
            double t79 = mu[i]*t2*t7*u3;
            double t80 = mu[i]*t2*t7*u2*(2.0/3.0);
            double t81 = mu[i]*t6*t7*(4.0/3.0);
            double t82 = mu[i]*param1*t6*t7*t31;

            f_udg[0*ng+i] = 0.0;
            f_udg[1*ng+i] = -t2*t3-t10*t14-mu[i]*t7*t28;
            f_udg[2*ng+i] = t43;
            f_udg[3*ng+i] = -t48*u2-mu[i]*t6*t7*t23*u3-mu[i]*t6*t7*t28*u2-mu[i]*t2*t7*t51*u3-mu[i]*t2*t7*t62*u2-mu[i]*param1*t2*t7*t31*t32*(t10*t38-t10*t14*u5+t10*u1*(t34+t35-u1*(t4*t15*u3*2.0+t4*t17*u2*2.0-t3*t40*u5-t5*t40*u5)-t13*u5))+mu[i]*param1*t4*t7*t31*t32*(t10*t30*u5+t10*t38*u1)*2.0;
            f_udg[4*ng+i] = 0.0;
            f_udg[5*ng+i] = t43;
            f_udg[6*ng+i] = -t2*t5-t10*t14+mu[i]*t7*t54;
            f_udg[7*ng+i] = -t48*u3-mu[i]*t6*t7*t23*u2-mu[i]*t2*t7*t51*u2+mu[i]*t6*t7*t54*u3+mu[i]*t2*t7*t71*u3-mu[i]*param1*t2*t7*t31*t32*(t10*t60-t10*t14*u9+t10*u1*(t56+t57-u1*(t4*t16*u2*2.0+t4*t18*u3*2.0-t3*t40*u9-t5*t40*u9)-t13*u9))+mu[i]*param1*t4*t7*t31*t32*(t10*t30*u9+t10*t60*u1)*2.0;
            f_udg[8*ng+i] = 1.0;
            f_udg[9*ng+i] = t6*u2*2.0-t6*t10*u2-mu[i]*t2*t7*u5*(4.0/3.0);
            f_udg[10*ng+i] = t64;
            f_udg[11*ng+i] = t68+t69-t2*t3*t10+mu[i]*t6*t7*t62-mu[i]*t4*t7*u2*u5*(4.0/3.0)-mu[i]*t4*t7*u3*u9-mu[i]*param1*t2*t7*t31*t32*(t10*u1*(u1*(t2*t17-t4*u2*u5)+t2*u2*u5)-t6*t10*u2*u5);
            f_udg[12*ng+i] = 0.0;
            f_udg[13*ng+i] = t64;
            f_udg[14*ng+i] = -t6*t10*u2+mu[i]*t2*t7*u5*(2.0/3.0);
            f_udg[15*ng+i] = t65-t2*t10*u2*u3+mu[i]*t4*t7*u3*u5*(2.0/3.0)-mu[i]*t4*t7*u2*u9-mu[i]*param1*t2*t7*t31*t32*(t10*u1*(u1*(t22-t42)+t2*u2*u9)-t6*t10*u2*u9);
            f_udg[16*ng+i] = 0.0;
            f_udg[17*ng+i] = -t6*t10*u3+mu[i]*t2*t7*u9*(2.0/3.0);
            f_udg[18*ng+i] = t67;
            f_udg[19*ng+i] = t65-t2*t10*u2*u3-mu[i]*t4*t7*u3*u5+mu[i]*t4*t7*u2*u9*(2.0/3.0)-mu[i]*param1*t2*t7*t31*t32*(t10*u1*(u1*(t20-t41)+t2*u3*u5)-t6*t10*u3*u5);
            f_udg[20*ng+i] = 1.0;
            f_udg[21*ng+i] = t67;
            f_udg[22*ng+i] = t6*u3*2.0-t6*t10*u3-mu[i]*t2*t7*u9*(4.0/3.0);
            f_udg[23*ng+i] = t68+t69-t2*t5*t10-mu[i]*t6*t7*t71-mu[i]*t4*t7*u2*u5-mu[i]*t4*t7*u3*u9*(4.0/3.0)-mu[i]*param1*t2*t7*t31*t32*(t10*u1*(u1*(t2*t18-t4*u3*u9)+t2*u3*u9)-t6*t10*u3*u9);
            f_udg[24*ng+i] = 0.0;
            f_udg[25*ng+i] = t10;
            f_udg[26*ng+i] = 0.0;
            f_udg[27*ng+i] = t73*u2-mu[i]*param1*t2*t7*t31*u5;
            f_udg[28*ng+i] = 0.0;
            f_udg[29*ng+i] = 0.0;
            f_udg[30*ng+i] = t10;
            f_udg[31*ng+i] = t73*u3-mu[i]*param1*t2*t7*t31*u9;
            f_udg[32*ng+i] = 0.0;
            f_udg[33*ng+i] = mu[i]*t2*t7*u2*(-4.0/3.0);
            f_udg[34*ng+i] = -mu[i]*t2*t7*u3;
            f_udg[35*ng+i] = mu[i]*t3*t4*t7*(-4.0/3.0)-mu[i]*t4*t5*t7-mu[i]*param1*t2*t7*t31*t32*t78;
            f_udg[36*ng+i] = 0.0;
            f_udg[37*ng+i] = -mu[i]*t2*t7*u3;
            f_udg[38*ng+i] = t80;
            f_udg[39*ng+i] = mu[i]*t4*t7*u2*u3*(-1.0/3.0);
            f_udg[40*ng+i] = 0.0;
            f_udg[41*ng+i] = t81;
            f_udg[42*ng+i] = 0.0;
            f_udg[43*ng+i] = mu[i]*t2*t7*u2*(4.0/3.0)-mu[i]*param1*t2*t7*t31*u2;
            f_udg[44*ng+i] = 0.0;
            f_udg[45*ng+i] = 0.0;
            f_udg[46*ng+i] = mu[i]*t6*t7*(-2.0/3.0);
            f_udg[47*ng+i] = mu[i]*t2*t7*u3*(-2.0/3.0);
            f_udg[48*ng+i] = 0.0;
            f_udg[49*ng+i] = 0.0;
            f_udg[50*ng+i] = t74;
            f_udg[51*ng+i] = t79-mu[i]*param1*t2*t7*t31*u3;
            f_udg[52*ng+i] = 0.0;
            f_udg[53*ng+i] = t74;
            f_udg[54*ng+i] = 0.0;
            f_udg[55*ng+i] = t75;
            f_udg[56*ng+i] = 0.0;
            f_udg[57*ng+i] = 0.0;
            f_udg[58*ng+i] = 0.0;
            f_udg[59*ng+i] = t82;
            f_udg[60*ng+i] = 0.0;
            f_udg[61*ng+i] = 0.0;
            f_udg[62*ng+i] = 0.0;
            f_udg[63*ng+i] = 0.0;
            f_udg[64*ng+i] = 0.0;
            f_udg[65*ng+i] = mu[i]*t2*t7*u3*(2.0/3.0);
            f_udg[66*ng+i] = -t75;
            f_udg[67*ng+i] = mu[i]*t4*t7*u2*u3*(-1.0/3.0);
            f_udg[68*ng+i] = 0.0;
            f_udg[69*ng+i] = -t75;
            f_udg[70*ng+i] = mu[i]*t2*t7*u3*(-4.0/3.0);
            f_udg[71*ng+i] = -mu[i]*t3*t4*t7-mu[i]*t4*t5*t7*(4.0/3.0)-mu[i]*param1*t2*t7*t31*t32*t78;
            f_udg[72*ng+i] = 0.0;
            f_udg[73*ng+i] = 0.0;
            f_udg[74*ng+i] = t74;
            f_udg[75*ng+i] = t79;
            f_udg[76*ng+i] = 0.0;
            f_udg[77*ng+i] = t74;
            f_udg[78*ng+i] = 0.0;
            f_udg[79*ng+i] = t75-mu[i]*param1*t2*t7*t31*u2;
            f_udg[80*ng+i] = 0.0;
            f_udg[81*ng+i] = mu[i]*t6*t7*(-2.0/3.0);
            f_udg[82*ng+i] = 0.0;
            f_udg[83*ng+i] = -t80;
            f_udg[84*ng+i] = 0.0;
            f_udg[85*ng+i] = 0.0;
            f_udg[86*ng+i] = t81;
            f_udg[87*ng+i] = mu[i]*t2*t7*u3*(4.0/3.0)-mu[i]*param1*t2*t7*t31*u3;
            f_udg[88*ng+i] = 0.0;
            f_udg[89*ng+i] = 0.0;
            f_udg[90*ng+i] = 0.0;
            f_udg[91*ng+i] = 0.0;
            f_udg[92*ng+i] = 0.0;
            f_udg[93*ng+i] = 0.0;
            f_udg[94*ng+i] = 0.0;
            f_udg[95*ng+i] = t82;
        }

        if (app.viscosityModel != 0) {

            for (int i = 0; i <ng; i++) {
                double x1 = pg[0*ng+i];
                double x2 = pg[1*ng+i];

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

                double t2 = 1.0/u1;
                double t3 = 1.0/param3;
                double t15 = t2*u3*u5;
                double t4 = -t15+u7;
                double t5 = t2*t4;
                double t22 = t2*u2*u9;
                double t6 = -t22+u10;
                double t7 = t2*t6;
                double t8 = t5+t7;
                double t14 = t2*u2*u5;
                double t9 = -t14+u6;
                double t10 = t2*t9*(4.0/3.0);
                double t24 = t2*u3*u9;
                double t11 = -t24+u11;
                double t12 = t10-t2*t11*(2.0/3.0);
                double t13 = 1.0/(u1*u1);
                double t16 = u2*u2;
                double t17 = t13*t16*(1.0/2.0);
                double t18 = u3*u3;
                double t19 = t13*t18*(1.0/2.0);
                double t20 = t17+t19;
                double t21 = param1-1.0;
                double t23 = t3*t8;
                double t25 = t2*t9*(2.0/3.0);
                double t26 = t25-t2*t11*(4.0/3.0);
                double t27 = 1.0/param4;
                double t28 = u4-t20*u1;
                double t29 = 1.0/t21;

                f_mu[0*ng+i] = 0.0;
                f_mu[1*ng+i] = t3*t12;
                f_mu[2*ng+i] = t23;
                f_mu[3*ng+i] = t2*t3*t8*u3+t2*t3*t12*u2-param1*t3*t13*t27*t29*(t21*u1*(-u8+t20*u5+u1*(t4*t13*u3+t9*t13*u2))+t21*t28*u5);
                f_mu[4*ng+i] = 0.0;
                f_mu[5*ng+i] = t23;
                f_mu[6*ng+i] = -t3*t26;
                f_mu[7*ng+i] = t2*t3*t8*u2-t2*t3*t26*u3-param1*t3*t13*t27*t29*(t21*u1*(-u12+t20*u9+u1*(t6*t13*u2+t11*t13*u3))+t21*t28*u9);
            }

            chainRule_f_udg(f_udg, f_mu, mu_udg, ng, nc, ncu, nd);
        }
    }
    
	/* Deallocate dynamic memory */
	delete[] mu;
	if (computeJacobian == 1 && app.viscosityModel != 0) {
		delete[] mu_udg;
		delete[] f_mu;
	}
}
