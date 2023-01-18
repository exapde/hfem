
#include "../../../utilities/surrogateMaximum.cpp"

// Written by: C. Nguyen & P. Fernandez

void flux_isotropicAV_ns2d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
                           double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
	double r, ru, rv, rE, rx, ry, rux, rvy, u, v, ux, vy, p, H, rH, h;
	double p_reg, D_p_reg_D_p, r_reg, H_reg, D_H_reg_D_H, D_r_reg_D_r, trash;
	double c_cr, x, ff, D_ff_D_x;

	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

    double alpha = 1.0e4;
    double beta = 1.0e-2;
    double k_h = 1.5;
    double eps_v = 1.0e-8;
    double rRegMin = 1.0e-2;
    double pRegMin = 1.0e-3;
    double HregMin = 1.0e-2;
	double rampFactor = app.rampFactor;
    int porder = app.porder[0];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];

		h = pg[2*ng+i];
		error("Need to determine the initial index of element size (h)\n");

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
        rH = rE + p;
        H = rH/r;
        surrogateMaximum(&H_reg, &trash, &D_H_reg_D_H, HregMin, H, alpha);

        c_cr = sqrt((2*(param1-1)*H_reg)/(param1+1));
        x = ((k_h*h/porder)*(ux+vy))/c_cr;
        surrogateMaximum(&ff, &trash, &D_ff_D_x, 0, x-beta, alpha);

		double t2 = 1.0/(u1*u1);
		double t3 = 1.0/porder;
		double t4 = 1.0/r_reg;
		double t5 = p_reg*param1*t4;
		double t6 = sqrt(t5);
		double t7 = u2*u2;
		double t8 = t2*t7;
		double t9 = u3*u3;
		double t10 = t2*t9;
		double t11 = eps_v*eps_v;
		double t12 = t8+t10+t11;
		double t13 = sqrt(t12);
		double t14 = t6+t13;
		double t15 = 1.0/u1;
		double t16 = param1*(1.0/2.0);
		double t17 = t16-1.0/2.0;

		f[0*ng+i] += ff*h*k_h*rampFactor*t3*t14*u5;
		f[1*ng+i] += ff*h*k_h*rampFactor*t3*t14*u6;
		f[2*ng+i] += ff*h*k_h*rampFactor*t3*t14*u7;
		f[3*ng+i] += -ff*h*k_h*rampFactor*t3*t14*(t17*(t12*u5+t15*u2*(u6-t15*u2*u5)*2.0+t15*u3*(u7-t15*u3*u5)*2.0)-param1*u8);
		f[4*ng+i] += ff*h*k_h*rampFactor*t3*t14*u9;
		f[5*ng+i] += ff*h*k_h*rampFactor*t3*t14*u10;
		f[6*ng+i] += ff*h*k_h*rampFactor*t3*t14*u11;
		f[7*ng+i] += -ff*h*k_h*rampFactor*t3*t14*(t17*(t12*u9+t15*u2*(u10-t15*u2*u9)*2.0+t15*u3*(u11-t15*u3*u9)*2.0)-param1*u12);
	}

	if (computeJacobian == 1) {

		for (int i = 0; i <ng; i++) {
			double x1 = pg[0*ng+i];
			double x2 = pg[1*ng+i];

			h = pg[2*ng+i];
			error("Need to determine the initial index of element size (h)\n");

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
            rH = rE + p;
            H = rH/r;
            surrogateMaximum(&H_reg, &trash, &D_H_reg_D_H, HregMin, H, alpha);

            c_cr = sqrt((2*(param1-1)*H_reg)/(param1+1));
            x = ((k_h*h/porder)*(ux+vy))/c_cr;
            surrogateMaximum(&ff, &trash, &D_ff_D_x, 0, x-beta, alpha);

            double t2 = 1.0/(u1*u1);
            double t3 = 1.0/porder;
            double t4 = 1.0/u1;
            double t5 = 1.0/(u1*u1*u1);
            double t6 = param1*2.0;
            double t7 = t6-2.0;
            double t8 = u2*u2;
            double t9 = t2*t8;
            double t10 = u3*u3;
            double t11 = t2*t10;
            double t12 = eps_v*eps_v;
            double t13 = t9+t11+t12;
            double t14 = param1-1.0;
            double t34 = t4*u2*u5;
            double t15 = -t34+u6;
            double t36 = t4*u3*u9;
            double t16 = -t36+u11;
            double t17 = param1+1.0;
            double t18 = 1.0/t17;
            double t19 = H_reg*t7*t18;
            double t20 = t5*t8*2.0;
            double t21 = t5*t10*2.0;
            double t22 = t20+t21;
            double t23 = 1.0/r_reg;
            double t24 = p_reg*param1*t23;
            double t25 = 1.0/sqrt(t24);
            double t26 = t2*t8*(1.0/2.0);
            double t27 = t2*t10*(1.0/2.0);
            double t28 = t12*(1.0/2.0);
            double t44 = t22*u1*(1.0/2.0);
            double t29 = t26+t27+t28-t44;
            double t30 = sqrt(t24);
            double t31 = sqrt(t13);
            double t32 = t30+t31;
            double t33 = 1.0/sqrt(t19);
            double t35 = t2*t15;
            double t37 = t2*t16;
            double t54 = t5*u2*u5;
            double t55 = t5*u3*u9;
            double t38 = t35+t37-t54-t55;
            double t39 = h*k_h*t3*t33*t38;
            double t56 = t13*u1*(1.0/2.0);
            double t40 = -t56+u4;
            double t41 = t14*t40;
            double t42 = t41+u4;
            double t43 = t2*t42;
            double t45 = t4*t14*t29;
            double t46 = t43+t45;
            double t47 = t4*t15;
            double t48 = t4*t16;
            double t49 = t47+t48;
            double t50 = 1.0/pow(t19,3.0/2.0);
            double t57 = D_H_reg_D_H*h*k_h*t3*t7*t18*t46*t49*t50*(1.0/2.0);
            double t51 = t39-t57;
            double t52 = 1.0/sqrt(t13);
            double t53 = 1.0/(r_reg*r_reg);
            double t58 = param1*(1.0/2.0);
            double t59 = t58-1.0/2.0;
            double t63 = t4*u3*u5;
            double t60 = -t63+u7;
            double t61 = t13*u5;
            double t62 = t4*t15*u2*2.0;
            double t64 = t4*t60*u3*2.0;
            double t65 = t61+t62+t64;
            double t66 = param1*u8;
            double t68 = t59*t65;
            double t67 = t66-t68;
            double t71 = t4*u2*u9;
            double t69 = -t71+u10;
            double t70 = t13*u9;
            double t72 = t4*t69*u2*2.0;
            double t73 = t4*t16*u3*2.0;
            double t74 = t70+t72+t73;
            double t75 = param1*u12;
            double t77 = t59*t74;
            double t76 = t75-t77;
            double t78 = h*k_h*t2*t3*t33*u5;
            double t80 = D_H_reg_D_H*h*k_h*t2*t3*t7*t14*t18*t49*t50*u2*(1.0/2.0);
            double t79 = t78-t80;
            double t81 = h*k_h*t2*t3*t33*u9;
            double t83 = D_H_reg_D_H*h*k_h*t2*t3*t7*t14*t18*t49*t50*u3*(1.0/2.0);
            double t82 = t81-t83;
            double t84 = h*h;
            double t85 = k_h*k_h;
            double t86 = 1.0/(porder*porder);
            double t87 = ff*h*k_h*rampFactor*t3*t32;
            double t88 = t9+t11-t12;
            double t89 = ff*h*k_h*rampFactor*t3*t32*t59*t88;
            double t90 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u5;
            double t91 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u6;
            double t92 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u7;
            double t93 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*(t66-t68);
            double t94 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u9;
            double t95 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u10;
            double t96 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*u11;
            double t97 = D_ff_D_x*rampFactor*t4*t32*t33*t84*t85*t86*(t75-t77);
            double t98 = ff*h*k_h*param1*rampFactor*t3*t32;

            f_udg[0*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u5-ff*h*k_h*rampFactor*t3*t22*t52*u5*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u5*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u5*(1.0/2.0);
            f_udg[1*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u6-ff*h*k_h*rampFactor*t3*t22*t52*u6*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u6*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u6*(1.0/2.0);
            f_udg[2*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u7-ff*h*k_h*rampFactor*t3*t22*t52*u7*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u7*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u7*(1.0/2.0);
            f_udg[3*ng+i] += ff*h*k_h*rampFactor*t3*t32*t59*(t22*u5-t5*t8*u5*2.0+t2*t15*u2*2.0-t5*t10*u5*2.0+t2*t60*u3*2.0)-D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*t67-ff*h*k_h*rampFactor*t3*t22*t52*t67*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*t67*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*t67*(1.0/2.0);
            f_udg[4*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u9-ff*h*k_h*rampFactor*t3*t22*t52*u9*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u9*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u9*(1.0/2.0);
            f_udg[5*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u10-ff*h*k_h*rampFactor*t3*t22*t52*u10*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u10*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u10*(1.0/2.0);
            f_udg[6*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*u11-ff*h*k_h*rampFactor*t3*t22*t52*u11*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*u11*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*u11*(1.0/2.0);
            f_udg[7*ng+i] += ff*h*k_h*rampFactor*t3*t32*t59*(t22*u9+t2*t16*u3*2.0-t5*t8*u9*2.0-t5*t10*u9*2.0+t2*t69*u2*2.0)-D_ff_D_x*h*k_h*rampFactor*t3*t32*t51*t76-ff*h*k_h*rampFactor*t3*t22*t52*t76*(1.0/2.0)-D_r_reg_D_r*ff*h*k_h*p_reg*param1*rampFactor*t3*t25*t53*t76*(1.0/2.0)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*t29*t76*(1.0/2.0);
            f_udg[8*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u5+ff*h*k_h*rampFactor*t2*t3*t52*u2*u5-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u5*(1.0/2.0);
            f_udg[9*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u6+ff*h*k_h*rampFactor*t2*t3*t52*u2*u6-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u6*(1.0/2.0);
            f_udg[10*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u7+ff*h*k_h*rampFactor*t2*t3*t52*u2*u7-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u7*(1.0/2.0);
            f_udg[11*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t67*t79-ff*h*k_h*rampFactor*t3*t4*t15*t32*t59*2.0+ff*h*k_h*rampFactor*t2*t3*t52*u2*(t66-t68)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*t67*u2*(1.0/2.0);
            f_udg[12*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u9+ff*h*k_h*rampFactor*t2*t3*t52*u2*u9-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u9*(1.0/2.0);
            f_udg[13*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u10+ff*h*k_h*rampFactor*t2*t3*t52*u2*u10-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u10*(1.0/2.0);
            f_udg[14*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t79*u11+ff*h*k_h*rampFactor*t2*t3*t52*u2*u11-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u2*u11*(1.0/2.0);
            f_udg[15*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t76*t79-ff*h*k_h*rampFactor*t3*t4*t32*t59*t69*2.0+ff*h*k_h*rampFactor*t2*t3*t52*u2*(t75-t77)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*t76*u2*(1.0/2.0);
            f_udg[16*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u5+ff*h*k_h*rampFactor*t2*t3*t52*u3*u5-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u5*(1.0/2.0);
            f_udg[17*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u6+ff*h*k_h*rampFactor*t2*t3*t52*u3*u6-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u6*(1.0/2.0);
            f_udg[18*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u7+ff*h*k_h*rampFactor*t2*t3*t52*u3*u7-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u7*(1.0/2.0);
            f_udg[19*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t67*t82-ff*h*k_h*rampFactor*t3*t4*t32*t59*t60*2.0+ff*h*k_h*rampFactor*t2*t3*t52*u3*(t66-t68)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*t67*u3*(1.0/2.0);
            f_udg[20*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u9+ff*h*k_h*rampFactor*t2*t3*t52*u3*u9-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u9*(1.0/2.0);
            f_udg[21*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u10+ff*h*k_h*rampFactor*t2*t3*t52*u3*u10-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u10*(1.0/2.0);
            f_udg[22*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t82*u11+ff*h*k_h*rampFactor*t2*t3*t52*u3*u11-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*u3*u11*(1.0/2.0);
            f_udg[23*ng+i] += -D_ff_D_x*h*k_h*rampFactor*t3*t32*t76*t82-ff*h*k_h*rampFactor*t3*t4*t16*t32*t59*2.0+ff*h*k_h*rampFactor*t2*t3*t52*u3*(t75-t77)-D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t4*t14*t23*t25*t76*u3*(1.0/2.0);
            f_udg[24*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u5*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u5*(1.0/2.0);
            f_udg[25*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u6*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u6*(1.0/2.0);
            f_udg[26*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u7*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u7*(1.0/2.0);
            f_udg[27*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*(t66-t68)*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t67*t84*t85*t86*(1.0/2.0);
            f_udg[28*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u9*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u9*(1.0/2.0);
            f_udg[29*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u10*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u10*(1.0/2.0);
            f_udg[30*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*u11*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t84*t85*t86*u11*(1.0/2.0);
            f_udg[31*ng+i] += D_p_reg_D_p*ff*h*k_h*param1*rampFactor*t3*t14*t23*t25*(t75-t77)*(1.0/2.0)-D_H_reg_D_H*D_ff_D_x*param1*rampFactor*t4*t7*t18*t32*t49*t50*t76*t84*t85*t86*(1.0/2.0);
            f_udg[32*ng+i] += t87-D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u5;
            f_udg[33*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u6;
            f_udg[34*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u7;
            f_udg[35*ng+i] += t89-D_ff_D_x*rampFactor*t2*t32*t33*t67*t84*t85*t86*u2;
            f_udg[36*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u9;
            f_udg[37*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u10;
            f_udg[38*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u2*u11;
            f_udg[39*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t76*t84*t85*t86*u2;
            f_udg[40*ng+i] += t90;
            f_udg[41*ng+i] += t87+t91;
            f_udg[42*ng+i] += t92;
            f_udg[43*ng+i] += t93-ff*h*k_h*rampFactor*t3*t4*t32*t59*u2*2.0;
            f_udg[44*ng+i] += t94;
            f_udg[45*ng+i] += t95;
            f_udg[46*ng+i] += t96;
            f_udg[47*ng+i] += t97;
            f_udg[48*ng+i] += 0.0;
            f_udg[49*ng+i] += 0.0;
            f_udg[50*ng+i] += t87;
            f_udg[51*ng+i] += ff*h*k_h*rampFactor*t3*t4*t32*t59*u3*-2.0;
            f_udg[52*ng+i] += 0.0;
            f_udg[53*ng+i] += 0.0;
            f_udg[54*ng+i] += 0.0;
            f_udg[55*ng+i] += 0.0;
            f_udg[56*ng+i] += 0.0;
            f_udg[57*ng+i] += 0.0;
            f_udg[58*ng+i] += 0.0;
            f_udg[59*ng+i] += t98;
            f_udg[60*ng+i] += 0.0;
            f_udg[61*ng+i] += 0.0;
            f_udg[62*ng+i] += 0.0;
            f_udg[63*ng+i] += 0.0;
            f_udg[64*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u5;
            f_udg[65*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u6;
            f_udg[66*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u7;
            f_udg[67*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t67*t84*t85*t86*u3;
            f_udg[68*ng+i] += t87-D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u9;
            f_udg[69*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u10;
            f_udg[70*ng+i] += -D_ff_D_x*rampFactor*t2*t32*t33*t84*t85*t86*u3*u11;
            f_udg[71*ng+i] += t89-D_ff_D_x*rampFactor*t2*t32*t33*t76*t84*t85*t86*u3;
            f_udg[72*ng+i] += 0.0;
            f_udg[73*ng+i] += 0.0;
            f_udg[74*ng+i] += 0.0;
            f_udg[75*ng+i] += 0.0;
            f_udg[76*ng+i] += 0.0;
            f_udg[77*ng+i] += t87;
            f_udg[78*ng+i] += 0.0;
            f_udg[79*ng+i] += ff*h*k_h*rampFactor*t3*t4*t32*t59*u2*-2.0;
            f_udg[80*ng+i] += t90;
            f_udg[81*ng+i] += t91;
            f_udg[82*ng+i] += t92;
            f_udg[83*ng+i] += t93;
            f_udg[84*ng+i] += t94;
            f_udg[85*ng+i] += t95;
            f_udg[86*ng+i] += t87+t96;
            f_udg[87*ng+i] += t97-ff*h*k_h*rampFactor*t3*t4*t32*t59*u3*2.0;
            f_udg[88*ng+i] += 0.0;
            f_udg[89*ng+i] += 0.0;
            f_udg[90*ng+i] += 0.0;
            f_udg[91*ng+i] += 0.0;
            f_udg[92*ng+i] += 0.0;
            f_udg[93*ng+i] += 0.0;
            f_udg[94*ng+i] += 0.0;
            f_udg[95*ng+i] += t98;
		}
	}
}
