
// Written by: C. Nguyen & P. Fernandez

void flux_euler3dNEW(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
                     double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];

		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];

		double t2 = u2*u2;
		double t3 = 1.0/(u1*u1);
		double t4 = 1.0/u1;
		double t5 = t2*t3*(1.0/2.0);
		double t6 = u3*u3;
		double t7 = t3*t6*(1.0/2.0);
		double t8 = u4*u4;
		double t9 = t3*t8*(1.0/2.0);
		double t10 = t5+t7+t9;
		double t14 = t10*u1;
		double t11 = -t14+u5;
		double t12 = param1-1.0;
		double t13 = t4*u2*u3;
		double t15 = t11*t12;
		double t16 = t4*u5;
		double t17 = t4*t11*t12;
		double t18 = t16+t17;
		double t19 = t4*u2*u4;
		double t20 = t4*u3*u4;

		f[0*ng+i] = u2;
		f[1*ng+i] = t15+t2*t4;
		f[2*ng+i] = t13;
		f[3*ng+i] = t19;
		f[4*ng+i] = t18*u2;
		f[5*ng+i] = u3;
		f[6*ng+i] = t13;
		f[7*ng+i] = t15+t4*t6;
		f[8*ng+i] = t20;
		f[9*ng+i] = t18*u3;
		f[10*ng+i] = u4;
		f[11*ng+i] = t19;
		f[12*ng+i] = t20;
		f[13*ng+i] = t15+t4*t8;
		f[14*ng+i] = t18*u4;
	}

	if (computeJacobian == 1) {

		for (int i = 0; i <ng; i++) {
			double x1 = pg[0*ng+i];
			double x2 = pg[1*ng+i];
			double x3 = pg[2*ng+i];

			double u1 = udg[0*ng+i];
			double u2 = udg[1*ng+i];
			double u3 = udg[2*ng+i];
			double u4 = udg[3*ng+i];
			double u5 = udg[4*ng+i];

			double t2 = 1.0/(u1*u1);
			double t3 = u2*u2;
			double t4 = 1.0/(u1*u1*u1);
			double t5 = u3*u3;
			double t6 = u4*u4;
			double t7 = t2*t3*(1.0/2.0);
			double t8 = t2*t5*(1.0/2.0);
			double t9 = t2*t6*(1.0/2.0);
			double t10 = param1-1.0;
			double t11 = t3*t4;
			double t12 = t4*t5;
			double t13 = t4*t6;
			double t14 = t11+t12+t13;
			double t16 = t14*u1;
			double t15 = t7+t8+t9-t16;
			double t17 = t2*u5;
			double t18 = t7+t8+t9;
			double t24 = t18*u1;
			double t19 = -t24+u5;
			double t20 = t2*t10*t19;
			double t21 = 1.0/u1;
			double t22 = t10*t15*t21;
			double t23 = t17+t20+t22;
			double t25 = t21*u3;
			double t26 = t21*u4;
			double t27 = t21*u2;
			double t28 = t21*u5;
			double t29 = t10*t19*t21;
			double t30 = t10*t21;
			double t31 = t21+t30;

			f_udg[0*ng+i] = 0.0;
			f_udg[1*ng+i] = -t2*t3-t10*t15;
			f_udg[2*ng+i] = -t2*u2*u3;
			f_udg[3*ng+i] = -t2*u2*u4;
			f_udg[4*ng+i] = -t23*u2;
			f_udg[5*ng+i] = 0.0;
			f_udg[6*ng+i] = -t2*u2*u3;
			f_udg[7*ng+i] = -t2*t5-t10*t15;
			f_udg[8*ng+i] = -t2*u3*u4;
			f_udg[9*ng+i] = -t23*u3;
			f_udg[10*ng+i] = 0.0;
			f_udg[11*ng+i] = -t2*u2*u4;
			f_udg[12*ng+i] = -t2*u3*u4;
			f_udg[13*ng+i] = -t2*t6-t10*t15;
			f_udg[14*ng+i] = -t23*u4;
			f_udg[15*ng+i] = 1.0;
			f_udg[16*ng+i] = t21*u2*2.0-t10*t21*u2;
			f_udg[17*ng+i] = t25;
			f_udg[18*ng+i] = t26;
			f_udg[19*ng+i] = t28+t29-t2*t3*t10;
			f_udg[20*ng+i] = 0.0;
			f_udg[21*ng+i] = t25;
			f_udg[22*ng+i] = -t10*t21*u2;
			f_udg[23*ng+i] = 0.0;
			f_udg[24*ng+i] = -t2*t10*u2*u3;
			f_udg[25*ng+i] = 0.0;
			f_udg[26*ng+i] = t26;
			f_udg[27*ng+i] = 0.0;
			f_udg[28*ng+i] = -t10*t21*u2;
			f_udg[29*ng+i] = -t2*t10*u2*u4;
			f_udg[30*ng+i] = 0.0;
			f_udg[31*ng+i] = -t10*t21*u3;
			f_udg[32*ng+i] = t27;
			f_udg[33*ng+i] = 0.0;
			f_udg[34*ng+i] = -t2*t10*u2*u3;
			f_udg[35*ng+i] = 1.0;
			f_udg[36*ng+i] = t27;
			f_udg[37*ng+i] = t21*u3*2.0-t10*t21*u3;
			f_udg[38*ng+i] = t26;
			f_udg[39*ng+i] = t28+t29-t2*t5*t10;
			f_udg[40*ng+i] = 0.0;
			f_udg[41*ng+i] = 0.0;
			f_udg[42*ng+i] = t26;
			f_udg[43*ng+i] = -t10*t21*u3;
			f_udg[44*ng+i] = -t2*t10*u3*u4;
			f_udg[45*ng+i] = 0.0;
			f_udg[46*ng+i] = -t10*t21*u4;
			f_udg[47*ng+i] = 0.0;
			f_udg[48*ng+i] = t27;
			f_udg[49*ng+i] = -t2*t10*u2*u4;
			f_udg[50*ng+i] = 0.0;
			f_udg[51*ng+i] = 0.0;
			f_udg[52*ng+i] = -t10*t21*u4;
			f_udg[53*ng+i] = t25;
			f_udg[54*ng+i] = -t2*t10*u3*u4;
			f_udg[55*ng+i] = 1.0;
			f_udg[56*ng+i] = t27;
			f_udg[57*ng+i] = t25;
			f_udg[58*ng+i] = t21*u4*2.0-t10*t21*u4;
			f_udg[59*ng+i] = t28+t29-t2*t6*t10;
			f_udg[60*ng+i] = 0.0;
			f_udg[61*ng+i] = t10;
			f_udg[62*ng+i] = 0.0;
			f_udg[63*ng+i] = 0.0;
			f_udg[64*ng+i] = t31*u2;
			f_udg[65*ng+i] = 0.0;
			f_udg[66*ng+i] = 0.0;
			f_udg[67*ng+i] = t10;
			f_udg[68*ng+i] = 0.0;
			f_udg[69*ng+i] = t31*u3;
			f_udg[70*ng+i] = 0.0;
			f_udg[71*ng+i] = 0.0;
			f_udg[72*ng+i] = 0.0;
			f_udg[73*ng+i] = t10;
			f_udg[74*ng+i] = t31*u4;
		}
	}
}
