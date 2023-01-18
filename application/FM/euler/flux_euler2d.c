void flux_euler2d(double *f, double *f_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];

		double t2 = u2*u2;
		double t3 = 1.0/(u1*u1);
		double t4 = 1.0/u1;
		double t5 = t2*t3*(1.0/2.0);
		double t6 = u3*u3;
		double t7 = t3*t6*(1.0/2.0);
		double t8 = t5+t7;
		double t12 = t8*u1;
		double t9 = -t12+u4;
		double t10 = param1-1.0;
		double t11 = t4*u2*u3;
		double t13 = t9*t10;
		double t14 = t4*u4;
		double t15 = t4*t9*t10;
		double t16 = t14+t15;
		f[0*ng+i] = u2;
		f[1*ng+i] = t13+t2*t4;
		f[2*ng+i] = t11;
		f[3*ng+i] = t16*u2;
		f[4*ng+i] = u3;
		f[5*ng+i] = t11;
		f[6*ng+i] = t13+t4*t6;
		f[7*ng+i] = t16*u3;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];

		double t2 = 1.0/(u1*u1);
		double t3 = u2*u2;
		double t4 = 1.0/(u1*u1*u1);
		double t5 = u3*u3;
		double t6 = t2*t3*(1.0/2.0);
		double t7 = t2*t5*(1.0/2.0);
		double t8 = param1-1.0;
		double t9 = t3*t4;
		double t10 = t4*t5;
		double t11 = t9+t10;
		double t13 = t11*u1;
		double t12 = t6+t7-t13;
		double t14 = t2*u4;
		double t15 = t6+t7;
		double t21 = t15*u1;
		double t16 = -t21+u4;
		double t17 = t2*t8*t16;
		double t18 = 1.0/u1;
		double t19 = t8*t12*t18;
		double t20 = t14+t17+t19;
		double t22 = t18*u3;
		double t23 = t18*u2;
		double t24 = t18*u4;
		double t25 = t8*t16*t18;
		double t26 = t8*t18;
		double t27 = t18+t26;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = -t2*t3-t8*t12;
		f_udg[2*ng+i] = -t2*u2*u3;
		f_udg[3*ng+i] = -t20*u2;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = -t2*u2*u3;
		f_udg[6*ng+i] = -t2*t5-t8*t12;
		f_udg[7*ng+i] = -t20*u3;
		f_udg[8*ng+i] = 1.0;
		f_udg[9*ng+i] = t18*u2*2.0-t8*t18*u2;
		f_udg[10*ng+i] = t22;
		f_udg[11*ng+i] = t24+t25-t2*t3*t8;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = t22;
		f_udg[14*ng+i] = -t8*t18*u2;
		f_udg[15*ng+i] = -t2*t8*u2*u3;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = -t8*t18*u3;
		f_udg[18*ng+i] = t23;
		f_udg[19*ng+i] = -t2*t8*u2*u3;
		f_udg[20*ng+i] = 1.0;
		f_udg[21*ng+i] = t23;
		f_udg[22*ng+i] = t18*u3*2.0-t8*t18*u3;
		f_udg[23*ng+i] = t24+t25-t2*t5*t8;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t8;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = t27*u2;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = t8;
		f_udg[31*ng+i] = t27*u3;

	}
}

void fluxonly_euler2d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];

		double t2 = u2*u2;
		double t3 = 1.0/(u1*u1);
		double t4 = 1.0/u1;
		double t5 = t2*t3*(1.0/2.0);
		double t6 = u3*u3;
		double t7 = t3*t6*(1.0/2.0);
		double t8 = t5+t7;
		double t12 = t8*u1;
		double t9 = -t12+u4;
		double t10 = param1-1.0;
		double t11 = t4*u2*u3;
		double t13 = t9*t10;
		double t14 = t4*u4;
		double t15 = t4*t9*t10;
		double t16 = t14+t15;
		f[0*ng+i] = u2;
		f[1*ng+i] = t13+t2*t4;
		f[2*ng+i] = t11;
		f[3*ng+i] = t16*u2;
		f[4*ng+i] = u3;
		f[5*ng+i] = t11;
		f[6*ng+i] = t13+t4*t6;
		f[7*ng+i] = t16*u3;

	}
}

