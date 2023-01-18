void flux_ledisp2d(double *f, double *f_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];

		double t2 = u4+u5;
		double t3 = param1*t2;
		double t4 = u3+u6;
		double t5 = param2*t4;
		f[0*ng+i] = t5+param1*u3*2.0;
		f[1*ng+i] = t3;
		f[2*ng+i] = t3;
		f[3*ng+i] = t5+param1*u6*2.0;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];

		double t2 = param1*2.0;
		double t3 = param2+t2;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = t3;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = param2;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = param1;
		f_udg[14*ng+i] = param1;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = param1;
		f_udg[18*ng+i] = param1;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = param2;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = t3;

	}
}

void fluxonly_ledisp2d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];

		double t2 = u4+u5;
		double t3 = param1*t2;
		double t4 = u3+u6;
		double t5 = param2*t4;
		f[0*ng+i] = t5+param1*u3*2.0;
		f[1*ng+i] = t3;
		f[2*ng+i] = t3;
		f[3*ng+i] = t5+param1*u6*2.0;

	}
}

