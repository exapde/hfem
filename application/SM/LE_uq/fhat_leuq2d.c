void fhat_leuq2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = u3+u6+2.0;
		double t3 = param2*t2;
		double t4 = u4+u5;
		fh[0*ng+i] = nl1*(t3+param1*(u3+1.0)*2.0)+param3*(u1-uh1)+nl2*param1*t4;
		fh[1*ng+i] = nl2*(t3+param1*(u6+1.0)*2.0)+param3*(u2-uh2)+nl1*param1*t4;

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = nl2*param1;
		double t3 = nl1*param1;
		double t4 = param1*2.0;
		double t5 = param2+t4;
		fh_udg[0*ng+i] = param3;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = param3;
		fh_udg[4*ng+i] = nl1*t5;
		fh_udg[5*ng+i] = nl2*param2;
		fh_udg[6*ng+i] = t2;
		fh_udg[7*ng+i] = t3;
		fh_udg[8*ng+i] = t2;
		fh_udg[9*ng+i] = t3;
		fh_udg[10*ng+i] = nl1*param2;
		fh_udg[11*ng+i] = nl2*t5;

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		fh_uh[0*ng+i] = -param3;
		fh_uh[1*ng+i] = 0.0;
		fh_uh[2*ng+i] = 0.0;
		fh_uh[3*ng+i] = -param3;

	}
}

void fhatonly_leuq2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = u3+u6+2.0;
		double t3 = param2*t2;
		double t4 = u4+u5;
		fh[0*ng+i] = nl1*(t3+param1*(u3+1.0)*2.0)+param3*(u1-uh1)+nl2*param1*t4;
		fh[1*ng+i] = nl2*(t3+param1*(u6+1.0)*2.0)+param3*(u2-uh2)+nl1*param1*t4;

	}
}

