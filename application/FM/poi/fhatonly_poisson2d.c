void fhatonly_poisson2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double uh1 = uh[0*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		fh[i] = param2*(u1-uh1)+nl1*param1*u2+nl2*param1*u3;
	}
}

