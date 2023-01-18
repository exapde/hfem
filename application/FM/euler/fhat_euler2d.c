void fhat_euler2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = 1.0/uh1;
		double t5 = uh3*uh3;
		double t6 = t2*t3*(1.0/2.0);
		double t7 = t3*t5*(1.0/2.0);
		double t8 = t6+t7;
		double t9 = uh4-t8*uh1;
		double t10 = param1-1.0;
		double t11 = t9*t10;
		double t12 = param1*uh1*uh4*2.0;
		double t13 = t2+t5+t12-param1*t2-param1*t5;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)+nl1*(t11+t2*t4)+nl2*t4*uh2*uh3;
		fh[2*ng+i] = param6*(u3-uh3)+nl2*(t11+t4*t5)+nl1*t4*uh2*uh3;
		fh[3*ng+i] = param6*(u4-uh4)+nl1*t3*t13*uh2*(1.0/2.0)+nl2*t3*t13*uh3*(1.0/2.0);

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		fh_udg[0*ng+i] = param6;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = 0.0;
		fh_udg[5*ng+i] = param6;
		fh_udg[6*ng+i] = 0.0;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = 0.0;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = param6;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = 0.0;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = param6;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/(uh1*uh1);
		double t3 = uh2*uh2;
		double t4 = 1.0/(uh1*uh1*uh1);
		double t5 = uh3*uh3;
		double t6 = param1-1.0;
		double t7 = t2*t3*(1.0/2.0);
		double t8 = t2*t5*(1.0/2.0);
		double t9 = t3*t4;
		double t10 = t4*t5;
		double t11 = t9+t10;
		double t12 = t7+t8-t11*uh1;
		double t13 = t6*t12;
		double t14 = param1*uh1*uh4*2.0;
		double t17 = param1*t3;
		double t18 = param1*t5;
		double t15 = t3+t5+t14-t17-t18;
		double t16 = 1.0/uh1;
		double t19 = uh2*2.0;
		double t20 = t19-param1*uh2*2.0;
		double t21 = uh3*2.0;
		double t22 = t21-param1*uh3*2.0;
		fh_uh[0*ng+i] = -param6;
		fh_uh[1*ng+i] = -nl1*(t13+t2*t3)-nl2*t2*uh2*uh3;
		fh_uh[2*ng+i] = -nl2*(t13+t2*t5)-nl1*t2*uh2*uh3;
		fh_uh[3*ng+i] = -nl1*t4*t15*uh2-nl2*t4*t15*uh3+nl1*param1*t2*uh2*uh4+nl2*param1*t2*uh3*uh4;
		fh_uh[4*ng+i] = nl1;
		fh_uh[5*ng+i] = -param6+nl1*(t16*uh2*2.0-t6*t16*uh2)+nl2*t16*uh3;
		fh_uh[6*ng+i] = nl1*t16*uh3-nl2*t6*t16*uh2;
		fh_uh[7*ng+i] = nl1*t2*t15*(1.0/2.0)+nl1*t2*t20*uh2*(1.0/2.0)+nl2*t2*t20*uh3*(1.0/2.0);
		fh_uh[8*ng+i] = nl2;
		fh_uh[9*ng+i] = nl2*t16*uh2-nl1*t6*t16*uh3;
		fh_uh[10*ng+i] = -param6+nl2*(t16*uh3*2.0-t6*t16*uh3)+nl1*t16*uh2;
		fh_uh[11*ng+i] = nl2*t2*t15*(1.0/2.0)+nl1*t2*t22*uh2*(1.0/2.0)+nl2*t2*t22*uh3*(1.0/2.0);
		fh_uh[12*ng+i] = 0.0;
		fh_uh[13*ng+i] = nl1*t6;
		fh_uh[14*ng+i] = nl2*t6;
		fh_uh[15*ng+i] = -param6+nl1*param1*t16*uh2+nl2*param1*t16*uh3;

	}
}

void fhatonly_euler2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = 1.0/uh1;
		double t5 = uh3*uh3;
		double t6 = t2*t3*(1.0/2.0);
		double t7 = t3*t5*(1.0/2.0);
		double t8 = t6+t7;
		double t9 = uh4-t8*uh1;
		double t10 = param1-1.0;
		double t11 = t9*t10;
		double t12 = param1*uh1*uh4*2.0;
		double t13 = t2+t5+t12-param1*t2-param1*t5;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)+nl1*(t11+t2*t4)+nl2*t4*uh2*uh3;
		fh[2*ng+i] = param6*(u3-uh3)+nl2*(t11+t4*t5)+nl1*t4*uh2*uh3;
		fh[3*ng+i] = param6*(u4-uh4)+nl1*t3*t13*uh2*(1.0/2.0)+nl2*t3*t13*uh3*(1.0/2.0);

	}
}

