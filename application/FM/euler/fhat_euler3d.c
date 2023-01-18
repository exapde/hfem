void fhat_euler3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = 1.0/uh1;
		double t5 = uh3*uh3;
		double t6 = t2*t3*(1.0/2.0);
		double t7 = t3*t5*(1.0/2.0);
		double t8 = uh4*uh4;
		double t9 = t3*t8*(1.0/2.0);
		double t10 = t6+t7+t9;
		double t14 = t10*uh1;
		double t11 = -t14+uh5;
		double t12 = param1-1.0;
		double t13 = t11*t12;
		double t15 = param1*uh1*uh5*2.0;
		double t17 = param1*t2;
		double t18 = param1*t5;
		double t19 = param1*t8;
		double t16 = t2+t5+t8+t15-t17-t18-t19;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)+nl1*(t13+t2*t4)+nl2*t4*uh2*uh3+nl3*t4*uh2*uh4;
		fh[2*ng+i] = param6*(u3-uh3)+nl2*(t13+t4*t5)+nl1*t4*uh2*uh3+nl3*t4*uh3*uh4;
		fh[3*ng+i] = param6*(u4-uh4)+nl3*(t13+t4*t8)+nl1*t4*uh2*uh4+nl2*t4*uh3*uh4;
		fh[4*ng+i] = param6*(u5-uh5)+nl1*t3*t16*uh2*(1.0/2.0)+nl2*t3*t16*uh3*(1.0/2.0)+nl3*t3*t16*uh4*(1.0/2.0);

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		fh_udg[0*ng+i] = param6;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = 0.0;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = param6;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = 0.0;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = 0.0;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = param6;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = 0.0;
		fh_udg[16*ng+i] = 0.0;
		fh_udg[17*ng+i] = 0.0;
		fh_udg[18*ng+i] = param6;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = 0.0;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = param6;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = 1.0/(uh1*uh1);
		double t3 = uh2*uh2;
		double t4 = 1.0/(uh1*uh1*uh1);
		double t5 = uh3*uh3;
		double t6 = uh4*uh4;
		double t7 = param1-1.0;
		double t8 = t2*t3*(1.0/2.0);
		double t9 = t2*t5*(1.0/2.0);
		double t10 = t2*t6*(1.0/2.0);
		double t11 = t3*t4;
		double t12 = t4*t5;
		double t13 = t4*t6;
		double t14 = t11+t12+t13;
		double t17 = t14*uh1;
		double t15 = t8+t9+t10-t17;
		double t16 = t7*t15;
		double t18 = param1*uh1*uh5*2.0;
		double t20 = param1*t3;
		double t21 = param1*t5;
		double t22 = param1*t6;
		double t19 = t3+t5+t6+t18-t20-t21-t22;
		double t23 = 1.0/uh1;
		double t24 = uh2*2.0;
		double t26 = param1*uh2*2.0;
		double t25 = t24-t26;
		double t27 = nl3*t23*uh4;
		double t28 = uh3*2.0;
		double t30 = param1*uh3*2.0;
		double t29 = t28-t30;
		double t31 = nl1*t23*uh2;
		double t32 = nl2*t23*uh3;
		double t33 = uh4*2.0;
		double t35 = param1*uh4*2.0;
		double t34 = t33-t35;
		fh_uh[0*ng+i] = -param6;
		fh_uh[1*ng+i] = -nl1*(t16+t2*t3)-nl2*t2*uh2*uh3-nl3*t2*uh2*uh4;
		fh_uh[2*ng+i] = -nl2*(t16+t2*t5)-nl1*t2*uh2*uh3-nl3*t2*uh3*uh4;
		fh_uh[3*ng+i] = -nl3*(t16+t2*t6)-nl1*t2*uh2*uh4-nl2*t2*uh3*uh4;
		fh_uh[4*ng+i] = -nl1*t4*t19*uh2-nl2*t4*t19*uh3-nl3*t4*t19*uh4+nl1*param1*t2*uh2*uh5+nl2*param1*t2*uh3*uh5+nl3*param1*t2*uh4*uh5;
		fh_uh[5*ng+i] = nl1;
		fh_uh[6*ng+i] = -param6+t27+t32+nl1*(t23*uh2*2.0-t7*t23*uh2);
		fh_uh[7*ng+i] = nl1*t23*uh3-nl2*t7*t23*uh2;
		fh_uh[8*ng+i] = nl1*t23*uh4-nl3*t7*t23*uh2;
		fh_uh[9*ng+i] = nl1*t2*t19*(1.0/2.0)+nl1*t2*t25*uh2*(1.0/2.0)+nl2*t2*t25*uh3*(1.0/2.0)+nl3*t2*t25*uh4*(1.0/2.0);
		fh_uh[10*ng+i] = nl2;
		fh_uh[11*ng+i] = nl2*t23*uh2-nl1*t7*t23*uh3;
		fh_uh[12*ng+i] = -param6+t27+t31+nl2*(t23*uh3*2.0-t7*t23*uh3);
		fh_uh[13*ng+i] = nl2*t23*uh4-nl3*t7*t23*uh3;
		fh_uh[14*ng+i] = nl2*t2*t19*(1.0/2.0)+nl1*t2*t29*uh2*(1.0/2.0)+nl2*t2*t29*uh3*(1.0/2.0)+nl3*t2*t29*uh4*(1.0/2.0);
		fh_uh[15*ng+i] = nl3;
		fh_uh[16*ng+i] = nl3*t23*uh2-nl1*t7*t23*uh4;
		fh_uh[17*ng+i] = nl3*t23*uh3-nl2*t7*t23*uh4;
		fh_uh[18*ng+i] = -param6+t31+t32+nl3*(t23*uh4*2.0-t7*t23*uh4);
		fh_uh[19*ng+i] = nl3*t2*t19*(1.0/2.0)+nl1*t2*t34*uh2*(1.0/2.0)+nl2*t2*t34*uh3*(1.0/2.0)+nl3*t2*t34*uh4*(1.0/2.0);
		fh_uh[20*ng+i] = 0.0;
		fh_uh[21*ng+i] = nl1*t7;
		fh_uh[22*ng+i] = nl2*t7;
		fh_uh[23*ng+i] = nl3*t7;
		fh_uh[24*ng+i] = -param6+nl1*param1*t23*uh2+nl2*param1*t23*uh3+nl3*param1*t23*uh4;

	}
}

void fhatonly_euler3d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = 1.0/uh1;
		double t5 = uh3*uh3;
		double t6 = t2*t3*(1.0/2.0);
		double t7 = t3*t5*(1.0/2.0);
		double t8 = uh4*uh4;
		double t9 = t3*t8*(1.0/2.0);
		double t10 = t6+t7+t9;
		double t14 = t10*uh1;
		double t11 = -t14+uh5;
		double t12 = param1-1.0;
		double t13 = t11*t12;
		double t15 = param1*uh1*uh5*2.0;
		double t17 = param1*t2;
		double t18 = param1*t5;
		double t19 = param1*t8;
		double t16 = t2+t5+t8+t15-t17-t18-t19;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)+nl1*(t13+t2*t4)+nl2*t4*uh2*uh3+nl3*t4*uh2*uh4;
		fh[2*ng+i] = param6*(u3-uh3)+nl2*(t13+t4*t5)+nl1*t4*uh2*uh3+nl3*t4*uh3*uh4;
		fh[3*ng+i] = param6*(u4-uh4)+nl3*(t13+t4*t8)+nl1*t4*uh2*uh4+nl2*t4*uh3*uh4;
		fh[4*ng+i] = param6*(u5-uh5)+nl1*t3*t16*uh2*(1.0/2.0)+nl2*t3*t16*uh3*(1.0/2.0)+nl3*t3*t16*uh4*(1.0/2.0);

	}
}

