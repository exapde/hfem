void fhat_mhd2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/uh1;
		double t3 = uh2*uh2;
		double t4 = uh3*uh3;
		double t5 = uh6*uh6;
		double t6 = uh2*uh3;
		double t7 = t6-uh1*uh5*uh6;
		double t8 = uh1*uh4*2.0;
		double t9 = param1*t3;
		double t10 = param1*t4;
		double t11 = uh5*uh5;
		double t12 = param1*t11*uh1;
		double t13 = param1*t5*uh1;
		double t14 = 1.0/(uh1*uh1);
		double t15 = t2*uh2*uh6;
		double t16 = t15-t2*uh3*uh5;
		double t17 = param2*param2;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = param2*(u2-uh2)-nl1*t2*(t3*-3.0-t4+t8+t9+t10+t12+t13-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl2*t2*t7;
		fh[2*ng+i] = param2*(u3-uh3)-nl2*t2*(-t3-t4*3.0+t8+t9+t10+t12+t13-t11*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl1*t2*t7;
		fh[3*ng+i] = param2*(u4-uh4)-nl1*t14*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t11*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t14*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t11*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t11*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh[4*ng+i] = -nl2*t16+nl1*uh7+param2*(u5-uh5);
		fh[5*ng+i] = nl1*t16+nl2*uh7+param2*(u6-uh6);
		fh[6*ng+i] = param2*(u7-uh7)+nl1*t17*uh5+nl2*t17*uh6;

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
		double u7 = udg[6*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		fh_udg[0*ng+i] = param2;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = 0.0;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = 0.0;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = param2;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = 0.0;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = 0.0;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = 0.0;
		fh_udg[16*ng+i] = param2;
		fh_udg[17*ng+i] = 0.0;
		fh_udg[18*ng+i] = 0.0;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = 0.0;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = param2;
		fh_udg[25*ng+i] = 0.0;
		fh_udg[26*ng+i] = 0.0;
		fh_udg[27*ng+i] = 0.0;
		fh_udg[28*ng+i] = 0.0;
		fh_udg[29*ng+i] = 0.0;
		fh_udg[30*ng+i] = 0.0;
		fh_udg[31*ng+i] = 0.0;
		fh_udg[32*ng+i] = param2;
		fh_udg[33*ng+i] = 0.0;
		fh_udg[34*ng+i] = 0.0;
		fh_udg[35*ng+i] = 0.0;
		fh_udg[36*ng+i] = 0.0;
		fh_udg[37*ng+i] = 0.0;
		fh_udg[38*ng+i] = 0.0;
		fh_udg[39*ng+i] = 0.0;
		fh_udg[40*ng+i] = param2;
		fh_udg[41*ng+i] = 0.0;
		fh_udg[42*ng+i] = 0.0;
		fh_udg[43*ng+i] = 0.0;
		fh_udg[44*ng+i] = 0.0;
		fh_udg[45*ng+i] = 0.0;
		fh_udg[46*ng+i] = 0.0;
		fh_udg[47*ng+i] = 0.0;
		fh_udg[48*ng+i] = param2;

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
		double u7 = udg[6*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/(uh1*uh1);
		double t3 = uh2*uh2;
		double t4 = uh3*uh3;
		double t5 = uh6*uh6;
		double t6 = uh5*uh5;
		double t7 = 1.0/uh1;
		double t8 = uh2*uh3;
		double t9 = t8-uh1*uh5*uh6;
		double t10 = uh1*uh4*2.0;
		double t11 = param1*t3;
		double t12 = param1*t4;
		double t13 = param1*t6*uh1;
		double t14 = param1*t5*uh1;
		double t15 = uh4*2.0;
		double t16 = param1*t6;
		double t17 = param1*t5;
		double t18 = 1.0/(uh1*uh1*uh1);
		double t19 = t2*uh2*uh6;
		double t20 = t19-t2*uh3*uh5;
		double t21 = param1*uh1*uh4*2.0;
		double t22 = param1*uh2*uh3*2.0;
		double t23 = uh1*uh5*uh6*2.0;
		double t24 = t22+t23-uh2*uh3*2.0;
		double t25 = uh1*2.0;
		double t26 = t25-param1*uh1*2.0;
		double t27 = nl2*t7*uh3;
		double t28 = nl1*t7*uh3;
		double t29 = nl2*t7*uh2;
		double t30 = nl1*t7*uh2;
		double t31 = param2*param2;
		fh_uh[0*ng+i] = -param2;
		fh_uh[1*ng+i] = nl1*t7*(t5*-2.0+t15+t16+t17-param1*uh4*2.0)*(-1.0/2.0)+nl1*t2*(t3*-3.0-t4+t10+t11+t12+t13+t14-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)-nl2*t2*t9-nl2*t7*uh5*uh6;
		fh_uh[2*ng+i] = nl2*t7*(t6*-2.0+t15+t16+t17-param1*uh4*2.0)*(-1.0/2.0)+nl2*t2*(-t3-t4*3.0+t10+t11+t12+t13+t14-t6*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)-nl1*t2*t9-nl1*t7*uh5*uh6;
		fh_uh[3*ng+i] = nl1*t18*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t6*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)+nl2*t18*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t6*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t6*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)-nl1*t2*(t5*uh2*-2.0+param1*t5*uh2+param1*t6*uh2-param1*uh2*uh4*2.0+uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t2*(t6*uh3*-2.0+param1*t5*uh3+param1*t6*uh3-param1*uh3*uh4*2.0+uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh_uh[4*ng+i] = nl2*t20;
		fh_uh[5*ng+i] = -nl1*t20;
		fh_uh[6*ng+i] = 0.0;
		fh_uh[7*ng+i] = nl1;
		fh_uh[8*ng+i] = -param2+t27+nl1*t7*(uh2*6.0-param1*uh2*2.0)*(1.0/2.0);
		fh_uh[9*ng+i] = t28+nl2*t7*(uh2*2.0-param1*uh2*2.0)*(1.0/2.0);
		fh_uh[10*ng+i] = nl1*t2*(t3*3.0+t4-t12-t13-t14+t21-param1*t3*3.0+t5*uh1*2.0)*(1.0/2.0)-nl2*t2*t24*(1.0/2.0);
		fh_uh[11*ng+i] = -nl2*t7*uh6;
		fh_uh[12*ng+i] = nl1*t7*uh6;
		fh_uh[13*ng+i] = 0.0;
		fh_uh[14*ng+i] = nl2;
		fh_uh[15*ng+i] = t29+nl1*t7*(uh3*2.0-param1*uh3*2.0)*(1.0/2.0);
		fh_uh[16*ng+i] = -param2+t30+nl2*t7*(uh3*6.0-param1*uh3*2.0)*(1.0/2.0);
		fh_uh[17*ng+i] = nl2*t2*(t3+t4*3.0-t11-t13-t14+t21-param1*t4*3.0+t6*uh1*2.0)*(1.0/2.0)-nl1*t2*t24*(1.0/2.0);
		fh_uh[18*ng+i] = nl2*t7*uh5;
		fh_uh[19*ng+i] = -nl1*t7*uh5;
		fh_uh[20*ng+i] = 0.0;
		fh_uh[21*ng+i] = 0.0;
		fh_uh[22*ng+i] = nl1*t7*t26*(-1.0/2.0);
		fh_uh[23*ng+i] = nl2*t7*t26*(-1.0/2.0);
		fh_uh[24*ng+i] = -param2+nl1*param1*t7*uh2+nl2*param1*t7*uh3;
		fh_uh[25*ng+i] = 0.0;
		fh_uh[26*ng+i] = 0.0;
		fh_uh[27*ng+i] = 0.0;
		fh_uh[28*ng+i] = 0.0;
		fh_uh[29*ng+i] = -nl2*uh6-nl1*param1*uh5;
		fh_uh[30*ng+i] = -nl1*uh6+nl2*t7*(uh1*uh5*4.0-param1*uh1*uh5*2.0)*(1.0/2.0);
		fh_uh[31*ng+i] = nl2*t2*(uh1*uh2*uh6*2.0-uh1*uh3*uh5*4.0+param1*uh1*uh3*uh5*2.0)*(-1.0/2.0)-nl1*t2*(uh1*uh3*uh6*2.0+param1*uh1*uh2*uh5*2.0)*(1.0/2.0);
		fh_uh[32*ng+i] = -param2+t27;
		fh_uh[33*ng+i] = -t28;
		fh_uh[34*ng+i] = nl1*t31;
		fh_uh[35*ng+i] = 0.0;
		fh_uh[36*ng+i] = -nl2*uh5+nl1*t7*(uh1*uh6*4.0-param1*uh1*uh6*2.0)*(1.0/2.0);
		fh_uh[37*ng+i] = -nl1*uh5-nl2*param1*uh6;
		fh_uh[38*ng+i] = nl1*t2*(uh1*uh2*uh6*-4.0+uh1*uh3*uh5*2.0+param1*uh1*uh2*uh6*2.0)*(-1.0/2.0)-nl2*t2*(uh1*uh2*uh5*2.0+param1*uh1*uh3*uh6*2.0)*(1.0/2.0);
		fh_uh[39*ng+i] = -t29;
		fh_uh[40*ng+i] = -param2+t30;
		fh_uh[41*ng+i] = nl2*t31;
		fh_uh[42*ng+i] = 0.0;
		fh_uh[43*ng+i] = 0.0;
		fh_uh[44*ng+i] = 0.0;
		fh_uh[45*ng+i] = 0.0;
		fh_uh[46*ng+i] = nl1;
		fh_uh[47*ng+i] = nl2;
		fh_uh[48*ng+i] = -param2;

	}
}

void fhatonly_mhd2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = 1.0/uh1;
		double t3 = uh2*uh2;
		double t4 = uh3*uh3;
		double t5 = uh6*uh6;
		double t6 = uh2*uh3;
		double t7 = t6-uh1*uh5*uh6;
		double t8 = uh1*uh4*2.0;
		double t9 = param1*t3;
		double t10 = param1*t4;
		double t11 = uh5*uh5;
		double t12 = param1*t11*uh1;
		double t13 = param1*t5*uh1;
		double t14 = 1.0/(uh1*uh1);
		double t15 = t2*uh2*uh6;
		double t16 = t15-t2*uh3*uh5;
		double t17 = param2*param2;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = param2*(u2-uh2)-nl1*t2*(t3*-3.0-t4+t8+t9+t10+t12+t13-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl2*t2*t7;
		fh[2*ng+i] = param2*(u3-uh3)-nl2*t2*(-t3-t4*3.0+t8+t9+t10+t12+t13-t11*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl1*t2*t7;
		fh[3*ng+i] = param2*(u4-uh4)-nl1*t14*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t11*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t14*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t11*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t11*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh[4*ng+i] = -nl2*t16+nl1*uh7+param2*(u5-uh5);
		fh[5*ng+i] = nl1*t16+nl2*uh7+param2*(u6-uh6);
		fh[6*ng+i] = param2*(u7-uh7)+nl1*t17*uh5+nl2*t17*uh6;

	}
}

