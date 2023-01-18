void fhat_mhd23d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double uh8 = uh[7*ng+i];
		double uh9 = uh[8*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = uh6*uh6;
		double t5 = t4*(1.0/2.0);
		double t6 = uh7*uh7;
		double t7 = t6*(1.0/2.0);
		double t8 = uh8*uh8;
		double t9 = t8*(1.0/2.0);
		double t10 = 1.0/uh1;
		double t11 = uh3*uh3;
		double t12 = param1-1.0;
		double t13 = t2*t3*(1.0/2.0);
		double t14 = t3*t11*(1.0/2.0);
		double t15 = uh4*uh4;
		double t16 = t3*t15*(1.0/2.0);
		double t17 = t13+t14+t16;
		double t18 = t17*uh1;
		double t19 = t5+t7+t9+t18-uh5;
		double t20 = uh6*uh7;
		double t21 = t20-t10*uh2*uh3;
		double t22 = t10*uh2*uh7;
		double t23 = t22-t10*uh3*uh6;
		double t24 = param3*param3;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = -nl2*t21+nl1*(-t5+t7+t9+t2*t10-t12*t19)+param2*(u2-uh2);
		fh[2*ng+i] = -nl1*t21+nl2*(t5-t7+t9+t10*t11-t12*t19)+param2*(u3-uh3);
		fh[3*ng+i] = -nl1*(uh6*uh8-t10*uh2*uh4)-nl2*(uh7*uh8-t10*uh3*uh4)+param2*(u4-uh4);
		fh[4*ng+i] = param2*(u5-uh5)-nl1*t3*(-t2*uh2-t11*uh2-t15*uh2+param1*t2*uh2+param1*t11*uh2+param1*t15*uh2-t6*uh1*uh2*2.0-t8*uh1*uh2*2.0+param1*t4*uh1*uh2+param1*t6*uh1*uh2+param1*t8*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t3*(-t2*uh3-t11*uh3-t15*uh3+param1*t2*uh3+param1*t11*uh3+param1*t15*uh3-t4*uh1*uh3*2.0-t8*uh1*uh3*2.0+param1*t4*uh1*uh3+param1*t6*uh1*uh3+param1*t8*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)*(1.0/2.0);
		fh[5*ng+i] = -nl2*t23+nl1*uh9+param2*(u6-uh6);
		fh[6*ng+i] = nl1*t23+nl2*uh9+param2*(u7-uh7);
		fh[7*ng+i] = param2*(u8-uh8)+nl1*(t10*uh2*uh8-t10*uh4*uh6)+nl2*(t10*uh3*uh8-t10*uh4*uh7);
		fh[8*ng+i] = param2*(u9-uh9)+nl1*t24*uh6+nl2*t24*uh7;

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
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double uh8 = uh[7*ng+i];
		double uh9 = uh[8*ng+i];
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
		fh_udg[8*ng+i] = 0.0;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = param2;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = 0.0;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = 0.0;
		fh_udg[16*ng+i] = 0.0;
		fh_udg[17*ng+i] = 0.0;
		fh_udg[18*ng+i] = 0.0;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = param2;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = 0.0;
		fh_udg[25*ng+i] = 0.0;
		fh_udg[26*ng+i] = 0.0;
		fh_udg[27*ng+i] = 0.0;
		fh_udg[28*ng+i] = 0.0;
		fh_udg[29*ng+i] = 0.0;
		fh_udg[30*ng+i] = param2;
		fh_udg[31*ng+i] = 0.0;
		fh_udg[32*ng+i] = 0.0;
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
		fh_udg[48*ng+i] = 0.0;
		fh_udg[49*ng+i] = 0.0;
		fh_udg[50*ng+i] = param2;
		fh_udg[51*ng+i] = 0.0;
		fh_udg[52*ng+i] = 0.0;
		fh_udg[53*ng+i] = 0.0;
		fh_udg[54*ng+i] = 0.0;
		fh_udg[55*ng+i] = 0.0;
		fh_udg[56*ng+i] = 0.0;
		fh_udg[57*ng+i] = 0.0;
		fh_udg[58*ng+i] = 0.0;
		fh_udg[59*ng+i] = 0.0;
		fh_udg[60*ng+i] = param2;
		fh_udg[61*ng+i] = 0.0;
		fh_udg[62*ng+i] = 0.0;
		fh_udg[63*ng+i] = 0.0;
		fh_udg[64*ng+i] = 0.0;
		fh_udg[65*ng+i] = 0.0;
		fh_udg[66*ng+i] = 0.0;
		fh_udg[67*ng+i] = 0.0;
		fh_udg[68*ng+i] = 0.0;
		fh_udg[69*ng+i] = 0.0;
		fh_udg[70*ng+i] = param2;
		fh_udg[71*ng+i] = 0.0;
		fh_udg[72*ng+i] = 0.0;
		fh_udg[73*ng+i] = 0.0;
		fh_udg[74*ng+i] = 0.0;
		fh_udg[75*ng+i] = 0.0;
		fh_udg[76*ng+i] = 0.0;
		fh_udg[77*ng+i] = 0.0;
		fh_udg[78*ng+i] = 0.0;
		fh_udg[79*ng+i] = 0.0;
		fh_udg[80*ng+i] = param2;

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
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double uh8 = uh[7*ng+i];
		double uh9 = uh[8*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

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
		double t15 = t8+t9+t10-t14*uh1;
		double t16 = t7*t15;
		double t17 = uh7*uh7;
		double t18 = uh8*uh8;
		double t19 = uh6*uh6;
		double t20 = t2*uh2*uh7;
		double t21 = t20-t2*uh3*uh6;
		double t22 = 1.0/uh1;
		double t23 = param1*uh2*uh3*2.0;
		double t24 = uh1*uh6*uh7*2.0;
		double t25 = t23+t24-uh2*uh3*2.0;
		double t26 = t18*uh1*2.0;
		double t27 = param1*uh1*uh5*2.0;
		double t28 = nl1*t22*uh2;
		double t29 = nl2*t22*uh3;
		double t30 = t7*uh6;
		double t31 = nl1*t22*uh3;
		double t32 = nl1*t22*uh4;
		double t33 = uh1*uh4*uh8*2.0;
		double t34 = nl2*t22*uh2;
		double t35 = nl2*t22*uh4;
		double t36 = param3*param3;
		double t37 = uh8-t7*uh8;
		double t38 = -param2+t28+t29;
		fh_uh[0*ng+i] = -param2;
		fh_uh[1*ng+i] = -nl1*(t16+t2*t3)-nl2*t2*uh2*uh3;
		fh_uh[2*ng+i] = -nl2*(t16+t2*t5)-nl1*t2*uh2*uh3;
		fh_uh[3*ng+i] = -nl1*t2*uh2*uh4-nl2*t2*uh3*uh4;
		fh_uh[4*ng+i] = nl1*t4*(-t3*uh2-t5*uh2-t6*uh2+param1*t3*uh2+param1*t5*uh2+param1*t6*uh2-t17*uh1*uh2*2.0-t18*uh1*uh2*2.0+param1*t17*uh1*uh2+param1*t18*uh1*uh2+param1*t19*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)+nl2*t4*(-t3*uh3-t5*uh3-t6*uh3+param1*t3*uh3+param1*t5*uh3+param1*t6*uh3-t18*uh1*uh3*2.0-t19*uh1*uh3*2.0+param1*t17*uh1*uh3+param1*t18*uh1*uh3+param1*t19*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)-nl1*t2*(t17*uh2*-2.0-t18*uh2*2.0+param1*t17*uh2+param1*t18*uh2+param1*t19*uh2-param1*uh2*uh5*2.0+uh3*uh6*uh7*2.0+uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t2*(t18*uh3*-2.0-t19*uh3*2.0+param1*t17*uh3+param1*t18*uh3+param1*t19*uh3-param1*uh3*uh5*2.0+uh2*uh6*uh7*2.0+uh4*uh7*uh8*2.0)*(1.0/2.0);
		fh_uh[5*ng+i] = nl2*t21;
		fh_uh[6*ng+i] = -nl1*t21;
		fh_uh[7*ng+i] = -nl1*(t2*uh2*uh8-t2*uh4*uh6)-nl2*(t2*uh3*uh8-t2*uh4*uh7);
		fh_uh[8*ng+i] = 0.0;
		fh_uh[9*ng+i] = nl1;
		fh_uh[10*ng+i] = -param2+t29+nl1*(t22*uh2*2.0-t7*t22*uh2);
		fh_uh[11*ng+i] = t31-nl2*t7*t22*uh2;
		fh_uh[12*ng+i] = t32;
		fh_uh[13*ng+i] = nl1*t2*(t3*3.0+t5+t6+t26+t27-param1*t3*3.0-param1*t5-param1*t6+t17*uh1*2.0-param1*t17*uh1-param1*t18*uh1-param1*t19*uh1)*(1.0/2.0)-nl2*t2*t25*(1.0/2.0);
		fh_uh[14*ng+i] = -nl2*t22*uh7;
		fh_uh[15*ng+i] = nl1*t22*uh7;
		fh_uh[16*ng+i] = nl1*t22*uh8;
		fh_uh[17*ng+i] = 0.0;
		fh_uh[18*ng+i] = nl2;
		fh_uh[19*ng+i] = t34-nl1*t7*t22*uh3;
		fh_uh[20*ng+i] = -param2+t28+nl2*(t22*uh3*2.0-t7*t22*uh3);
		fh_uh[21*ng+i] = t35;
		fh_uh[22*ng+i] = nl2*t2*(t3+t5*3.0+t6+t26+t27-param1*t3-param1*t5*3.0-param1*t6+t19*uh1*2.0-param1*t17*uh1-param1*t18*uh1-param1*t19*uh1)*(1.0/2.0)-nl1*t2*t25*(1.0/2.0);
		fh_uh[23*ng+i] = nl2*t22*uh6;
		fh_uh[24*ng+i] = -nl1*t22*uh6;
		fh_uh[25*ng+i] = nl2*t22*uh8;
		fh_uh[26*ng+i] = 0.0;
		fh_uh[27*ng+i] = 0.0;
		fh_uh[28*ng+i] = -nl1*t7*t22*uh4;
		fh_uh[29*ng+i] = -nl2*t7*t22*uh4;
		fh_uh[30*ng+i] = t38;
		fh_uh[31*ng+i] = nl1*t2*(uh2*uh4*-2.0+param1*uh2*uh4*2.0+uh1*uh6*uh8*2.0)*(-1.0/2.0)-nl2*t2*(uh3*uh4*-2.0+param1*uh3*uh4*2.0+uh1*uh7*uh8*2.0)*(1.0/2.0);
		fh_uh[32*ng+i] = 0.0;
		fh_uh[33*ng+i] = 0.0;
		fh_uh[34*ng+i] = -nl1*t22*uh6-nl2*t22*uh7;
		fh_uh[35*ng+i] = 0.0;
		fh_uh[36*ng+i] = 0.0;
		fh_uh[37*ng+i] = nl1*t7;
		fh_uh[38*ng+i] = nl2*t7;
		fh_uh[39*ng+i] = 0.0;
		fh_uh[40*ng+i] = -param2+nl1*param1*t22*uh2+nl2*param1*t22*uh3;
		fh_uh[41*ng+i] = 0.0;
		fh_uh[42*ng+i] = 0.0;
		fh_uh[43*ng+i] = 0.0;
		fh_uh[44*ng+i] = 0.0;
		fh_uh[45*ng+i] = 0.0;
		fh_uh[46*ng+i] = -nl1*(t30+uh6)-nl2*uh7;
		fh_uh[47*ng+i] = -nl1*uh7-nl2*(t30-uh6);
		fh_uh[48*ng+i] = -nl1*uh8;
		fh_uh[49*ng+i] = nl1*t2*(t33+uh1*uh3*uh7*2.0+param1*uh1*uh2*uh6*2.0)*(-1.0/2.0)-nl2*t2*(uh1*uh2*uh7*2.0-uh1*uh3*uh6*4.0+param1*uh1*uh3*uh6*2.0)*(1.0/2.0);
		fh_uh[50*ng+i] = -param2+t29;
		fh_uh[51*ng+i] = -t31;
		fh_uh[52*ng+i] = -t32;
		fh_uh[53*ng+i] = nl1*t36;
		fh_uh[54*ng+i] = 0.0;
		fh_uh[55*ng+i] = -nl2*uh6+nl1*(uh7-t7*uh7);
		fh_uh[56*ng+i] = -nl1*uh6-nl2*(uh7+t7*uh7);
		fh_uh[57*ng+i] = -nl2*uh8;
		fh_uh[58*ng+i] = nl2*t2*(t33+uh1*uh2*uh6*2.0+param1*uh1*uh3*uh7*2.0)*(-1.0/2.0)-nl1*t2*(uh1*uh2*uh7*-4.0+uh1*uh3*uh6*2.0+param1*uh1*uh2*uh7*2.0)*(1.0/2.0);
		fh_uh[59*ng+i] = -t34;
		fh_uh[60*ng+i] = -param2+t28;
		fh_uh[61*ng+i] = -t35;
		fh_uh[62*ng+i] = nl2*t36;
		fh_uh[63*ng+i] = 0.0;
		fh_uh[64*ng+i] = nl1*t37;
		fh_uh[65*ng+i] = nl2*t37;
		fh_uh[66*ng+i] = -nl1*uh6-nl2*uh7;
		fh_uh[67*ng+i] = nl1*t2*(uh1*uh2*uh8*-4.0+uh1*uh4*uh6*2.0+param1*uh1*uh2*uh8*2.0)*(-1.0/2.0)-nl2*t2*(uh1*uh3*uh8*-4.0+uh1*uh4*uh7*2.0+param1*uh1*uh3*uh8*2.0)*(1.0/2.0);
		fh_uh[68*ng+i] = 0.0;
		fh_uh[69*ng+i] = 0.0;
		fh_uh[70*ng+i] = t38;
		fh_uh[71*ng+i] = 0.0;
		fh_uh[72*ng+i] = 0.0;
		fh_uh[73*ng+i] = 0.0;
		fh_uh[74*ng+i] = 0.0;
		fh_uh[75*ng+i] = 0.0;
		fh_uh[76*ng+i] = 0.0;
		fh_uh[77*ng+i] = nl1;
		fh_uh[78*ng+i] = nl2;
		fh_uh[79*ng+i] = 0.0;
		fh_uh[80*ng+i] = -param2;

	}
}

void fhatonly_mhd2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double uh6 = uh[5*ng+i];
		double uh7 = uh[6*ng+i];
		double uh8 = uh[7*ng+i];
		double uh9 = uh[8*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = uh2*uh2;
		double t3 = 1.0/(uh1*uh1);
		double t4 = uh6*uh6;
		double t5 = t4*(1.0/2.0);
		double t6 = uh7*uh7;
		double t7 = t6*(1.0/2.0);
		double t8 = uh8*uh8;
		double t9 = t8*(1.0/2.0);
		double t10 = 1.0/uh1;
		double t11 = uh3*uh3;
		double t12 = param1-1.0;
		double t13 = t2*t3*(1.0/2.0);
		double t14 = t3*t11*(1.0/2.0);
		double t15 = uh4*uh4;
		double t16 = t3*t15*(1.0/2.0);
		double t17 = t13+t14+t16;
		double t18 = t17*uh1;
		double t19 = t5+t7+t9+t18-uh5;
		double t20 = uh6*uh7;
		double t21 = t20-t10*uh2*uh3;
		double t22 = t10*uh2*uh7;
		double t23 = t22-t10*uh3*uh6;
		double t24 = param3*param3;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = -nl2*t21+nl1*(-t5+t7+t9+t2*t10-t12*t19)+param2*(u2-uh2);
		fh[2*ng+i] = -nl1*t21+nl2*(t5-t7+t9+t10*t11-t12*t19)+param2*(u3-uh3);
		fh[3*ng+i] = -nl1*(uh6*uh8-t10*uh2*uh4)-nl2*(uh7*uh8-t10*uh3*uh4)+param2*(u4-uh4);
		fh[4*ng+i] = param2*(u5-uh5)-nl1*t3*(-t2*uh2-t11*uh2-t15*uh2+param1*t2*uh2+param1*t11*uh2+param1*t15*uh2-t6*uh1*uh2*2.0-t8*uh1*uh2*2.0+param1*t4*uh1*uh2+param1*t6*uh1*uh2+param1*t8*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t3*(-t2*uh3-t11*uh3-t15*uh3+param1*t2*uh3+param1*t11*uh3+param1*t15*uh3-t4*uh1*uh3*2.0-t8*uh1*uh3*2.0+param1*t4*uh1*uh3+param1*t6*uh1*uh3+param1*t8*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)*(1.0/2.0);
		fh[5*ng+i] = -nl2*t23+nl1*uh9+param2*(u6-uh6);
		fh[6*ng+i] = nl1*t23+nl2*uh9+param2*(u7-uh7);
		fh[7*ng+i] = param2*(u8-uh8)+nl1*(t10*uh2*uh8-t10*uh4*uh6)+nl2*(t10*uh3*uh8-t10*uh4*uh7);
		fh[8*ng+i] = param2*(u9-uh9)+nl1*t24*uh6+nl2*t24*uh7;

	}
}
