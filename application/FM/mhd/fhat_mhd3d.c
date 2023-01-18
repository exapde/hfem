void fhat_mhd3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double nl3 = nl[2*ng+i];

		double t2 = uh2*uh2;
		double t3 = uh3*uh3;
		double t4 = uh4*uh4;
		double t5 = uh7*uh7;
		double t6 = uh8*uh8;
		double t7 = 1.0/uh1;
		double t8 = uh6*uh6;
		double t9 = 1.0/(uh1*uh1);
		double t10 = t7*uh2*uh7;
		double t11 = t10-t7*uh3*uh6;
		double t12 = t7*uh2*uh8;
		double t13 = t12-t7*uh4*uh6;
		double t14 = t7*uh3*uh8;
		double t15 = t14-t7*uh4*uh7;
		double t16 = param2*param2;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param2*(u1-uh1);
		fh[1*ng+i] = t7*(nl1*t2*-3.0-nl1*t3-nl1*t4+nl1*param1*t2+nl1*param1*t3+nl1*param1*t4-nl1*t5*uh1*2.0-nl1*t6*uh1*2.0+nl1*uh1*uh5*2.0-nl2*uh2*uh3*2.0-nl3*uh2*uh4*2.0-param2*u2*uh1*2.0+param2*uh1*uh2*2.0+nl1*param1*t5*uh1+nl1*param1*t6*uh1+nl1*param1*t8*uh1-nl1*param1*uh1*uh5*2.0+nl2*uh1*uh6*uh7*2.0+nl3*uh1*uh6*uh8*2.0)*(-1.0/2.0);
		fh[2*ng+i] = t7*(-nl2*t2-nl2*t3*3.0-nl2*t4+nl2*param1*t2+nl2*param1*t3+nl2*param1*t4-nl2*t6*uh1*2.0-nl2*t8*uh1*2.0-nl1*uh2*uh3*2.0+nl2*uh1*uh5*2.0-nl3*uh3*uh4*2.0-param2*u3*uh1*2.0+param2*uh1*uh3*2.0+nl2*param1*t5*uh1+nl2*param1*t6*uh1+nl2*param1*t8*uh1-nl2*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh7*2.0+nl3*uh1*uh7*uh8*2.0)*(-1.0/2.0);
		fh[3*ng+i] = t7*(-nl3*t2-nl3*t3-nl3*t4*3.0+nl3*param1*t2+nl3*param1*t3+nl3*param1*t4-nl3*t5*uh1*2.0-nl3*t8*uh1*2.0-nl1*uh2*uh4*2.0-nl2*uh3*uh4*2.0+nl3*uh1*uh5*2.0-param2*u4*uh1*2.0+param2*uh1*uh4*2.0+nl3*param1*t5*uh1+nl3*param1*t6*uh1+nl3*param1*t8*uh1-nl3*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh8*2.0+nl2*uh1*uh7*uh8*2.0)*(-1.0/2.0);
		fh[4*ng+i] = param2*(u5-uh5)-nl1*t9*(-t2*uh2-t3*uh2-t4*uh2+param1*t2*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0-t6*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t6*uh1*uh2+param1*t8*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t9*(-t2*uh3-t3*uh3-t4*uh3+param1*t2*uh3+param1*t3*uh3+param1*t4*uh3-t6*uh1*uh3*2.0-t8*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t6*uh1*uh3+param1*t8*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)*(1.0/2.0)-nl3*t9*(-t2*uh4-t3*uh4-t4*uh4+param1*t2*uh4+param1*t3*uh4+param1*t4*uh4-t5*uh1*uh4*2.0-t8*uh1*uh4*2.0+param1*t5*uh1*uh4+param1*t6*uh1*uh4+param1*t8*uh1*uh4-param1*uh1*uh4*uh5*2.0+uh1*uh2*uh6*uh8*2.0+uh1*uh3*uh7*uh8*2.0)*(1.0/2.0);
		fh[5*ng+i] = nl2*t11+nl3*t13+nl1*uh9+param2*(u6-uh6);
		fh[6*ng+i] = -nl1*t11+nl3*t15+nl2*uh9+param2*(u7-uh7);
		fh[7*ng+i] = -nl1*t13-nl2*t15+nl3*uh9+param2*(u8-uh8);
		fh[8*ng+i] = param2*(u9-uh9)+nl1*t16*uh6+nl2*t16*uh7+nl3*t16*uh8;

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
		double nl3 = nl[2*ng+i];

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
		double x3 = pg[2*ng+i];
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
		double nl3 = nl[2*ng+i];

		double t2 = uh7*uh7;
		double t3 = uh8*uh8;
		double t4 = uh2*uh2;
		double t5 = uh3*uh3;
		double t6 = uh4*uh4;
		double t7 = uh6*uh6;
		double t8 = 1.0/uh1;
		double t9 = 1.0/(uh1*uh1);
		double t10 = 1.0/(uh1*uh1*uh1);
		double t11 = t9*uh2*uh7;
		double t12 = t11-t9*uh3*uh6;
		double t13 = t9*uh2*uh8;
		double t14 = t13-t9*uh4*uh6;
		double t15 = t9*uh3*uh8;
		double t16 = t15-t9*uh4*uh7;
		double t17 = nl1*uh3*2.0;
		double t18 = nl2*uh2*2.0;
		double t19 = nl3*uh4*2.0;
		double t20 = param1*uh2*uh3*2.0;
		double t21 = uh1*uh6*uh7*2.0;
		double t22 = t20+t21-uh2*uh3*2.0;
		double t23 = t3*uh1*2.0;
		double t24 = param1*uh1*uh5*2.0;
		double t25 = nl3*t8*uh8;
		double t26 = nl1*uh4*2.0;
		double t27 = nl3*uh2*2.0;
		double t28 = nl2*uh4*2.0;
		double t29 = nl3*uh3*2.0;
		double t30 = nl1*uh2*2.0;
		double t31 = nl2*uh3*2.0;
		double t32 = param1*uh2*uh4*2.0;
		double t33 = uh1*uh6*uh8*2.0;
		double t34 = t32+t33-uh2*uh4*2.0;
		double t35 = param1*uh3*uh4*2.0;
		double t36 = uh1*uh7*uh8*2.0;
		double t37 = t35+t36-uh3*uh4*2.0;
		double t38 = t7*uh1*2.0;
		double t39 = t2*uh1*2.0;
		double t40 = nl1*t8*uh6;
		double t41 = nl2*t8*uh7;
		double t42 = nl3*uh1*uh8*2.0;
		double t43 = uh1*uh4*uh8*2.0;
		double t44 = param2*param2;
		double t45 = nl1*uh1*uh6*2.0;
		double t46 = nl2*uh1*uh7*2.0;
		double t47 = uh1*uh2*uh6*2.0;
		double t48 = uh1*uh3*uh7*2.0;
		fh_uh[0*ng+i] = -param2;
		fh_uh[1*ng+i] = t9*(nl1*t4*-3.0-nl1*t5-nl1*t6+nl1*param1*t4+nl1*param1*t5+nl1*param1*t6-nl1*t2*uh1*2.0-nl1*t3*uh1*2.0+nl1*uh1*uh5*2.0-nl2*uh2*uh3*2.0-nl3*uh2*uh4*2.0-param2*u2*uh1*2.0+param2*uh1*uh2*2.0+nl1*param1*t2*uh1+nl1*param1*t3*uh1+nl1*param1*t7*uh1-nl1*param1*uh1*uh5*2.0+nl2*uh1*uh6*uh7*2.0+nl3*uh1*uh6*uh8*2.0)*(1.0/2.0)-t8*(nl1*t2*-2.0-nl1*t3*2.0+nl1*uh5*2.0-param2*u2*2.0+param2*uh2*2.0+nl1*param1*t2+nl1*param1*t3+nl1*param1*t7-nl1*param1*uh5*2.0+nl2*uh6*uh7*2.0+nl3*uh6*uh8*2.0)*(1.0/2.0);
		fh_uh[2*ng+i] = t9*(-nl2*t4-nl2*t5*3.0-nl2*t6+nl2*param1*t4+nl2*param1*t5+nl2*param1*t6-nl2*t3*uh1*2.0-nl2*t7*uh1*2.0-nl1*uh2*uh3*2.0+nl2*uh1*uh5*2.0-nl3*uh3*uh4*2.0-param2*u3*uh1*2.0+param2*uh1*uh3*2.0+nl2*param1*t2*uh1+nl2*param1*t3*uh1+nl2*param1*t7*uh1-nl2*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh7*2.0+nl3*uh1*uh7*uh8*2.0)*(1.0/2.0)-t8*(nl2*t3*-2.0-nl2*t7*2.0+nl2*uh5*2.0-param2*u3*2.0+param2*uh3*2.0+nl2*param1*t2+nl2*param1*t3+nl2*param1*t7-nl2*param1*uh5*2.0+nl1*uh6*uh7*2.0+nl3*uh7*uh8*2.0)*(1.0/2.0);
		fh_uh[3*ng+i] = t9*(-nl3*t4-nl3*t5-nl3*t6*3.0+nl3*param1*t4+nl3*param1*t5+nl3*param1*t6-nl3*t2*uh1*2.0-nl3*t7*uh1*2.0-nl1*uh2*uh4*2.0-nl2*uh3*uh4*2.0+nl3*uh1*uh5*2.0-param2*u4*uh1*2.0+param2*uh1*uh4*2.0+nl3*param1*t2*uh1+nl3*param1*t3*uh1+nl3*param1*t7*uh1-nl3*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh8*2.0+nl2*uh1*uh7*uh8*2.0)*(1.0/2.0)-t8*(nl3*t2*-2.0-nl3*t7*2.0+nl3*uh5*2.0-param2*u4*2.0+param2*uh4*2.0+nl3*param1*t2+nl3*param1*t3+nl3*param1*t7-nl3*param1*uh5*2.0+nl1*uh6*uh8*2.0+nl2*uh7*uh8*2.0)*(1.0/2.0);
		fh_uh[4*ng+i] = nl1*t10*(-t4*uh2-t5*uh2-t6*uh2+param1*t4*uh2+param1*t5*uh2+param1*t6*uh2-t2*uh1*uh2*2.0-t3*uh1*uh2*2.0+param1*t2*uh1*uh2+param1*t3*uh1*uh2+param1*t7*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)+nl2*t10*(-t4*uh3-t5*uh3-t6*uh3+param1*t4*uh3+param1*t5*uh3+param1*t6*uh3-t3*uh1*uh3*2.0-t7*uh1*uh3*2.0+param1*t2*uh1*uh3+param1*t3*uh1*uh3+param1*t7*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)+nl3*t10*(-t4*uh4-t5*uh4-t6*uh4+param1*t4*uh4+param1*t5*uh4+param1*t6*uh4-t2*uh1*uh4*2.0-t7*uh1*uh4*2.0+param1*t2*uh1*uh4+param1*t3*uh1*uh4+param1*t7*uh1*uh4-param1*uh1*uh4*uh5*2.0+uh1*uh2*uh6*uh8*2.0+uh1*uh3*uh7*uh8*2.0)-nl1*t9*(t2*uh2*-2.0-t3*uh2*2.0+param1*t2*uh2+param1*t3*uh2+param1*t7*uh2-param1*uh2*uh5*2.0+uh3*uh6*uh7*2.0+uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t9*(t3*uh3*-2.0-t7*uh3*2.0+param1*t2*uh3+param1*t3*uh3+param1*t7*uh3-param1*uh3*uh5*2.0+uh2*uh6*uh7*2.0+uh4*uh7*uh8*2.0)*(1.0/2.0)-nl3*t9*(t2*uh4*-2.0-t7*uh4*2.0+param1*t2*uh4+param1*t3*uh4+param1*t7*uh4-param1*uh4*uh5*2.0+uh2*uh6*uh8*2.0+uh3*uh7*uh8*2.0)*(1.0/2.0);
		fh_uh[5*ng+i] = -nl2*t12-nl3*t14;
		fh_uh[6*ng+i] = nl1*t12-nl3*t16;
		fh_uh[7*ng+i] = nl1*t14+nl2*t16;
		fh_uh[8*ng+i] = 0.0;
		fh_uh[9*ng+i] = nl1;
		fh_uh[10*ng+i] = t8*(t19+t31+nl1*uh2*6.0-param2*uh1*2.0-nl1*param1*uh2*2.0)*(1.0/2.0);
		fh_uh[11*ng+i] = t8*(t17+t18-nl2*param1*uh2*2.0)*(1.0/2.0);
		fh_uh[12*ng+i] = t8*(t26+t27-nl3*param1*uh2*2.0)*(1.0/2.0);
		fh_uh[13*ng+i] = nl1*t9*(t4*3.0+t5+t6+t23+t24+t39-param1*t4*3.0-param1*t5-param1*t6-param1*t2*uh1-param1*t3*uh1-param1*t7*uh1)*(1.0/2.0)-nl2*t9*t22*(1.0/2.0)-nl3*t9*t34*(1.0/2.0);
		fh_uh[14*ng+i] = t25+t41;
		fh_uh[15*ng+i] = -nl1*t8*uh7;
		fh_uh[16*ng+i] = -nl1*t8*uh8;
		fh_uh[17*ng+i] = 0.0;
		fh_uh[18*ng+i] = nl2;
		fh_uh[19*ng+i] = t8*(t17+t18-nl1*param1*uh3*2.0)*(1.0/2.0);
		fh_uh[20*ng+i] = t8*(t19+t30+nl2*uh3*6.0-param2*uh1*2.0-nl2*param1*uh3*2.0)*(1.0/2.0);
		fh_uh[21*ng+i] = t8*(t28+t29-nl3*param1*uh3*2.0)*(1.0/2.0);
		fh_uh[22*ng+i] = nl2*t9*(t4+t5*3.0+t6+t23+t24+t38-param1*t4-param1*t5*3.0-param1*t6-param1*t2*uh1-param1*t3*uh1-param1*t7*uh1)*(1.0/2.0)-nl1*t9*t22*(1.0/2.0)-nl3*t9*t37*(1.0/2.0);
		fh_uh[23*ng+i] = -nl2*t8*uh6;
		fh_uh[24*ng+i] = t25+t40;
		fh_uh[25*ng+i] = -nl2*t8*uh8;
		fh_uh[26*ng+i] = 0.0;
		fh_uh[27*ng+i] = nl3;
		fh_uh[28*ng+i] = t8*(t26+t27-nl1*param1*uh4*2.0)*(1.0/2.0);
		fh_uh[29*ng+i] = t8*(t28+t29-nl2*param1*uh4*2.0)*(1.0/2.0);
		fh_uh[30*ng+i] = t8*(t30+t31+nl3*uh4*6.0-param2*uh1*2.0-nl3*param1*uh4*2.0)*(1.0/2.0);
		fh_uh[31*ng+i] = nl3*t9*(t4+t5+t6*3.0+t24+t38+t39-param1*t4-param1*t5-param1*t6*3.0-param1*t2*uh1-param1*t3*uh1-param1*t7*uh1)*(1.0/2.0)-nl1*t9*t34*(1.0/2.0)-nl2*t9*t37*(1.0/2.0);
		fh_uh[32*ng+i] = -nl3*t8*uh6;
		fh_uh[33*ng+i] = -nl3*t8*uh7;
		fh_uh[34*ng+i] = t40+t41;
		fh_uh[35*ng+i] = 0.0;
		fh_uh[36*ng+i] = 0.0;
		fh_uh[37*ng+i] = t8*(nl1*uh1*2.0-nl1*param1*uh1*2.0)*(-1.0/2.0);
		fh_uh[38*ng+i] = t8*(nl2*uh1*2.0-nl2*param1*uh1*2.0)*(-1.0/2.0);
		fh_uh[39*ng+i] = t8*(nl3*uh1*2.0-nl3*param1*uh1*2.0)*(-1.0/2.0);
		fh_uh[40*ng+i] = -param2+nl1*param1*t8*uh2+nl2*param1*t8*uh3+nl3*param1*t8*uh4;
		fh_uh[41*ng+i] = 0.0;
		fh_uh[42*ng+i] = 0.0;
		fh_uh[43*ng+i] = 0.0;
		fh_uh[44*ng+i] = 0.0;
		fh_uh[45*ng+i] = 0.0;
		fh_uh[46*ng+i] = t8*(t42+t46+nl1*param1*uh1*uh6*2.0)*(-1.0/2.0);
		fh_uh[47*ng+i] = t8*(nl1*uh1*uh7*2.0-nl2*uh1*uh6*4.0+nl2*param1*uh1*uh6*2.0)*(-1.0/2.0);
		fh_uh[48*ng+i] = t8*(nl1*uh1*uh8*2.0-nl3*uh1*uh6*4.0+nl3*param1*uh1*uh6*2.0)*(-1.0/2.0);
		fh_uh[49*ng+i] = nl2*t9*(uh1*uh2*uh7*2.0-uh1*uh3*uh6*4.0+param1*uh1*uh3*uh6*2.0)*(-1.0/2.0)-nl3*t9*(uh1*uh2*uh8*2.0-uh1*uh4*uh6*4.0+param1*uh1*uh4*uh6*2.0)*(1.0/2.0)-nl1*t9*(t43+t48+param1*uh1*uh2*uh6*2.0)*(1.0/2.0);
		fh_uh[50*ng+i] = -param2-nl2*t8*uh3-nl3*t8*uh4;
		fh_uh[51*ng+i] = nl1*t8*uh3;
		fh_uh[52*ng+i] = nl1*t8*uh4;
		fh_uh[53*ng+i] = nl1*t44;
		fh_uh[54*ng+i] = 0.0;
		fh_uh[55*ng+i] = t8*(nl1*uh1*uh7*-4.0+nl2*uh1*uh6*2.0+nl1*param1*uh1*uh7*2.0)*(-1.0/2.0);
		fh_uh[56*ng+i] = t8*(t42+t45+nl2*param1*uh1*uh7*2.0)*(-1.0/2.0);
		fh_uh[57*ng+i] = t8*(nl2*uh1*uh8*2.0-nl3*uh1*uh7*4.0+nl3*param1*uh1*uh7*2.0)*(-1.0/2.0);
		fh_uh[58*ng+i] = nl1*t9*(uh1*uh2*uh7*-4.0+uh1*uh3*uh6*2.0+param1*uh1*uh2*uh7*2.0)*(-1.0/2.0)-nl3*t9*(uh1*uh3*uh8*2.0-uh1*uh4*uh7*4.0+param1*uh1*uh4*uh7*2.0)*(1.0/2.0)-nl2*t9*(t43+t47+param1*uh1*uh3*uh7*2.0)*(1.0/2.0);
		fh_uh[59*ng+i] = nl2*t8*uh2;
		fh_uh[60*ng+i] = -param2-nl1*t8*uh2-nl3*t8*uh4;
		fh_uh[61*ng+i] = nl2*t8*uh4;
		fh_uh[62*ng+i] = nl2*t44;
		fh_uh[63*ng+i] = 0.0;
		fh_uh[64*ng+i] = t8*(nl1*uh1*uh8*-4.0+nl3*uh1*uh6*2.0+nl1*param1*uh1*uh8*2.0)*(-1.0/2.0);
		fh_uh[65*ng+i] = t8*(nl2*uh1*uh8*-4.0+nl3*uh1*uh7*2.0+nl2*param1*uh1*uh8*2.0)*(-1.0/2.0);
		fh_uh[66*ng+i] = t8*(t45+t46+nl3*param1*uh1*uh8*2.0)*(-1.0/2.0);
		fh_uh[67*ng+i] = nl1*t9*(uh1*uh2*uh8*-4.0+uh1*uh4*uh6*2.0+param1*uh1*uh2*uh8*2.0)*(-1.0/2.0)-nl2*t9*(uh1*uh3*uh8*-4.0+uh1*uh4*uh7*2.0+param1*uh1*uh3*uh8*2.0)*(1.0/2.0)-nl3*t9*(t47+t48+param1*uh1*uh4*uh8*2.0)*(1.0/2.0);
		fh_uh[68*ng+i] = nl3*t8*uh2;
		fh_uh[69*ng+i] = nl3*t8*uh3;
		fh_uh[70*ng+i] = -param2-nl1*t8*uh2-nl2*t8*uh3;
		fh_uh[71*ng+i] = nl3*t44;
		fh_uh[72*ng+i] = 0.0;
		fh_uh[73*ng+i] = 0.0;
		fh_uh[74*ng+i] = 0.0;
		fh_uh[75*ng+i] = 0.0;
		fh_uh[76*ng+i] = 0.0;
		fh_uh[77*ng+i] = nl1;
		fh_uh[78*ng+i] = nl2;
		fh_uh[79*ng+i] = nl3;
		fh_uh[80*ng+i] = -param2;

	}
}

void fhatonly_mhd3d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double nl3 = nl[2*ng+i];

		double t2 = uh2*uh2;
		double t3 = uh3*uh3;
		double t4 = uh4*uh4;
		double t5 = uh7*uh7;
		double t6 = uh8*uh8;
		double t7 = 1.0/uh1;
		double t8 = uh6*uh6;
		double t9 = 1.0/(uh1*uh1);
		double t10 = t7*uh2*uh7;
		double t11 = t10-t7*uh3*uh6;
		double t12 = t7*uh2*uh8;
		double t13 = t12-t7*uh4*uh6;
		double t14 = t7*uh3*uh8;
		double t15 = t14-t7*uh4*uh7;
		double t16 = param2*param2;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param2*(u1-uh1);
		fh[1*ng+i] = t7*(nl1*t2*-3.0-nl1*t3-nl1*t4+nl1*param1*t2+nl1*param1*t3+nl1*param1*t4-nl1*t5*uh1*2.0-nl1*t6*uh1*2.0+nl1*uh1*uh5*2.0-nl2*uh2*uh3*2.0-nl3*uh2*uh4*2.0-param2*u2*uh1*2.0+param2*uh1*uh2*2.0+nl1*param1*t5*uh1+nl1*param1*t6*uh1+nl1*param1*t8*uh1-nl1*param1*uh1*uh5*2.0+nl2*uh1*uh6*uh7*2.0+nl3*uh1*uh6*uh8*2.0)*(-1.0/2.0);
		fh[2*ng+i] = t7*(-nl2*t2-nl2*t3*3.0-nl2*t4+nl2*param1*t2+nl2*param1*t3+nl2*param1*t4-nl2*t6*uh1*2.0-nl2*t8*uh1*2.0-nl1*uh2*uh3*2.0+nl2*uh1*uh5*2.0-nl3*uh3*uh4*2.0-param2*u3*uh1*2.0+param2*uh1*uh3*2.0+nl2*param1*t5*uh1+nl2*param1*t6*uh1+nl2*param1*t8*uh1-nl2*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh7*2.0+nl3*uh1*uh7*uh8*2.0)*(-1.0/2.0);
		fh[3*ng+i] = t7*(-nl3*t2-nl3*t3-nl3*t4*3.0+nl3*param1*t2+nl3*param1*t3+nl3*param1*t4-nl3*t5*uh1*2.0-nl3*t8*uh1*2.0-nl1*uh2*uh4*2.0-nl2*uh3*uh4*2.0+nl3*uh1*uh5*2.0-param2*u4*uh1*2.0+param2*uh1*uh4*2.0+nl3*param1*t5*uh1+nl3*param1*t6*uh1+nl3*param1*t8*uh1-nl3*param1*uh1*uh5*2.0+nl1*uh1*uh6*uh8*2.0+nl2*uh1*uh7*uh8*2.0)*(-1.0/2.0);
		fh[4*ng+i] = param2*(u5-uh5)-nl1*t9*(-t2*uh2-t3*uh2-t4*uh2+param1*t2*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0-t6*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t6*uh1*uh2+param1*t8*uh1*uh2-param1*uh1*uh2*uh5*2.0+uh1*uh3*uh6*uh7*2.0+uh1*uh4*uh6*uh8*2.0)*(1.0/2.0)-nl2*t9*(-t2*uh3-t3*uh3-t4*uh3+param1*t2*uh3+param1*t3*uh3+param1*t4*uh3-t6*uh1*uh3*2.0-t8*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t6*uh1*uh3+param1*t8*uh1*uh3-param1*uh1*uh3*uh5*2.0+uh1*uh2*uh6*uh7*2.0+uh1*uh4*uh7*uh8*2.0)*(1.0/2.0)-nl3*t9*(-t2*uh4-t3*uh4-t4*uh4+param1*t2*uh4+param1*t3*uh4+param1*t4*uh4-t5*uh1*uh4*2.0-t8*uh1*uh4*2.0+param1*t5*uh1*uh4+param1*t6*uh1*uh4+param1*t8*uh1*uh4-param1*uh1*uh4*uh5*2.0+uh1*uh2*uh6*uh8*2.0+uh1*uh3*uh7*uh8*2.0)*(1.0/2.0);
		fh[5*ng+i] = nl2*t11+nl3*t13+nl1*uh9+param2*(u6-uh6);
		fh[6*ng+i] = -nl1*t11+nl3*t15+nl2*uh9+param2*(u7-uh7);
		fh[7*ng+i] = -nl1*t13-nl2*t15+nl3*uh9+param2*(u8-uh8);
		fh[8*ng+i] = param2*(u9-uh9)+nl1*t16*uh6+nl2*t16*uh7+nl3*t16*uh8;

	}
}
