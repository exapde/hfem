void flux_mhd23di(double *f, double *f_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u8 = udg[7*ng+i];

		double t2 = u2*u2;
		double t3 = 1.0/(u1*u1);
		double t4 = u6*u6;
		double t5 = t4*(1.0/2.0);
		double t6 = u7*u7;
		double t7 = t6*(1.0/2.0);
		double t8 = u8*u8;
		double t9 = t8*(1.0/2.0);
		double t10 = 1.0/u1;
		double t11 = param1-1.0;
		double t12 = t2*t3*(1.0/2.0);
		double t13 = u3*u3;
		double t14 = t3*t13*(1.0/2.0);
		double t15 = u4*u4;
		double t16 = t3*t15*(1.0/2.0);
		double t17 = t12+t14+t16;
		double t18 = t17*u1;
		double t19 = t5+t7+t9+t18-u5;
		double t20 = t10*u2*u3;
		double t21 = t20-u6*u7;
		double t22 = t5+t7+t9;
		double t23 = t10*t22;
		double t24 = t10*u5;
		double t25 = t23+t24-t10*t11*t19;
		double t26 = t10*u2*u6;
		double t27 = t10*u3*u7;
		double t28 = t10*u4*u8;
		double t29 = t26+t27+t28;
		double t30 = t10*u2*u7;
		f[0*ng+i] = u2;
		f[1*ng+i] = -t5+t7+t9+t2*t10-t11*t19;
		f[2*ng+i] = t21;
		f[3*ng+i] = -u6*u8+t10*u2*u4;
		f[4*ng+i] = t25*u2-t29*u6;
		f[5*ng+i] = 0.0;
		f[6*ng+i] = t30-t10*u3*u6;
		f[7*ng+i] = t10*u2*u8-t10*u4*u6;
		f[8*ng+i] = u3;
		f[9*ng+i] = t21;
		f[10*ng+i] = t5-t7+t9+t10*t13-t11*t19;
		f[11*ng+i] = -u7*u8+t10*u3*u4;
		f[12*ng+i] = t25*u3-t29*u7;
		f[13*ng+i] = -t30+t10*u3*u6;
		f[14*ng+i] = 0.0;
		f[15*ng+i] = t10*u3*u8-t10*u4*u7;

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

		double t2 = 1.0/(u1*u1);
		double t3 = u2*u2;
		double t4 = 1.0/(u1*u1*u1);
		double t5 = u3*u3;
		double t6 = u4*u4;
		double t7 = param1-1.0;
		double t8 = t2*t3*(1.0/2.0);
		double t9 = t2*t5*(1.0/2.0);
		double t10 = t2*t6*(1.0/2.0);
		double t11 = u6*u6;
		double t12 = t11*(1.0/2.0);
		double t13 = u7*u7;
		double t14 = t13*(1.0/2.0);
		double t15 = u8*u8;
		double t16 = t15*(1.0/2.0);
		double t17 = t3*t4;
		double t18 = t4*t5;
		double t19 = t4*t6;
		double t20 = t17+t18+t19;
		double t22 = t20*u1;
		double t21 = t8+t9+t10-t22;
		double t23 = t2*u2*u6;
		double t24 = t2*u3*u7;
		double t25 = t2*u4*u8;
		double t26 = t23+t24+t25;
		double t27 = t12+t14+t16;
		double t28 = t2*t27;
		double t29 = t2*u5;
		double t30 = t8+t9+t10;
		double t31 = t30*u1;
		double t32 = t12+t14+t16+t31-u5;
		double t33 = 1.0/u1;
		double t34 = t7*t21*t33;
		double t35 = t28+t29+t34-t2*t7*t32;
		double t36 = t2*u3*u6;
		double t37 = t33*u3;
		double t38 = t33*u7;
		double t39 = -t33*u6*u7-t2*t7*u2*u3;
		double t40 = t33*u2;
		double t41 = t33*u4;
		double t42 = t27*t33;
		double t43 = t33*u5;
		double t44 = t33*u8;
		double t45 = t33*u6;
		double t46 = t7*t33;
		double t47 = t33+t46;
		double t48 = t45-t7*t33*u6;
		double t49 = t38-t7*t33*u7;
		double t50 = u8-t7*u8;
		double t51 = t44-t7*t33*u8;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = -t2*t3-t7*t21;
		f_udg[2*ng+i] = -t2*u2*u3;
		f_udg[3*ng+i] = -t2*u2*u4;
		f_udg[4*ng+i] = t26*u6-t35*u2;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = t36-t2*u2*u7;
		f_udg[7*ng+i] = -t2*u2*u8+t2*u4*u6;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = -t2*u2*u3;
		f_udg[10*ng+i] = -t2*t5-t7*t21;
		f_udg[11*ng+i] = -t2*u3*u4;
		f_udg[12*ng+i] = t26*u7-t35*u3;
		f_udg[13*ng+i] = -t36+t2*u2*u7;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = -t2*u3*u8+t2*u4*u7;
		f_udg[16*ng+i] = 1.0;
		f_udg[17*ng+i] = t33*u2*2.0-t7*t33*u2;
		f_udg[18*ng+i] = t37;
		f_udg[19*ng+i] = t41;
		f_udg[20*ng+i] = t42+t43-t11*t33-t2*t3*t7-t7*t32*t33;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = t38;
		f_udg[23*ng+i] = t44;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t37;
		f_udg[26*ng+i] = -t7*t33*u2;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = t39;
		f_udg[29*ng+i] = -t38;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = -t7*t33*u3;
		f_udg[34*ng+i] = t40;
		f_udg[35*ng+i] = 0.0;
		f_udg[36*ng+i] = t39;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = -t33*u6;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 1.0;
		f_udg[41*ng+i] = t40;
		f_udg[42*ng+i] = t33*u3*2.0-t7*t33*u3;
		f_udg[43*ng+i] = t41;
		f_udg[44*ng+i] = t42+t43-t13*t33-t2*t5*t7-t7*t32*t33;
		f_udg[45*ng+i] = t45;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = t44;
		f_udg[48*ng+i] = 0.0;
		f_udg[49*ng+i] = -t7*t33*u4;
		f_udg[50*ng+i] = 0.0;
		f_udg[51*ng+i] = t40;
		f_udg[52*ng+i] = -t33*u6*u8-t2*t7*u2*u4;
		f_udg[53*ng+i] = 0.0;
		f_udg[54*ng+i] = 0.0;
		f_udg[55*ng+i] = -t45;
		f_udg[56*ng+i] = 0.0;
		f_udg[57*ng+i] = 0.0;
		f_udg[58*ng+i] = -t7*t33*u4;
		f_udg[59*ng+i] = t37;
		f_udg[60*ng+i] = -t33*u7*u8-t2*t7*u3*u4;
		f_udg[61*ng+i] = 0.0;
		f_udg[62*ng+i] = 0.0;
		f_udg[63*ng+i] = -t38;
		f_udg[64*ng+i] = 0.0;
		f_udg[65*ng+i] = t7;
		f_udg[66*ng+i] = 0.0;
		f_udg[67*ng+i] = 0.0;
		f_udg[68*ng+i] = t47*u2;
		f_udg[69*ng+i] = 0.0;
		f_udg[70*ng+i] = 0.0;
		f_udg[71*ng+i] = 0.0;
		f_udg[72*ng+i] = 0.0;
		f_udg[73*ng+i] = 0.0;
		f_udg[74*ng+i] = t7;
		f_udg[75*ng+i] = 0.0;
		f_udg[76*ng+i] = t47*u3;
		f_udg[77*ng+i] = 0.0;
		f_udg[78*ng+i] = 0.0;
		f_udg[79*ng+i] = 0.0;
		f_udg[80*ng+i] = 0.0;
		f_udg[81*ng+i] = -u6-t7*u6;
		f_udg[82*ng+i] = -u7;
		f_udg[83*ng+i] = -u8;
		f_udg[84*ng+i] = t48*u2-t33*u2*u6*2.0-t33*u3*u7-t33*u4*u8;
		f_udg[85*ng+i] = 0.0;
		f_udg[86*ng+i] = -t37;
		f_udg[87*ng+i] = -t41;
		f_udg[88*ng+i] = 0.0;
		f_udg[89*ng+i] = -u7;
		f_udg[90*ng+i] = u6-t7*u6;
		f_udg[91*ng+i] = 0.0;
		f_udg[92*ng+i] = t48*u3-t33*u2*u7;
		f_udg[93*ng+i] = t37;
		f_udg[94*ng+i] = 0.0;
		f_udg[95*ng+i] = 0.0;
		f_udg[96*ng+i] = 0.0;
		f_udg[97*ng+i] = u7-t7*u7;
		f_udg[98*ng+i] = -u6;
		f_udg[99*ng+i] = 0.0;
		f_udg[100*ng+i] = t49*u2-t33*u3*u6;
		f_udg[101*ng+i] = 0.0;
		f_udg[102*ng+i] = t40;
		f_udg[103*ng+i] = 0.0;
		f_udg[104*ng+i] = 0.0;
		f_udg[105*ng+i] = -u6;
		f_udg[106*ng+i] = -u7-t7*u7;
		f_udg[107*ng+i] = -u8;
		f_udg[108*ng+i] = t49*u3-t33*u2*u6-t33*u3*u7*2.0-t33*u4*u8;
		f_udg[109*ng+i] = -t40;
		f_udg[110*ng+i] = 0.0;
		f_udg[111*ng+i] = -t41;
		f_udg[112*ng+i] = 0.0;
		f_udg[113*ng+i] = t50;
		f_udg[114*ng+i] = 0.0;
		f_udg[115*ng+i] = -u6;
		f_udg[116*ng+i] = t51*u2-t33*u4*u6;
		f_udg[117*ng+i] = 0.0;
		f_udg[118*ng+i] = 0.0;
		f_udg[119*ng+i] = t40;
		f_udg[120*ng+i] = 0.0;
		f_udg[121*ng+i] = 0.0;
		f_udg[122*ng+i] = t50;
		f_udg[123*ng+i] = -u7;
		f_udg[124*ng+i] = t51*u3-t33*u4*u7;
		f_udg[125*ng+i] = 0.0;
		f_udg[126*ng+i] = 0.0;
		f_udg[127*ng+i] = t37;

	}
}

void fluxonly_mhd2d(double *f, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u8 = udg[7*ng+i];

		double t2 = u2*u2;
		double t3 = 1.0/(u1*u1);
		double t4 = u6*u6;
		double t5 = t4*(1.0/2.0);
		double t6 = u7*u7;
		double t7 = t6*(1.0/2.0);
		double t8 = u8*u8;
		double t9 = t8*(1.0/2.0);
		double t10 = 1.0/u1;
		double t11 = param1-1.0;
		double t12 = t2*t3*(1.0/2.0);
		double t13 = u3*u3;
		double t14 = t3*t13*(1.0/2.0);
		double t15 = u4*u4;
		double t16 = t3*t15*(1.0/2.0);
		double t17 = t12+t14+t16;
		double t18 = t17*u1;
		double t19 = t5+t7+t9+t18-u5;
		double t20 = t10*u2*u3;
		double t21 = t20-u6*u7;
		double t22 = t5+t7+t9;
		double t23 = t10*t22;
		double t24 = t10*u5;
		double t25 = t23+t24-t10*t11*t19;
		double t26 = t10*u2*u6;
		double t27 = t10*u3*u7;
		double t28 = t10*u4*u8;
		double t29 = t26+t27+t28;
		double t30 = t10*u2*u7;
		f[0*ng+i] = u2;
		f[1*ng+i] = -t5+t7+t9+t2*t10-t11*t19;
		f[2*ng+i] = t21;
		f[3*ng+i] = -u6*u8+t10*u2*u4;
		f[4*ng+i] = t25*u2-t29*u6;
		f[5*ng+i] = 0.0;
		f[6*ng+i] = t30-t10*u3*u6;
		f[7*ng+i] = t10*u2*u8-t10*u4*u6;
		f[8*ng+i] = u3;
		f[9*ng+i] = t21;
		f[10*ng+i] = t5-t7+t9+t10*t13-t11*t19;
		f[11*ng+i] = -u7*u8+t10*u3*u4;
		f[12*ng+i] = t25*u3-t29*u7;
		f[13*ng+i] = -t30+t10*u3*u6;
		f[14*ng+i] = 0.0;
		f[15*ng+i] = t10*u3*u8-t10*u4*u7;

	}
}
