void source_mhd23dq(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];
		double u13 = udg[12*ng+i];
		double u14 = udg[13*ng+i];
		double u15 = udg[14*ng+i];
		double u16 = udg[15*ng+i];
		double u17 = udg[16*ng+i];
		double u18 = udg[17*ng+i];
		double u19 = udg[18*ng+i];
		double u20 = udg[19*ng+i];
		double u21 = udg[20*ng+i];
		double u22 = udg[21*ng+i];
		double u23 = udg[22*ng+i];
		double u24 = udg[23*ng+i];
		double u25 = udg[24*ng+i];
		double u26 = udg[25*ng+i];
		double u27 = udg[26*ng+i];

		s[0*ng+i] = 0.0;
		s[1*ng+i] = 0.0;
		s[2*ng+i] = 0.0;
		s[3*ng+i] = 0.0;
		s[4*ng+i] = 0.0;
		s[5*ng+i] = 0.0;
		s[6*ng+i] = 0.0;
		s[7*ng+i] = 0.0;
		s[8*ng+i] = u9*-2.0;

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
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];
		double u13 = udg[12*ng+i];
		double u14 = udg[13*ng+i];
		double u15 = udg[14*ng+i];
		double u16 = udg[15*ng+i];
		double u17 = udg[16*ng+i];
		double u18 = udg[17*ng+i];
		double u19 = udg[18*ng+i];
		double u20 = udg[19*ng+i];
		double u21 = udg[20*ng+i];
		double u22 = udg[21*ng+i];
		double u23 = udg[22*ng+i];
		double u24 = udg[23*ng+i];
		double u25 = udg[24*ng+i];
		double u26 = udg[25*ng+i];
		double u27 = udg[26*ng+i];

		s_udg[0*ng+i] = 0.0;
		s_udg[1*ng+i] = 0.0;
		s_udg[2*ng+i] = 0.0;
		s_udg[3*ng+i] = 0.0;
		s_udg[4*ng+i] = 0.0;
		s_udg[5*ng+i] = 0.0;
		s_udg[6*ng+i] = 0.0;
		s_udg[7*ng+i] = 0.0;
		s_udg[8*ng+i] = 0.0;
		s_udg[9*ng+i] = 0.0;
		s_udg[10*ng+i] = 0.0;
		s_udg[11*ng+i] = 0.0;
		s_udg[12*ng+i] = 0.0;
		s_udg[13*ng+i] = 0.0;
		s_udg[14*ng+i] = 0.0;
		s_udg[15*ng+i] = 0.0;
		s_udg[16*ng+i] = 0.0;
		s_udg[17*ng+i] = 0.0;
		s_udg[18*ng+i] = 0.0;
		s_udg[19*ng+i] = 0.0;
		s_udg[20*ng+i] = 0.0;
		s_udg[21*ng+i] = 0.0;
		s_udg[22*ng+i] = 0.0;
		s_udg[23*ng+i] = 0.0;
		s_udg[24*ng+i] = 0.0;
		s_udg[25*ng+i] = 0.0;
		s_udg[26*ng+i] = 0.0;
		s_udg[27*ng+i] = 0.0;
		s_udg[28*ng+i] = 0.0;
		s_udg[29*ng+i] = 0.0;
		s_udg[30*ng+i] = 0.0;
		s_udg[31*ng+i] = 0.0;
		s_udg[32*ng+i] = 0.0;
		s_udg[33*ng+i] = 0.0;
		s_udg[34*ng+i] = 0.0;
		s_udg[35*ng+i] = 0.0;
		s_udg[36*ng+i] = 0.0;
		s_udg[37*ng+i] = 0.0;
		s_udg[38*ng+i] = 0.0;
		s_udg[39*ng+i] = 0.0;
		s_udg[40*ng+i] = 0.0;
		s_udg[41*ng+i] = 0.0;
		s_udg[42*ng+i] = 0.0;
		s_udg[43*ng+i] = 0.0;
		s_udg[44*ng+i] = 0.0;
		s_udg[45*ng+i] = 0.0;
		s_udg[46*ng+i] = 0.0;
		s_udg[47*ng+i] = 0.0;
		s_udg[48*ng+i] = 0.0;
		s_udg[49*ng+i] = 0.0;
		s_udg[50*ng+i] = 0.0;
		s_udg[51*ng+i] = 0.0;
		s_udg[52*ng+i] = 0.0;
		s_udg[53*ng+i] = 0.0;
		s_udg[54*ng+i] = 0.0;
		s_udg[55*ng+i] = 0.0;
		s_udg[56*ng+i] = 0.0;
		s_udg[57*ng+i] = 0.0;
		s_udg[58*ng+i] = 0.0;
		s_udg[59*ng+i] = 0.0;
		s_udg[60*ng+i] = 0.0;
		s_udg[61*ng+i] = 0.0;
		s_udg[62*ng+i] = 0.0;
		s_udg[63*ng+i] = 0.0;
		s_udg[64*ng+i] = 0.0;
		s_udg[65*ng+i] = 0.0;
		s_udg[66*ng+i] = 0.0;
		s_udg[67*ng+i] = 0.0;
		s_udg[68*ng+i] = 0.0;
		s_udg[69*ng+i] = 0.0;
		s_udg[70*ng+i] = 0.0;
		s_udg[71*ng+i] = 0.0;
		s_udg[72*ng+i] = 0.0;
		s_udg[73*ng+i] = 0.0;
		s_udg[74*ng+i] = 0.0;
		s_udg[75*ng+i] = 0.0;
		s_udg[76*ng+i] = 0.0;
		s_udg[77*ng+i] = 0.0;
		s_udg[78*ng+i] = 0.0;
		s_udg[79*ng+i] = 0.0;
		s_udg[80*ng+i] = -2.0;

	}
}

void sourceonly_mhd23dq(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];
		double u13 = udg[12*ng+i];
		double u14 = udg[13*ng+i];
		double u15 = udg[14*ng+i];
		double u16 = udg[15*ng+i];
		double u17 = udg[16*ng+i];
		double u18 = udg[17*ng+i];
		double u19 = udg[18*ng+i];
		double u20 = udg[19*ng+i];
		double u21 = udg[20*ng+i];
		double u22 = udg[21*ng+i];
		double u23 = udg[22*ng+i];
		double u24 = udg[23*ng+i];
		double u25 = udg[24*ng+i];
		double u26 = udg[25*ng+i];
		double u27 = udg[26*ng+i];

		s[0*ng+i] = 0.0;
		s[1*ng+i] = 0.0;
		s[2*ng+i] = 0.0;
		s[3*ng+i] = 0.0;
		s[4*ng+i] = 0.0;
		s[5*ng+i] = 0.0;
		s[6*ng+i] = 0.0;
		s[7*ng+i] = 0.0;
		s[8*ng+i] = u9*-2.0;

	}
}

