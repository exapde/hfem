void fbou_ledisp3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *ui, double *param, double time, int ib, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double uinf1 = ui[0];
	double uinf2 = ui[1];
	double uinf3 = ui[2];
    double t2, t3, t4, t5, t6, t7, t8, t9;
    
    for (int i = 0; i <ng*ncu*nc; i++)
		fh_udg[i] = 0.0;
    
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
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		if (ib ==1) {
			fh[0*ng+i] = -uh1+uinf1;
			fh[1*ng+i] = -uh2+uinf2;
			fh[2*ng+i] = -uh3+uinf3;
			fh_uh[0*ng+i] = -1.0;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -1.0;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -1.0;
		}
		else if (ib ==2) {
			t2 = u4+u8+u12;
			t3 = param2*t2;
			t4 = u5+u7;
			t5 = u6+u10;
			t6 = u9+u11;
			fh[0*ng+i] = -uinf1+param3*(u1-uh1)+nl1*(t3+param1*u4*2.0)+nl2*param1*t4+nl3*param1*t5;
			fh[1*ng+i] = -uinf2+param3*(u2-uh2)+nl2*(t3+param1*u8*2.0)+nl1*param1*t4+nl3*param1*t6;
			fh[2*ng+i] = -uinf3+param3*(u3-uh3)+nl3*(t3+param1*u12*2.0)+nl1*param1*t5+nl2*param1*t6;
			t2 = nl1*param1;
			t3 = nl2*param1;
			t4 = param1*2.0;
			t5 = param2+t4;
			t6 = nl3*param2;
			t7 = nl3*param1;
			t8 = nl1*param2;
			t9 = nl2*param2;
			fh_udg[0*ng+i] = param3;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = param3;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = param3;
			fh_udg[9*ng+i] = nl1*t5;
			fh_udg[10*ng+i] = t9;
			fh_udg[11*ng+i] = t6;
			fh_udg[12*ng+i] = t3;
			fh_udg[13*ng+i] = t2;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = t7;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = t2;
			fh_udg[18*ng+i] = t3;
			fh_udg[19*ng+i] = t2;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = t8;
			fh_udg[22*ng+i] = nl2*t5;
			fh_udg[23*ng+i] = t6;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = t7;
			fh_udg[26*ng+i] = t3;
			fh_udg[27*ng+i] = t7;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = t2;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = t7;
			fh_udg[32*ng+i] = t3;
			fh_udg[33*ng+i] = t8;
			fh_udg[34*ng+i] = t9;
			fh_udg[35*ng+i] = nl3*t5;
			fh_uh[0*ng+i] = -param3;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -param3;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -param3;
		}
		else if (ib ==3) {
			t2 = u4+u8+u12;
			t3 = param2*t2;
			t4 = u9+u11;
			fh[0*ng+i] = -uh1+uinf1;
			fh[1*ng+i] = -uinf2+param3*(u2-uh2)+nl2*(t3+param1*u8*2.0)+nl1*param1*(u5+u7)+nl3*param1*t4;
			fh[2*ng+i] = -uinf3+param3*(u3-uh3)+nl3*(t3+param1*u12*2.0)+nl1*param1*(u6+u10)+nl2*param1*t4;
			t2 = nl1*param1;
			t3 = nl3*param2;
			t4 = nl3*param1;
			t5 = nl2*param1;
			t6 = nl2*param2;
			t7 = param1*2.0;
			t8 = param2+t7;
			fh_udg[0*ng+i] = 0.0;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = param3;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = param3;
			fh_udg[9*ng+i] = 0.0;
			fh_udg[10*ng+i] = t6;
			fh_udg[11*ng+i] = t3;
			fh_udg[12*ng+i] = 0.0;
			fh_udg[13*ng+i] = t2;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = 0.0;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = t2;
			fh_udg[18*ng+i] = 0.0;
			fh_udg[19*ng+i] = t2;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = 0.0;
			fh_udg[22*ng+i] = nl2*t8;
			fh_udg[23*ng+i] = t3;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = t4;
			fh_udg[26*ng+i] = t5;
			fh_udg[27*ng+i] = 0.0;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = t2;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = t4;
			fh_udg[32*ng+i] = t5;
			fh_udg[33*ng+i] = 0.0;
			fh_udg[34*ng+i] = t6;
			fh_udg[35*ng+i] = nl3*t8;
			fh_uh[0*ng+i] = -1.0;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -param3;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -param3;
		}
		else if (ib ==4) {
			t2 = u4+u8+u12;
			t3 = param2*t2;
			t4 = u6+u10;
			fh[0*ng+i] = -uinf1+param3*(u1-uh1)+nl1*(t3+param1*u4*2.0)+nl2*param1*(u5+u7)+nl3*param1*t4;
			fh[1*ng+i] = -uh2+uinf2;
			fh[2*ng+i] = -uinf3+param3*(u3-uh3)+nl3*(t3+param1*u12*2.0)+nl2*param1*(u9+u11)+nl1*param1*t4;
			t2 = nl2*param1;
			t3 = nl3*param2;
			t4 = nl3*param1;
			t5 = nl1*param1;
			t6 = nl1*param2;
			t7 = param1*2.0;
			t8 = param2+t7;
			fh_udg[0*ng+i] = param3;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = 0.0;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = param3;
			fh_udg[9*ng+i] = nl1*t8;
			fh_udg[10*ng+i] = 0.0;
			fh_udg[11*ng+i] = t3;
			fh_udg[12*ng+i] = t2;
			fh_udg[13*ng+i] = 0.0;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = t4;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = t5;
			fh_udg[18*ng+i] = t2;
			fh_udg[19*ng+i] = 0.0;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = t6;
			fh_udg[22*ng+i] = 0.0;
			fh_udg[23*ng+i] = t3;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = 0.0;
			fh_udg[26*ng+i] = t2;
			fh_udg[27*ng+i] = t4;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = t5;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = 0.0;
			fh_udg[32*ng+i] = t2;
			fh_udg[33*ng+i] = t6;
			fh_udg[34*ng+i] = 0.0;
			fh_udg[35*ng+i] = nl3*t8;
			fh_uh[0*ng+i] = -param3;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -1.0;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -param3;
		}
		else if (ib ==5) {
			t2 = u4+u8+u12;
			t3 = param2*t2;
			t4 = u5+u7;
			fh[0*ng+i] = -uinf1+param3*(u1-uh1)+nl1*(t3+param1*u4*2.0)+nl3*param1*(u6+u10)+nl2*param1*t4;
			fh[1*ng+i] = -uinf2+param3*(u2-uh2)+nl2*(t3+param1*u8*2.0)+nl3*param1*(u9+u11)+nl1*param1*t4;
			fh[2*ng+i] = -uh3+uinf3;
			t2 = nl2*param1;
			t3 = nl1*param1;
			t4 = param1*2.0;
			t5 = param2+t4;
			t6 = nl3*param1;
			t7 = nl1*param2;
			t8 = nl2*param2;
			fh_udg[0*ng+i] = param3;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = param3;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = 0.0;
			fh_udg[9*ng+i] = nl1*t5;
			fh_udg[10*ng+i] = t8;
			fh_udg[11*ng+i] = 0.0;
			fh_udg[12*ng+i] = t2;
			fh_udg[13*ng+i] = t3;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = t6;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = 0.0;
			fh_udg[18*ng+i] = t2;
			fh_udg[19*ng+i] = t3;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = t7;
			fh_udg[22*ng+i] = nl2*t5;
			fh_udg[23*ng+i] = 0.0;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = t6;
			fh_udg[26*ng+i] = 0.0;
			fh_udg[27*ng+i] = t6;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = 0.0;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = t6;
			fh_udg[32*ng+i] = 0.0;
			fh_udg[33*ng+i] = t7;
			fh_udg[34*ng+i] = t8;
			fh_uh[0*ng+i] = -param3;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -param3;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -1.0;
		}
		else if (ib ==6) {
			fh[0*ng+i] = -uh1+uinf1;
			fh[1*ng+i] = -uh2+uinf2;
			fh[2*ng+i] = -uinf3+nl3*(param1*u12*2.0+param2*(u4+u8+u12))+param3*(u3-uh3)+nl1*param1*(u6+u10)+nl2*param1*(u9+u11);
			t2 = nl3*param2;
			t3 = nl1*param1;
			t4 = nl2*param1;
			fh_udg[0*ng+i] = 0.0;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = 0.0;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = param3;
			fh_udg[9*ng+i] = 0.0;
			fh_udg[10*ng+i] = 0.0;
			fh_udg[11*ng+i] = t2;
			fh_udg[12*ng+i] = 0.0;
			fh_udg[13*ng+i] = 0.0;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = 0.0;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = t3;
			fh_udg[18*ng+i] = 0.0;
			fh_udg[19*ng+i] = 0.0;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = 0.0;
			fh_udg[22*ng+i] = 0.0;
			fh_udg[23*ng+i] = t2;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = 0.0;
			fh_udg[26*ng+i] = t4;
			fh_udg[27*ng+i] = 0.0;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = t3;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = 0.0;
			fh_udg[32*ng+i] = t4;
			fh_udg[33*ng+i] = 0.0;
			fh_udg[34*ng+i] = 0.0;
			fh_udg[35*ng+i] = nl3*(param1*2.0+param2);
			fh_uh[0*ng+i] = -1.0;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -1.0;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -param3;
		}
		else if (ib ==7) {
			fh[0*ng+i] = -uinf1+nl1*(param1*u4*2.0+param2*(u4+u8+u12))+param3*(u1-uh1)+nl2*param1*(u5+u7)+nl3*param1*(u6+u10);
			fh[1*ng+i] = -uh2+uinf2;
			fh[2*ng+i] = -uh3+uinf3;
			t2 = nl2*param1;
			t3 = nl3*param1;
			t4 = nl1*param2;
			fh_udg[0*ng+i] = param3;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = 0.0;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = 0.0;
			fh_udg[9*ng+i] = nl1*(param1*2.0+param2);
			fh_udg[10*ng+i] = 0.0;
			fh_udg[11*ng+i] = 0.0;
			fh_udg[12*ng+i] = t2;
			fh_udg[13*ng+i] = 0.0;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = t3;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = 0.0;
			fh_udg[18*ng+i] = t2;
			fh_udg[19*ng+i] = 0.0;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = t4;
			fh_udg[22*ng+i] = 0.0;
			fh_udg[23*ng+i] = 0.0;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = 0.0;
			fh_udg[26*ng+i] = 0.0;
			fh_udg[27*ng+i] = t3;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = 0.0;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = 0.0;
			fh_udg[32*ng+i] = 0.0;
			fh_udg[33*ng+i] = t4;
			fh_uh[0*ng+i] = -param3;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -1.0;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -1.0;
		}
		else if (ib ==8) {
			fh[0*ng+i] = -uh1+uinf1;
			fh[1*ng+i] = -uinf2+nl2*(param1*u8*2.0+param2*(u4+u8+u12))+param3*(u2-uh2)+nl1*param1*(u5+u7)+nl3*param1*(u9+u11);
			fh[2*ng+i] = -uh3+uinf3;
			t2 = nl1*param1;
			t3 = nl3*param1;
			t4 = nl2*param2;
			fh_udg[0*ng+i] = 0.0;
			fh_udg[1*ng+i] = 0.0;
			fh_udg[2*ng+i] = 0.0;
			fh_udg[3*ng+i] = 0.0;
			fh_udg[4*ng+i] = param3;
			fh_udg[5*ng+i] = 0.0;
			fh_udg[6*ng+i] = 0.0;
			fh_udg[7*ng+i] = 0.0;
			fh_udg[8*ng+i] = 0.0;
			fh_udg[9*ng+i] = 0.0;
			fh_udg[10*ng+i] = t4;
			fh_udg[11*ng+i] = 0.0;
			fh_udg[12*ng+i] = 0.0;
			fh_udg[13*ng+i] = t2;
			fh_udg[14*ng+i] = 0.0;
			fh_udg[15*ng+i] = 0.0;
			fh_udg[16*ng+i] = 0.0;
			fh_udg[17*ng+i] = 0.0;
			fh_udg[18*ng+i] = 0.0;
			fh_udg[19*ng+i] = t2;
			fh_udg[20*ng+i] = 0.0;
			fh_udg[21*ng+i] = 0.0;
			fh_udg[22*ng+i] = nl2*(param1*2.0+param2);
			fh_udg[23*ng+i] = 0.0;
			fh_udg[24*ng+i] = 0.0;
			fh_udg[25*ng+i] = t3;
			fh_udg[26*ng+i] = 0.0;
			fh_udg[27*ng+i] = 0.0;
			fh_udg[28*ng+i] = 0.0;
			fh_udg[29*ng+i] = 0.0;
			fh_udg[30*ng+i] = 0.0;
			fh_udg[31*ng+i] = t3;
			fh_udg[32*ng+i] = 0.0;
			fh_udg[33*ng+i] = 0.0;
			fh_udg[34*ng+i] = t4;
			fh_uh[0*ng+i] = -1.0;
			fh_uh[1*ng+i] = 0.0;
			fh_uh[2*ng+i] = 0.0;
			fh_uh[3*ng+i] = 0.0;
			fh_uh[4*ng+i] = -param3;
			fh_uh[5*ng+i] = 0.0;
			fh_uh[6*ng+i] = 0.0;
			fh_uh[7*ng+i] = 0.0;
			fh_uh[8*ng+i] = -1.0;
		}

	}
}

