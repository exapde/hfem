
// Written by C. Nguyen and P. Fernandez

void flux_AVprecomputed_ns2d(double *f, double *f_udg, double *pg, double *avField, double *udg, appstruct &app, double *param, double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian, Int AVtype)
{
    // Description: Artificial viscosity fluxes with precomputed AV field.
    
	double rampFactor;
	double gam = param[0];
    double Re = param[2];
    double Pr = param[3];
    double gam1 = gam - 1.0;
    
    double Pr_AV = 0.9;     // kappa^* = gam * r * beta^* / Pr_AV using the beta sensor
    Int wallBlendingFlag = 1;       // 0: No wall blending. 1: Initial approach. 2:  Approach in note.tex
    double avMax = 1.0e10; //1.0e-2;
    double wallDistance, wallBlending;
    
	for (int i = 0; i <ng; i++) {
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

		wallDistance = pg[(nd+1)*ng+i];
        if (wallDistance < 0.0)
            wallDistance = 0.0;
        
        if (wallBlendingFlag == 0)            // No wall blending
            wallBlending = 1.0;
        else if (wallBlendingFlag == 1)       // Initial approach
            wallBlending = 1.0 - 1.0/(10.0*wallDistance+1.0);
        else if (wallBlendingFlag == 2)       // Approach in note.tex
            wallBlending = 1.0 - exp( - sqrt(Re * wallDistance / 40) );
        
		double av = wallBlending * app.rampFactor * min(avField[i], avMax);

		double r = u1;
		double ru = u2;
		double rv = u3;
		double rE = u4;
		double rx = u5;
		double rux = u6;
		double rvx = u7;
		double rEx = u8;
		double ry = u9;
		double ruy = u10;
		double rvy = u11;
		double rEy = u12;

        if (AVtype == 0) {              // Laplacian with total enthalpy gradient for energy flux
            double t2 = 1.0/(r*r);
            double t3 = 1.0/r;
            double t4 = gam*(1.0/2.0);
            double t5 = t4-1.0/2.0;
            double t6 = ru*ru;
            double t7 = rv*rv;
            
            f[0*ng+i] += av*rx;
            f[1*ng+i] += av*rux;
            f[2*ng+i] += av*rvx;
            f[3*ng+i] += av*(gam*rEx-t5*(rx*t2*t6+rx*t2*t7+ru*t3*(rux-ru*rx*t3)*2.0+rv*t3*(rvx-rv*rx*t3)*2.0));
            f[4*ng+i] += av*ry;
            f[5*ng+i] += av*ruy;
            f[6*ng+i] += av*rvy;
            f[7*ng+i] += av*(gam*rEy-t5*(ry*t2*t6+ry*t2*t7+ru*t3*(ruy-ru*ry*t3)*2.0+rv*t3*(rvy-rv*ry*t3)*2.0));
        }
        else if (AVtype == 1) {         // Bulk viscosity only
            double t2 = 1.0/r;
            double t14 = ru*rx*t2;
            double t3 = rux-t14;
            double t4 = t2*t3;
            double t15 = rv*ry*t2;
            double t5 = rvy-t15;
            double t6 = t2*t5;
            double t7 = t4+t6;
            double t8 = 1.0/(r*r);
            double t9 = ru*ru;
            double t10 = t8*t9*(1.0/2.0);
            double t11 = rv*rv;
            double t12 = t8*t11*(1.0/2.0);
            double t13 = t10+t12;
            double t16 = av*r*t7;
            double t17 = 1.0/Pr_AV;
            double t18 = 1.0/gam1;
            double t19 = rE-r*t13;
            
            f[0*ng+i] += 0.0;
            f[1*ng+i] += t16;
            f[2*ng+i] += 0.0;
            f[3*ng+i] += av*ru*t7-av*gam*t2*t17*t18*(gam1*r*(-rEx+rx*t13+r*(ru*t3*t8+rv*t8*(rvx-rv*rx*t2)))+gam1*rx*t19);
            f[4*ng+i] += 0.0;
            f[5*ng+i] += 0.0;
            f[6*ng+i] += t16;
            f[7*ng+i] += av*rv*t7-av*gam*t2*t17*t18*(gam1*r*(-rEy+ry*t13+r*(rv*t5*t8+ru*t8*(ruy-ru*ry*t2)))+gam1*ry*t19);
        }
        else if (AVtype == 2) {         // Bulk viscosity and thermal conductivity with Pr^* = Pr
            double t2 = 1.0/r;
            double t14 = ru*rx*t2;
            double t3 = rux-t14;
            double t4 = t2*t3;
            double t15 = rv*ry*t2;
            double t5 = rvy-t15;
            double t6 = t2*t5;
            double t7 = t4+t6;
            double t8 = 1.0/(r*r);
            double t9 = ru*ru;
            double t10 = t8*t9*(1.0/2.0);
            double t11 = rv*rv;
            double t12 = t8*t11*(1.0/2.0);
            double t13 = t10+t12;
            double t16 = av*r*t7;
            double t17 = 1.0/Pr;
            double t18 = 1.0/gam1;
            double t19 = rE-r*t13;
            
            f[0*ng+i] += 0.0;
            f[1*ng+i] += t16;
            f[2*ng+i] += 0.0;
            f[3*ng+i] += av*ru*t7-av*gam*t2*t17*t18*(gam1*r*(-rEx+rx*t13+r*(ru*t3*t8+rv*t8*(rvx-rv*rx*t2)))+gam1*rx*t19);
            f[4*ng+i] += 0.0;
            f[5*ng+i] += 0.0;
            f[6*ng+i] += t16;
            f[7*ng+i] += av*rv*t7-av*gam*t2*t17*t18*(gam1*r*(-rEy+ry*t13+r*(rv*t5*t8+ru*t8*(ruy-ru*ry*t2)))+gam1*ry*t19);
        }
        else if (AVtype == 3) {             // Physical model with thermal conductivity only
            double t2 = 1.0/(r*r);
            double t3 = ru*ru;
            double t4 = t2*t3*(1.0/2.0);
            double t5 = rv*rv;
            double t6 = t2*t5*(1.0/2.0);
            double t7 = t4+t6;
            double t8 = 1.0/r;
            double t9 = 1.0/gam1;
            double t10 = rE-r*t7;
            
            f[0*ng+i] += 0.0;
            f[1*ng+i] += 0.0;
            f[2*ng+i] += 0.0;
            f[3*ng+i] += -av*t8*t9*(gam1*r*(-rEx+rx*t7+r*(ru*t2*(rux-ru*rx*t8)+rv*t2*(rvx-rv*rx*t8)))+gam1*rx*t10);
            f[4*ng+i] += 0.0;
            f[5*ng+i] += 0.0;
            f[6*ng+i] += 0.0;
            f[7*ng+i] += -av*t8*t9*(gam1*r*(-rEy+ry*t7+r*(ru*t2*(ruy-ru*ry*t8)+rv*t2*(rvy-rv*ry*t8)))+gam1*ry*t10);
        }
        else if (AVtype == 4) {             // Physical model with molecular viscosity only
            double t2 = 1.0/r;
            double t15 = ru*rx*t2;
            double t3 = rux-t15;
            double t4 = t2*t3*(4.0/3.0);
            double t16 = rv*ry*t2;
            double t5 = rvy-t16;
            double t6 = t4-t2*t5*(2.0/3.0);
            double t12 = ru*ry*t2;
            double t7 = ruy-t12;
            double t8 = t2*t7;
            double t13 = rv*rx*t2;
            double t9 = rvx-t13;
            double t10 = t2*t9;
            double t11 = t8+t10;
            double t14 = av*r*t11;
            double t17 = t2*t3*(2.0/3.0);
            double t18 = t17-t2*t5*(4.0/3.0);

            f[0*ng+i] += 0.0;
            f[1*ng+i] += av*r*t6;
            f[2*ng+i] += t14;
            f[3*ng+i] += av*ru*t6+av*rv*t11;
            f[4*ng+i] += 0.0;
            f[5*ng+i] += t14;
            f[6*ng+i] += -av*r*t18;
            f[7*ng+i] += av*ru*t11-av*rv*t18;
        }
	}
    
	if (computeJacobian == 1) {
		for (int i = 0; i <ng; i++) {
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
            
			wallDistance = pg[(nd+1)*ng+i];
            if (wallDistance < 0.0)
                wallDistance = 0.0;
            
            if (wallBlendingFlag == 0)            // No wall blending
                wallBlending = 1.0;
            else if (wallBlendingFlag == 1)       // Initial approach
                wallBlending = 1.0 - 1.0/(10.0*wallDistance+1.0);
            else if (wallBlendingFlag == 2)       // Approach in note.tex
                wallBlending = 1.0 - exp( - sqrt(Re * wallDistance / 40) );
            
            double av = wallBlending * app.rampFactor * min(avField[i], avMax);

			double r = u1;
			double ru = u2;
			double rv = u3;
			double rE = u4;
			double rx = u5;
			double rux = u6;
			double rvx = u7;
			double rEx = u8;
			double ry = u9;
			double ruy = u10;
			double rvy = u11;
			double rEy = u12;

            if (AVtype == 0) {            // Laplacian with total enthalpy gradient for energy flux
                double t2 = 1.0/(r*r);
                double t3 = 1.0/r;
                double t4 = gam*(1.0/2.0);
                double t5 = t4-1.0/2.0;
                double t6 = rux-ru*rx*t3;
                double t7 = ruy-ru*ry*t3;
                double t8 = rvx-rv*rx*t3;
                double t9 = rvy-rv*ry*t3;
                double t10 = ru*ru;
                double t11 = t2*t10;
                double t12 = rv*rv;
                double t13 = t2*t12;
                double t14 = t11+t13;
                double t15 = av*t5*t14;
                double t16 = av*gam;
                
                f_udg[3*ng+i] += av*t5*(ru*t2*t6*2.0+rv*t2*t8*2.0);
                f_udg[7*ng+i] += av*t5*(ru*t2*t7*2.0+rv*t2*t9*2.0);
                f_udg[11*ng+i] += av*t3*t5*t6*-2.0;
                f_udg[15*ng+i] += av*t3*t5*t7*-2.0;
                f_udg[19*ng+i] += av*t3*t5*t8*-2.0;
                f_udg[23*ng+i] += av*t3*t5*t9*-2.0;
                f_udg[32*ng+i] += av;
                f_udg[35*ng+i] += t15;
                f_udg[41*ng+i] += av;
                f_udg[43*ng+i] += av*ru*t3*t5*-2.0;
                f_udg[50*ng+i] += av;
                f_udg[51*ng+i] += av*rv*t3*t5*-2.0;
                f_udg[59*ng+i] += t16;
                f_udg[68*ng+i] += av;
                f_udg[71*ng+i] += t15;
                f_udg[77*ng+i] += av;
                f_udg[79*ng+i] += av*ru*t3*t5*-2.0;
                f_udg[86*ng+i] += av;
                f_udg[87*ng+i] += av*rv*t3*t5*-2.0;
                f_udg[95*ng+i] += t16;
            }
            else if (AVtype == 1) {
                double t2 = 1.0/r;
                double t7 = ru*rx*t2;
                double t3 = rux-t7;
                double t4 = 1.0/(r*r);
                double t9 = rv*ry*t2;
                double t5 = rvy-t9;
                double t6 = 1.0/(r*r*r);
                double t8 = t3*t4;
                double t10 = t4*t5;
                double t35 = ru*rx*t6;
                double t36 = rv*ry*t6;
                double t11 = t8+t10-t35-t36;
                double t12 = ru*ru;
                double t13 = 1.0/(r*r*r*r);
                double t14 = rv*rv;
                double t17 = rv*rx*t2;
                double t15 = rvx-t17;
                double t16 = ru*t3*t4;
                double t18 = rv*t4*t15;
                double t19 = t4*t12*(1.0/2.0);
                double t20 = t4*t14*(1.0/2.0);
                double t21 = t6*t12;
                double t22 = t6*t14;
                double t23 = t21+t22;
                double t24 = 1.0/Pr_AV;
                double t25 = 1.0/gam1;
                double t26 = t19+t20;
                double t27 = rx*t26;
                double t28 = t16+t18;
                double t29 = r*t28;
                double t30 = -rEx+t27+t29;
                double t31 = t2*t3;
                double t32 = t2*t5;
                double t33 = t31+t32;
                double t34 = av*t33;
                double t37 = t34-av*r*t11;
                double t39 = ru*ry*t2;
                double t38 = ruy-t39;
                double t40 = ru*t4*t38;
                double t41 = rv*t4*t5;
                double t49 = r*t23;
                double t42 = t19+t20-t49;
                double t48 = r*t26;
                double t43 = rE-t48;
                double t44 = ry*t26;
                double t45 = t40+t41;
                double t46 = r*t45;
                double t47 = -rEy+t44+t46;
                double t50 = av*rv*t2;
                double t51 = gam1*t43;
                double t52 = gam1*r*t42;
                double t53 = t51+t52;
                double t54 = av*ru*t2;
                double t55 = av*gam*t24;
                
                f_udg[0*ng+i] += 0.0;
                f_udg[1*ng+i] += t37;
                f_udg[2*ng+i] += 0.0;
                f_udg[3*ng+i] += -av*ru*t11+av*gam*t4*t24*t25*(gam1*r*t30+gam1*rx*t43)-av*gam*t2*t24*t25*(gam1*t30-gam1*rx*t42+gam1*r*(t16+t18-r*(ru*t3*t6*2.0+rv*t6*t15*2.0-rx*t12*t13-rx*t13*t14)-rx*t23));
                f_udg[4*ng+i] += 0.0;
                f_udg[5*ng+i] += 0.0;
                f_udg[6*ng+i] += t37;
                f_udg[7*ng+i] += -av*rv*t11+av*gam*t4*t24*t25*(gam1*r*t47+gam1*ry*t43)-av*gam*t2*t24*t25*(gam1*t47-gam1*ry*t42+gam1*r*(t40+t41-r*(ru*t6*t38*2.0+rv*t5*t6*2.0-ry*t12*t13-ry*t13*t14)-ry*t23));
                f_udg[8*ng+i] += 0.0;
                f_udg[9*ng+i] += -av*rx*t2;
                f_udg[10*ng+i] += 0.0;
                f_udg[11*ng+i] += t34-av*ru*rx*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t8-t35)+ru*rx*t4)-gam1*ru*rx*t2);
                f_udg[12*ng+i] += 0.0;
                f_udg[13*ng+i] += 0.0;
                f_udg[14*ng+i] += -av*rx*t2;
                f_udg[15*ng+i] += -av*rv*rx*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t4*t38-ru*ry*t6)+ru*ry*t4)-gam1*ru*ry*t2);
                f_udg[16*ng+i] += 0.0;
                f_udg[17*ng+i] += -av*ry*t2;
                f_udg[18*ng+i] += 0.0;
                f_udg[19*ng+i] += -av*ru*ry*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t4*t15-rv*rx*t6)+rv*rx*t4)-gam1*rv*rx*t2);
                f_udg[20*ng+i] += 0.0;
                f_udg[21*ng+i] += 0.0;
                f_udg[22*ng+i] += -av*ry*t2;
                f_udg[23*ng+i] += t34-av*rv*ry*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t10-t36)+rv*ry*t4)-gam1*rv*ry*t2);
                f_udg[24*ng+i] += 0.0;
                f_udg[25*ng+i] += 0.0;
                f_udg[26*ng+i] += 0.0;
                f_udg[27*ng+i] += -av*gam*rx*t2*t24;
                f_udg[28*ng+i] += 0.0;
                f_udg[29*ng+i] += 0.0;
                f_udg[30*ng+i] += 0.0;
                f_udg[31*ng+i] += -av*gam*ry*t2*t24;
                f_udg[32*ng+i] += 0.0;
                f_udg[33*ng+i] += -av*ru*t2;
                f_udg[34*ng+i] += 0.0;
                f_udg[35*ng+i] += -av*t4*t12-av*gam*t2*t24*t25*t53;
                f_udg[36*ng+i] += 0.0;
                f_udg[37*ng+i] += 0.0;
                f_udg[38*ng+i] += -av*ru*t2;
                f_udg[39*ng+i] += -av*ru*rv*t4;
                f_udg[40*ng+i] += 0.0;
                f_udg[41*ng+i] += av;
                f_udg[42*ng+i] += 0.0;
                f_udg[43*ng+i] += t54-av*gam*ru*t2*t24;
                f_udg[44*ng+i] += 0.0;
                f_udg[45*ng+i] += 0.0;
                f_udg[46*ng+i] += av;
                f_udg[47*ng+i] += t50;
                f_udg[48*ng+i] += 0.0;
                f_udg[49*ng+i] += 0.0;
                f_udg[50*ng+i] += 0.0;
                f_udg[51*ng+i] += -av*gam*rv*t2*t24;
                f_udg[52*ng+i] += 0.0;
                f_udg[53*ng+i] += 0.0;
                f_udg[54*ng+i] += 0.0;
                f_udg[55*ng+i] += 0.0;
                f_udg[56*ng+i] += 0.0;
                f_udg[57*ng+i] += 0.0;
                f_udg[58*ng+i] += 0.0;
                f_udg[59*ng+i] += t55;
                f_udg[60*ng+i] += 0.0;
                f_udg[61*ng+i] += 0.0;
                f_udg[62*ng+i] += 0.0;
                f_udg[63*ng+i] += 0.0;
                f_udg[64*ng+i] += 0.0;
                f_udg[65*ng+i] += -t50;
                f_udg[66*ng+i] += 0.0;
                f_udg[67*ng+i] += -av*ru*rv*t4;
                f_udg[68*ng+i] += 0.0;
                f_udg[69*ng+i] += 0.0;
                f_udg[70*ng+i] += -t50;
                f_udg[71*ng+i] += -av*t4*t14-av*gam*t2*t24*t25*t53;
                f_udg[72*ng+i] += 0.0;
                f_udg[73*ng+i] += 0.0;
                f_udg[74*ng+i] += 0.0;
                f_udg[75*ng+i] += 0.0;
                f_udg[76*ng+i] += 0.0;
                f_udg[77*ng+i] += 0.0;
                f_udg[78*ng+i] += 0.0;
                f_udg[79*ng+i] += -av*gam*ru*t2*t24;
                f_udg[80*ng+i] += 0.0;
                f_udg[81*ng+i] += av;
                f_udg[82*ng+i] += 0.0;
                f_udg[83*ng+i] += t54;
                f_udg[84*ng+i] += 0.0;
                f_udg[85*ng+i] += 0.0;
                f_udg[86*ng+i] += av;
                f_udg[87*ng+i] += t50-av*gam*rv*t2*t24;
                f_udg[88*ng+i] += 0.0;
                f_udg[89*ng+i] += 0.0;
                f_udg[90*ng+i] += 0.0;
                f_udg[91*ng+i] += 0.0;
                f_udg[92*ng+i] += 0.0;
                f_udg[93*ng+i] += 0.0;
                f_udg[94*ng+i] += 0.0;
                f_udg[95*ng+i] += t55;
            }
            else if (AVtype == 2) {         // Bulk viscosity and thermal conductivity
                double t2 = 1.0/r;
                double t7 = ru*rx*t2;
                double t3 = rux-t7;
                double t4 = 1.0/(r*r);
                double t9 = rv*ry*t2;
                double t5 = rvy-t9;
                double t6 = 1.0/(r*r*r);
                double t8 = t3*t4;
                double t10 = t4*t5;
                double t35 = ru*rx*t6;
                double t36 = rv*ry*t6;
                double t11 = t8+t10-t35-t36;
                double t12 = ru*ru;
                double t13 = 1.0/(r*r*r*r);
                double t14 = rv*rv;
                double t17 = rv*rx*t2;
                double t15 = rvx-t17;
                double t16 = ru*t3*t4;
                double t18 = rv*t4*t15;
                double t19 = t4*t12*(1.0/2.0);
                double t20 = t4*t14*(1.0/2.0);
                double t21 = t6*t12;
                double t22 = t6*t14;
                double t23 = t21+t22;
                double t24 = 1.0/Pr;
                double t25 = 1.0/gam1;
                double t26 = t19+t20;
                double t27 = rx*t26;
                double t28 = t16+t18;
                double t29 = r*t28;
                double t30 = -rEx+t27+t29;
                double t31 = t2*t3;
                double t32 = t2*t5;
                double t33 = t31+t32;
                double t34 = av*t33;
                double t37 = t34-av*r*t11;
                double t39 = ru*ry*t2;
                double t38 = ruy-t39;
                double t40 = ru*t4*t38;
                double t41 = rv*t4*t5;
                double t49 = r*t23;
                double t42 = t19+t20-t49;
                double t48 = r*t26;
                double t43 = rE-t48;
                double t44 = ry*t26;
                double t45 = t40+t41;
                double t46 = r*t45;
                double t47 = -rEy+t44+t46;
                double t50 = av*rv*t2;
                double t51 = gam1*t43;
                double t52 = gam1*r*t42;
                double t53 = t51+t52;
                double t54 = av*ru*t2;
                double t55 = av*gam*t24;
                
                f_udg[0*ng+i] += 0.0;
                f_udg[1*ng+i] += t37;
                f_udg[2*ng+i] += 0.0;
                f_udg[3*ng+i] += -av*ru*t11+av*gam*t4*t24*t25*(gam1*r*t30+gam1*rx*t43)-av*gam*t2*t24*t25*(gam1*t30-gam1*rx*t42+gam1*r*(t16+t18-r*(ru*t3*t6*2.0+rv*t6*t15*2.0-rx*t12*t13-rx*t13*t14)-rx*t23));
                f_udg[4*ng+i] += 0.0;
                f_udg[5*ng+i] += 0.0;
                f_udg[6*ng+i] += t37;
                f_udg[7*ng+i] += -av*rv*t11+av*gam*t4*t24*t25*(gam1*r*t47+gam1*ry*t43)-av*gam*t2*t24*t25*(gam1*t47-gam1*ry*t42+gam1*r*(t40+t41-r*(ru*t6*t38*2.0+rv*t5*t6*2.0-ry*t12*t13-ry*t13*t14)-ry*t23));
                f_udg[8*ng+i] += 0.0;
                f_udg[9*ng+i] += -av*rx*t2;
                f_udg[10*ng+i] += 0.0;
                f_udg[11*ng+i] += t34-av*ru*rx*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t8-t35)+ru*rx*t4)-gam1*ru*rx*t2);
                f_udg[12*ng+i] += 0.0;
                f_udg[13*ng+i] += 0.0;
                f_udg[14*ng+i] += -av*rx*t2;
                f_udg[15*ng+i] += -av*rv*rx*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t4*t38-ru*ry*t6)+ru*ry*t4)-gam1*ru*ry*t2);
                f_udg[16*ng+i] += 0.0;
                f_udg[17*ng+i] += -av*ry*t2;
                f_udg[18*ng+i] += 0.0;
                f_udg[19*ng+i] += -av*ru*ry*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t4*t15-rv*rx*t6)+rv*rx*t4)-gam1*rv*rx*t2);
                f_udg[20*ng+i] += 0.0;
                f_udg[21*ng+i] += 0.0;
                f_udg[22*ng+i] += -av*ry*t2;
                f_udg[23*ng+i] += t34-av*rv*ry*t4-av*gam*t2*t24*t25*(gam1*r*(r*(t10-t36)+rv*ry*t4)-gam1*rv*ry*t2);
                f_udg[24*ng+i] += 0.0;
                f_udg[25*ng+i] += 0.0;
                f_udg[26*ng+i] += 0.0;
                f_udg[27*ng+i] += -av*gam*rx*t2*t24;
                f_udg[28*ng+i] += 0.0;
                f_udg[29*ng+i] += 0.0;
                f_udg[30*ng+i] += 0.0;
                f_udg[31*ng+i] += -av*gam*ry*t2*t24;
                f_udg[32*ng+i] += 0.0;
                f_udg[33*ng+i] += -av*ru*t2;
                f_udg[34*ng+i] += 0.0;
                f_udg[35*ng+i] += -av*t4*t12-av*gam*t2*t24*t25*t53;
                f_udg[36*ng+i] += 0.0;
                f_udg[37*ng+i] += 0.0;
                f_udg[38*ng+i] += -av*ru*t2;
                f_udg[39*ng+i] += -av*ru*rv*t4;
                f_udg[40*ng+i] += 0.0;
                f_udg[41*ng+i] += av;
                f_udg[42*ng+i] += 0.0;
                f_udg[43*ng+i] += t54-av*gam*ru*t2*t24;
                f_udg[44*ng+i] += 0.0;
                f_udg[45*ng+i] += 0.0;
                f_udg[46*ng+i] += av;
                f_udg[47*ng+i] += t50;
                f_udg[48*ng+i] += 0.0;
                f_udg[49*ng+i] += 0.0;
                f_udg[50*ng+i] += 0.0;
                f_udg[51*ng+i] += -av*gam*rv*t2*t24;
                f_udg[52*ng+i] += 0.0;
                f_udg[53*ng+i] += 0.0;
                f_udg[54*ng+i] += 0.0;
                f_udg[55*ng+i] += 0.0;
                f_udg[56*ng+i] += 0.0;
                f_udg[57*ng+i] += 0.0;
                f_udg[58*ng+i] += 0.0;
                f_udg[59*ng+i] += t55;
                f_udg[60*ng+i] += 0.0;
                f_udg[61*ng+i] += 0.0;
                f_udg[62*ng+i] += 0.0;
                f_udg[63*ng+i] += 0.0;
                f_udg[64*ng+i] += 0.0;
                f_udg[65*ng+i] += -t50;
                f_udg[66*ng+i] += 0.0;
                f_udg[67*ng+i] += -av*ru*rv*t4;
                f_udg[68*ng+i] += 0.0;
                f_udg[69*ng+i] += 0.0;
                f_udg[70*ng+i] += -t50;
                f_udg[71*ng+i] += -av*t4*t14-av*gam*t2*t24*t25*t53;
                f_udg[72*ng+i] += 0.0;
                f_udg[73*ng+i] += 0.0;
                f_udg[74*ng+i] += 0.0;
                f_udg[75*ng+i] += 0.0;
                f_udg[76*ng+i] += 0.0;
                f_udg[77*ng+i] += 0.0;
                f_udg[78*ng+i] += 0.0;
                f_udg[79*ng+i] += -av*gam*ru*t2*t24;
                f_udg[80*ng+i] += 0.0;
                f_udg[81*ng+i] += av;
                f_udg[82*ng+i] += 0.0;
                f_udg[83*ng+i] += t54;
                f_udg[84*ng+i] += 0.0;
                f_udg[85*ng+i] += 0.0;
                f_udg[86*ng+i] += av;
                f_udg[87*ng+i] += t50-av*gam*rv*t2*t24;
                f_udg[88*ng+i] += 0.0;
                f_udg[89*ng+i] += 0.0;
                f_udg[90*ng+i] += 0.0;
                f_udg[91*ng+i] += 0.0;
                f_udg[92*ng+i] += 0.0;
                f_udg[93*ng+i] += 0.0;
                f_udg[94*ng+i] += 0.0;
                f_udg[95*ng+i] += t55;
            }
            else if (AVtype == 3) {     // Physical model with thermal conductivity only
                double t2 = 1.0/(r*r);
                double t3 = 1.0/r;
                double t4 = ru*ru;
                double t5 = 1.0/(r*r*r*r);
                double t6 = rv*rv;
                double t10 = ru*rx*t3;
                double t7 = rux-t10;
                double t8 = 1.0/(r*r*r);
                double t12 = rv*rx*t3;
                double t9 = rvx-t12;
                double t11 = ru*t2*t7;
                double t13 = rv*t2*t9;
                double t14 = t2*t4*(1.0/2.0);
                double t15 = t2*t6*(1.0/2.0);
                double t16 = t4*t8;
                double t17 = t6*t8;
                double t18 = t16+t17;
                double t19 = 1.0/gam1;
                double t20 = t14+t15;
                double t21 = rx*t20;
                double t22 = t11+t13;
                double t23 = r*t22;
                double t24 = -rEx+t21+t23;
                double t27 = ru*ry*t3;
                double t25 = ruy-t27;
                double t29 = rv*ry*t3;
                double t26 = rvy-t29;
                double t28 = ru*t2*t25;
                double t30 = rv*t2*t26;
                double t38 = r*t18;
                double t31 = t14+t15-t38;
                double t37 = r*t20;
                double t32 = rE-t37;
                double t33 = ry*t20;
                double t34 = t28+t30;
                double t35 = r*t34;
                double t36 = -rEy+t33+t35;
                double t39 = gam1*t32;
                double t40 = gam1*r*t31;
                double t41 = t39+t40;
                
                f_udg[0*ng+i] += 0.0;
                f_udg[1*ng+i] += 0.0;
                f_udg[2*ng+i] += 0.0;
                f_udg[3*ng+i] += av*t2*t19*(gam1*r*t24+gam1*rx*t32)-av*t3*t19*(gam1*t24-gam1*rx*t31+gam1*r*(t11+t13-r*(ru*t7*t8*2.0+rv*t8*t9*2.0-rx*t4*t5-rx*t5*t6)-rx*t18));
                f_udg[4*ng+i] += 0.0;
                f_udg[5*ng+i] += 0.0;
                f_udg[6*ng+i] += 0.0;
                f_udg[7*ng+i] += av*t2*t19*(gam1*r*t36+gam1*ry*t32)-av*t3*t19*(gam1*t36-gam1*ry*t31+gam1*r*(t28+t30-r*(ru*t8*t25*2.0+rv*t8*t26*2.0-ry*t4*t5-ry*t5*t6)-ry*t18));
                f_udg[8*ng+i] += 0.0;
                f_udg[9*ng+i] += 0.0;
                f_udg[10*ng+i] += 0.0;
                f_udg[11*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t7-ru*rx*t8)+ru*rx*t2)-gam1*ru*rx*t3);
                f_udg[12*ng+i] += 0.0;
                f_udg[13*ng+i] += 0.0;
                f_udg[14*ng+i] += 0.0;
                f_udg[15*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t25-ru*ry*t8)+ru*ry*t2)-gam1*ru*ry*t3);
                f_udg[16*ng+i] += 0.0;
                f_udg[17*ng+i] += 0.0;
                f_udg[18*ng+i] += 0.0;
                f_udg[19*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t9-rv*rx*t8)+rv*rx*t2)-gam1*rv*rx*t3);
                f_udg[20*ng+i] += 0.0;
                f_udg[21*ng+i] += 0.0;
                f_udg[22*ng+i] += 0.0;
                f_udg[23*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t26-rv*ry*t8)+rv*ry*t2)-gam1*rv*ry*t3);
                f_udg[24*ng+i] += 0.0;
                f_udg[25*ng+i] += 0.0;
                f_udg[26*ng+i] += 0.0;
                f_udg[27*ng+i] += -av*rx*t3;
                f_udg[28*ng+i] += 0.0;
                f_udg[29*ng+i] += 0.0;
                f_udg[30*ng+i] += 0.0;
                f_udg[31*ng+i] += -av*ry*t3;
                f_udg[32*ng+i] += 0.0;
                f_udg[33*ng+i] += 0.0;
                f_udg[34*ng+i] += 0.0;
                f_udg[35*ng+i] += -av*t3*t19*t41;
                f_udg[36*ng+i] += 0.0;
                f_udg[37*ng+i] += 0.0;
                f_udg[38*ng+i] += 0.0;
                f_udg[39*ng+i] += 0.0;
                f_udg[40*ng+i] += 0.0;
                f_udg[41*ng+i] += 0.0;
                f_udg[42*ng+i] += 0.0;
                f_udg[43*ng+i] += -av*ru*t3;
                f_udg[44*ng+i] += 0.0;
                f_udg[45*ng+i] += 0.0;
                f_udg[46*ng+i] += 0.0;
                f_udg[47*ng+i] += 0.0;
                f_udg[48*ng+i] += 0.0;
                f_udg[49*ng+i] += 0.0;
                f_udg[50*ng+i] += 0.0;
                f_udg[51*ng+i] += -av*rv*t3;
                f_udg[52*ng+i] += 0.0;
                f_udg[53*ng+i] += 0.0;
                f_udg[54*ng+i] += 0.0;
                f_udg[55*ng+i] += 0.0;
                f_udg[56*ng+i] += 0.0;
                f_udg[57*ng+i] += 0.0;
                f_udg[58*ng+i] += 0.0;
                f_udg[59*ng+i] += av;
                f_udg[60*ng+i] += 0.0;
                f_udg[61*ng+i] += 0.0;
                f_udg[62*ng+i] += 0.0;
                f_udg[63*ng+i] += 0.0;
                f_udg[64*ng+i] += 0.0;
                f_udg[65*ng+i] += 0.0;
                f_udg[66*ng+i] += 0.0;
                f_udg[67*ng+i] += 0.0;
                f_udg[68*ng+i] += 0.0;
                f_udg[69*ng+i] += 0.0;
                f_udg[70*ng+i] += 0.0;
                f_udg[71*ng+i] += -av*t3*t19*t41;
                f_udg[72*ng+i] += 0.0;
                f_udg[73*ng+i] += 0.0;
                f_udg[74*ng+i] += 0.0;
                f_udg[75*ng+i] += 0.0;
                f_udg[76*ng+i] += 0.0;
                f_udg[77*ng+i] += 0.0;
                f_udg[78*ng+i] += 0.0;
                f_udg[79*ng+i] += -av*ru*t3;
                f_udg[80*ng+i] += 0.0;
                f_udg[81*ng+i] += 0.0;
                f_udg[82*ng+i] += 0.0;
                f_udg[83*ng+i] += 0.0;
                f_udg[84*ng+i] += 0.0;
                f_udg[85*ng+i] += 0.0;
                f_udg[86*ng+i] += 0.0;
                f_udg[87*ng+i] += -av*rv*t3;
                f_udg[88*ng+i] += 0.0;
                f_udg[89*ng+i] += 0.0;
                f_udg[90*ng+i] += 0.0;
                f_udg[91*ng+i] += 0.0;
                f_udg[92*ng+i] += 0.0;
                f_udg[93*ng+i] += 0.0;
                f_udg[94*ng+i] += 0.0;
                f_udg[95*ng+i] += av;
            }
            else if (AVtype == 4) {             // Physical model with molecular viscosity only
                double t2 = 1.0/r;
                double t9 = ru*rx*t2;
                double t3 = rux-t9;
                double t4 = 1.0/(r*r);
                double t11 = rv*ry*t2;
                double t5 = rvy-t11;
                double t6 = 1.0/(r*r*r);
                double t14 = ru*ry*t2;
                double t7 = ruy-t14;
                double t16 = rv*rx*t2;
                double t8 = rvx-t16;
                double t10 = t3*t4*(4.0/3.0);
                double t12 = rv*ry*t6*(2.0/3.0);
                double t13 = t10+t12-t4*t5*(2.0/3.0)-ru*rx*t6*(4.0/3.0);
                double t15 = t4*t7;
                double t17 = t4*t8;
                double t23 = ru*ry*t6;
                double t24 = rv*rx*t6;
                double t18 = t15+t17-t23-t24;
                double t19 = t2*t7;
                double t20 = t2*t8;
                double t21 = t19+t20;
                double t22 = av*t21;
                double t25 = t22-av*r*t18;
                double t26 = t3*t4*(2.0/3.0);
                double t27 = rv*ry*t6*(4.0/3.0);
                double t28 = t26+t27-t4*t5*(4.0/3.0)-ru*rx*t6*(2.0/3.0);
                double t29 = t2*t3*(4.0/3.0);
                double t30 = t29-t2*t5*(2.0/3.0);
                double t31 = av*t30;
                double t32 = t2*t3*(2.0/3.0);
                double t33 = t32-t2*t5*(4.0/3.0);
                double t34 = av*ru*t2;
                double t35 = ru*ru;
                double t36 = rv*rv;
                double t37 = av*rv*t2;
                double t38 = av*ru*t2*(2.0/3.0);
                double t39 = av*(4.0/3.0);
                
                f_udg[0*ng+i] += 0.0;
                f_udg[1*ng+i] += t31-av*r*t13;
                f_udg[2*ng+i] += t25;
                f_udg[3*ng+i] += -av*ru*t13-av*rv*t18;
                f_udg[4*ng+i] += 0.0;
                f_udg[5*ng+i] += t25;
                f_udg[6*ng+i] += -av*t33+av*r*t28;
                f_udg[7*ng+i] += -av*ru*t18+av*rv*t28;
                f_udg[8*ng+i] += 0.0;
                f_udg[9*ng+i] += av*rx*t2*(-4.0/3.0);
                f_udg[10*ng+i] += -av*ry*t2;
                f_udg[11*ng+i] += t31-av*ru*rx*t4*(4.0/3.0)-av*rv*ry*t4;
                f_udg[12*ng+i] += 0.0;
                f_udg[13*ng+i] += -av*ry*t2;
                f_udg[14*ng+i] += av*rx*t2*(2.0/3.0);
                f_udg[15*ng+i] += t22-av*ru*ry*t4+av*rv*rx*t4*(2.0/3.0);
                f_udg[16*ng+i] += 0.0;
                f_udg[17*ng+i] += av*ry*t2*(2.0/3.0);
                f_udg[18*ng+i] += -av*rx*t2;
                f_udg[19*ng+i] += t22+av*ru*ry*t4*(2.0/3.0)-av*rv*rx*t4;
                f_udg[20*ng+i] += 0.0;
                f_udg[21*ng+i] += -av*rx*t2;
                f_udg[22*ng+i] += av*ry*t2*(-4.0/3.0);
                f_udg[23*ng+i] += -av*t33-av*ru*rx*t4-av*rv*ry*t4*(4.0/3.0);
                f_udg[24*ng+i] += 0.0;
                f_udg[25*ng+i] += 0.0;
                f_udg[26*ng+i] += 0.0;
                f_udg[27*ng+i] += 0.0;
                f_udg[28*ng+i] += 0.0;
                f_udg[29*ng+i] += 0.0;
                f_udg[30*ng+i] += 0.0;
                f_udg[31*ng+i] += 0.0;
                f_udg[32*ng+i] += 0.0;
                f_udg[33*ng+i] += av*ru*t2*(-4.0/3.0);
                f_udg[34*ng+i] += -av*rv*t2;
                f_udg[35*ng+i] += av*t4*t35*(-4.0/3.0)-av*t4*t36;
                f_udg[36*ng+i] += 0.0;
                f_udg[37*ng+i] += -av*rv*t2;
                f_udg[38*ng+i] += t38;
                f_udg[39*ng+i] += av*ru*rv*t4*(-1.0/3.0);
                f_udg[40*ng+i] += 0.0;
                f_udg[41*ng+i] += t39;
                f_udg[42*ng+i] += 0.0;
                f_udg[43*ng+i] += av*ru*t2*(4.0/3.0);
                f_udg[44*ng+i] += 0.0;
                f_udg[45*ng+i] += 0.0;
                f_udg[46*ng+i] += av*(-2.0/3.0);
                f_udg[47*ng+i] += av*rv*t2*(-2.0/3.0);
                f_udg[48*ng+i] += 0.0;
                f_udg[49*ng+i] += 0.0;
                f_udg[50*ng+i] += av;
                f_udg[51*ng+i] += t37;
                f_udg[52*ng+i] += 0.0;
                f_udg[53*ng+i] += av;
                f_udg[54*ng+i] += 0.0;
                f_udg[55*ng+i] += t34;
                f_udg[56*ng+i] += 0.0;
                f_udg[57*ng+i] += 0.0;
                f_udg[58*ng+i] += 0.0;
                f_udg[59*ng+i] += 0.0;
                f_udg[60*ng+i] += 0.0;
                f_udg[61*ng+i] += 0.0;
                f_udg[62*ng+i] += 0.0;
                f_udg[63*ng+i] += 0.0;
                f_udg[64*ng+i] += 0.0;
                f_udg[65*ng+i] += av*rv*t2*(2.0/3.0);
                f_udg[66*ng+i] += -t34;
                f_udg[67*ng+i] += av*ru*rv*t4*(-1.0/3.0);
                f_udg[68*ng+i] += 0.0;
                f_udg[69*ng+i] += -t34;
                f_udg[70*ng+i] += av*rv*t2*(-4.0/3.0);
                f_udg[71*ng+i] += -av*t4*t35-av*t4*t36*(4.0/3.0);
                f_udg[72*ng+i] += 0.0;
                f_udg[73*ng+i] += 0.0;
                f_udg[74*ng+i] += av;
                f_udg[75*ng+i] += t37;
                f_udg[76*ng+i] += 0.0;
                f_udg[77*ng+i] += av;
                f_udg[78*ng+i] += 0.0;
                f_udg[79*ng+i] += t34;
                f_udg[80*ng+i] += 0.0;
                f_udg[81*ng+i] += av*(-2.0/3.0);
                f_udg[82*ng+i] += 0.0;
                f_udg[83*ng+i] += -t38;
                f_udg[84*ng+i] += 0.0;
                f_udg[85*ng+i] += 0.0;
                f_udg[86*ng+i] += t39;
                f_udg[87*ng+i] += av*rv*t2*(4.0/3.0);
            }
		}
	}
}
