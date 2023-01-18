void fhat_mhd2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];
	double param7 = param[6];
	double param8 = param[7];

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
		double t17 = 1.0/u1;
		double t18 = t17*u2*1.0E2;
		double t19 = exp(t18);
		double t20 = t19+1.0;
		double t21 = log(t20);
		double t22 = u5*u5;
		double t23 = u6*u6;
		double t24 = t22+t23;
		double t25 = t17*t24;
		double t27 = param1-1.0;
		double t28 = u1*u4*2.0;
		double t29 = t22*u1;
		double t30 = t23*u1;
		double t31 = u2*u2;
		double t32 = u3*u3;
		double t33 = -t28+t29+t30+t31+t32;
		double t34 = 1.0/(u1*u1);
		double t35 = param1*t27*t33*t34*(1.0/2.0);
		double t26 = t25-t35;
		double t36 = sqrt(2.0);
		double t37 = t26*t26;
		double t38 = 1.0/(u1*u1*u1);
		double t39 = exp(-t18);
		double t40 = t39+1.0;
		double t41 = log(t40);
		double t42 = t17*u3*1.0E2;
		double t43 = exp(t42);
		double t44 = t43+1.0;
		double t45 = log(t44);
		double t46 = param1*t22*t27*t33*t38*2.0;
		double t47 = t37+t46;
		double t48 = sqrt(t47);
		double t49 = t25-t35+t48;
		double t50 = sqrt(t49);
		double t51 = param1*t23*t27*t33*t38*2.0;
		double t52 = t37+t51;
		double t53 = sqrt(t52);
		double t54 = t25-t35+t53;
		double t55 = sqrt(t54);
		double t57 = t21*4.503599627370496E15;
		double t58 = 1.0/3.141592653589793;
		double t59 = t41*1.0E6;
		double t60 = t21*1.0E6;
		double t61 = exp(-t42);
		double t62 = t61+1.0;
		double t63 = log(t62);
		double t64 = t45*1.0E6;
		double t65 = t36*t50*5.0E7;
		double t66 = t36*t55*5.0E7;
		double t67 = t36*t50*5.0E1;
		double t68 = t36*t55*5.0E1;
		double t69 = t21+t41-t45-t63+t67-t68;
		double t70 = t36*t50*2.251799813685248E17;
		double t56 = t57+t70+log(exp(t17*u2*-1.0E2)+1.0)*4.503599627370496E15+t69*(t58*atan(t59+t60-t64+t65-t66-log(exp(t17*u3*-1.0E2)+1.0)*1.0E6)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		double t71 = t41*4.503599627370496E15+t57+t70+t69*(t58*atan(t59+t60-t63*1.0E6-t64+t65-t66)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = param3*(u2-uh2)-nl1*t2*(t3*-3.0-t4+t8+t9+t10+t12+t13-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl2*t2*t7;
		fh[2*ng+i] = param4*(u3-uh3)-nl2*t2*(-t3-t4*3.0+t8+t9+t10+t12+t13-t11*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl1*t2*t7;
		fh[3*ng+i] = param5*(u4-uh4)-nl1*t14*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t11*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t14*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t11*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t11*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh[4*ng+i] = -nl2*t16+nl1*uh7+param6*(u5-uh5);
		fh[5*ng+i] = nl1*t16+nl2*uh7+param7*(u6-uh6);
		fh[6*ng+i] = param8*(u7-uh7)+nl1*(t56*t56)*uh5*4.930380657631324E-36+nl2*(t71*t71)*uh6*4.930380657631324E-36;

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

		double t2 = 1.0/u1;
		double t3 = u5*u5;
		double t4 = u6*u6;
		double t5 = t3+t4;
		double t6 = t2*t5;
		double t8 = param1-1.0;
		double t9 = u1*u4*2.0;
		double t10 = t3*u1;
		double t11 = t4*u1;
		double t12 = u2*u2;
		double t13 = u3*u3;
		double t14 = -t9+t10+t11+t12+t13;
		double t15 = 1.0/(u1*u1);
		double t16 = param1*t8*t14*t15*(1.0/2.0);
		double t7 = t6-t16;
		double t17 = sqrt(2.0);
		double t18 = t7*t7;
		double t19 = 1.0/(u1*u1*u1);
		double t20 = param1*t3*t8*t14*t19*2.0;
		double t21 = t18+t20;
		double t22 = sqrt(t21);
		double t23 = t6-t16+t22;
		double t24 = t5*t15;
		double t26 = u4*2.0;
		double t25 = t3+t4-t26;
		double t27 = param1*t8*t15*t25*(1.0/2.0);
		double t28 = param1*t4*t8*t14*t19*2.0;
		double t29 = t18+t28;
		double t30 = sqrt(t29);
		double t31 = t6-t16+t30;
		double t35 = param1*t8*t14*t19;
		double t32 = t24+t27-t35;
		double t33 = t7*t32*2.0;
		double t34 = 1.0/(u1*u1*u1*u1);
		double t36 = t2*u2*1.0E2;
		double t37 = exp(-t36);
		double t38 = exp(t36);
		double t39 = t38+1.0;
		double t40 = t2*u3*1.0E2;
		double t41 = exp(-t40);
		double t42 = exp(t40);
		double t43 = t42+1.0;
		double t44 = t37+1.0;
		double t45 = t41+1.0;
		double t46 = 1.0/sqrt(t23);
		double t47 = 1.0/sqrt(t21);
		double t48 = param1*t3*t8*t14*t34*6.0;
		double t63 = param1*t3*t8*t19*t25*2.0;
		double t49 = t33+t48-t63;
		double t50 = t47*t49*(1.0/2.0);
		double t51 = t24+t27-t35+t50;
		double t52 = 1.0/3.141592653589793;
		double t53 = log(t44);
		double t54 = t53*1.0E6;
		double t55 = log(t39);
		double t56 = t55*1.0E6;
		double t57 = log(t45);
		double t58 = log(t43);
		double t59 = sqrt(t23);
		double t60 = t17*t59*5.0E7;
		double t61 = sqrt(t31);
		double t74 = t57*1.0E6;
		double t75 = t58*1.0E6;
		double t76 = t17*t61*5.0E7;
		double t62 = t54+t56+t60-t74-t75-t76;
		double t64 = 1.0/sqrt(t31);
		double t65 = 1.0/sqrt(t29);
		double t66 = param1*t4*t8*t14*t34*6.0;
		double t83 = param1*t4*t8*t19*t25*2.0;
		double t67 = t33+t66-t83;
		double t68 = t65*t67*(1.0/2.0);
		double t69 = t24+t27-t35+t68;
		double t70 = 1.0/t44;
		double t71 = 1.0/t39;
		double t72 = 1.0/t45;
		double t73 = 1.0/t43;
		double t77 = atan(t62);
		double t78 = t52*t77*2.0;
		double t79 = t78-1.0;
		double t80 = t17*t59*5.0E1;
		double t92 = t17*t61*5.0E1;
		double t81 = t53+t55-t57-t58+t80-t92;
		double t82 = t17*t46*t51*2.5E1;
		double t84 = t15*t38*t71*u2*1.0E2;
		double t85 = t15*t41*t72*u3*1.0E2;
		double t86 = t82+t84+t85-t17*t64*t69*2.5E1-t15*t37*t70*u2*1.0E2-t15*t42*t73*u3*1.0E2;
		double t87 = t79*t86*2.251799813685248E15;
		double t88 = t17*t46*t51*1.125899906842624E17;
		double t89 = t62*t62;
		double t90 = t89+1.0;
		double t91 = 1.0/t90;
		double t93 = t17*t46*t51*2.5E7;
		double t94 = t15*t38*t71*u2*1.0E8;
		double t95 = t15*t41*t72*u3*1.0E8;
		double t96 = t93+t94+t95-t17*t64*t69*2.5E7-t15*t37*t70*u2*1.0E8-t15*t42*t73*u3*1.0E8;
		double t97 = t52*t81*t91*t96*4.503599627370496E15;
		double t98 = t15*t38*t71*u2*4.503599627370496E17;
		double t99 = t87+t88+t97+t98-t15*t37*t70*u2*4.503599627370496E17;
		double t100 = t53*4.503599627370496E15;
		double t101 = t55*4.503599627370496E15;
		double t102 = t79*t81*2.251799813685248E15;
		double t103 = t17*t59*2.251799813685248E17;
		double t104 = t100+t101+t102+t103+1.433540275E9;
		double t105 = param1*t7*t8*t15*u2*2.0;
		double t106 = param1*t3*t8*t19*u2*4.0;
		double t107 = param1*t8*t15*u2;
		double t108 = t105-t106;
		double t109 = t47*t108*(1.0/2.0);
		double t110 = t107+t109;
		double t116 = param1*t4*t8*t19*u2*4.0;
		double t111 = t105-t116;
		double t112 = t65*t111*(1.0/2.0);
		double t113 = t107+t112;
		double t114 = t2*t37*t70*1.0E2;
		double t115 = t17*t46*t110*2.5E1;
		double t117 = t114+t115-t2*t38*t71*1.0E2-t17*t64*t113*2.5E1;
		double t118 = t79*t117*2.251799813685248E15;
		double t119 = t2*t37*t70*4.503599627370496E17;
		double t120 = t17*t46*t110*1.125899906842624E17;
		double t121 = t2*t37*t70*1.0E8;
		double t122 = t17*t46*t110*2.5E7;
		double t123 = t121+t122-t2*t38*t71*1.0E8-t17*t64*t113*2.5E7;
		double t124 = t52*t81*t91*t123*4.503599627370496E15;
		double t125 = t118+t119+t120+t124-t2*t38*t71*4.503599627370496E17;
		double t126 = param1*t7*t8*t15*u3*2.0;
		double t127 = param1*t3*t8*t19*u3*4.0;
		double t128 = param1*t8*t15*u3;
		double t129 = t126-t127;
		double t130 = t47*t129*(1.0/2.0);
		double t131 = t128+t130;
		double t136 = param1*t4*t8*t19*u3*4.0;
		double t132 = t126-t136;
		double t133 = t65*t132*(1.0/2.0);
		double t134 = t128+t133;
		double t135 = t2*t41*t72*1.0E2;
		double t137 = t17*t64*t134*2.5E1;
		double t138 = t2*t41*t72*1.0E8;
		double t139 = t17*t64*t134*2.5E7;
		double t140 = t138+t139-t2*t42*t73*1.0E8-t17*t46*t131*2.5E7;
		double t141 = t52*t81*t91*t140*4.503599627370496E15;
		double t142 = param1*t2*t7*t8*2.0;
		double t143 = param1*t2*t8;
		double t147 = param1*t3*t8*t15*4.0;
		double t144 = t142-t147;
		double t145 = t47*t144*(1.0/2.0);
		double t146 = t143+t145;
		double t152 = param1*t4*t8*t15*4.0;
		double t148 = t142-t152;
		double t149 = t65*t148*(1.0/2.0);
		double t150 = t143+t149;
		double t151 = t17*t46*t146*2.5E1;
		double t153 = t151-t17*t64*t150*2.5E1;
		double t154 = t79*t153*2.251799813685248E15;
		double t155 = t17*t46*t146*1.125899906842624E17;
		double t156 = t17*t46*t146*2.5E7;
		double t157 = t156-t17*t64*t150*2.5E7;
		double t158 = t52*t81*t91*t157*4.503599627370496E15;
		double t159 = t154+t155+t158;
		double t160 = t2*u5*2.0;
		double t163 = param1*t2*t8*u5;
		double t161 = t160-t163;
		double t162 = t7*t161*2.0;
		double t164 = param1*t3*t8*t15*u5*4.0;
		double t165 = param1*t8*t14*t19*u5*4.0;
		double t166 = t162+t164+t165;
		double t167 = t47*t166*(1.0/2.0);
		double t168 = t160-t163+t167;
		double t169 = param1*t4*t8*t15*u5*4.0;
		double t170 = t162+t169;
		double t171 = t65*t170*(1.0/2.0);
		double t172 = t160-t163+t171;
		double t173 = t17*t46*t168*2.5E1;
		double t174 = t17*t46*t168*1.125899906842624E17;
		double t175 = t17*t46*t168*2.5E7;
		double t176 = t175-t17*t64*t172*2.5E7;
		double t177 = t52*t81*t91*t176*4.503599627370496E15;
		double t178 = t2*u6*2.0;
		double t181 = param1*t2*t8*u6;
		double t179 = t178-t181;
		double t180 = t7*t179*2.0;
		double t182 = param1*t3*t8*t15*u6*4.0;
		double t183 = t180+t182;
		double t184 = t47*t183*(1.0/2.0);
		double t185 = t178-t181+t184;
		double t186 = param1*t4*t8*t15*u6*4.0;
		double t187 = param1*t8*t14*t19*u6*4.0;
		double t188 = t180+t186+t187;
		double t189 = t65*t188*(1.0/2.0);
		double t190 = t178-t181+t189;
		double t191 = t17*t64*t190*2.5E1;
		double t192 = t17*t46*t185*1.125899906842624E17;
		double t193 = t17*t46*t185*2.5E7;
		double t194 = t193-t17*t64*t190*2.5E7;
		double t195 = t52*t81*t91*t194*4.503599627370496E15;
		fh_udg[0*ng+i] = param2;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = 0.0;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = nl1*t99*t104*uh5*(-9.860761315262648E-36)-nl2*t99*t104*uh6*9.860761315262648E-36;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = param3;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = 0.0;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = 0.0;
		fh_udg[13*ng+i] = nl1*t104*t125*uh5*(-9.860761315262648E-36)-nl2*t104*t125*uh6*9.860761315262648E-36;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = 0.0;
		fh_udg[16*ng+i] = param4;
		fh_udg[17*ng+i] = 0.0;
		fh_udg[18*ng+i] = 0.0;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = nl1*t104*uh5*(t141+t79*(t135+t137+t17*t46*(t47*(t127-param1*t7*t8*t15*u3*2.0)*(1.0/2.0)-param1*t8*t15*u3)*2.5E1-t2*t42*t73*1.0E2)*2.251799813685248E15-t17*t46*t131*1.125899906842624E17)*9.860761315262648E-36+nl2*t104*uh6*(t141+t79*(t135+t137-t2*t42*t73*1.0E2-t17*t46*t131*2.5E1)*2.251799813685248E15-t17*t46*t131*1.125899906842624E17)*9.860761315262648E-36;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = param5;
		fh_udg[25*ng+i] = 0.0;
		fh_udg[26*ng+i] = 0.0;
		fh_udg[27*ng+i] = nl1*t104*t159*uh5*9.860761315262648E-36+nl2*t104*t159*uh6*9.860761315262648E-36;
		fh_udg[28*ng+i] = 0.0;
		fh_udg[29*ng+i] = 0.0;
		fh_udg[30*ng+i] = 0.0;
		fh_udg[31*ng+i] = 0.0;
		fh_udg[32*ng+i] = param6;
		fh_udg[33*ng+i] = 0.0;
		fh_udg[34*ng+i] = nl2*t104*uh6*(t174+t177+t79*(t173-t17*t64*t172*2.5E1)*2.251799813685248E15)*9.860761315262648E-36+nl1*t104*uh5*(t174+t177+t79*(t173-t17*t64*(t160+t171-param1*t2*t8*u5)*2.5E1)*2.251799813685248E15)*9.860761315262648E-36;
		fh_udg[35*ng+i] = 0.0;
		fh_udg[36*ng+i] = 0.0;
		fh_udg[37*ng+i] = 0.0;
		fh_udg[38*ng+i] = 0.0;
		fh_udg[39*ng+i] = 0.0;
		fh_udg[40*ng+i] = param7;
		fh_udg[41*ng+i] = nl2*t104*uh6*(t192+t195-t79*(t191-t17*t46*t185*2.5E1)*2.251799813685248E15)*9.860761315262648E-36+nl1*t104*uh5*(t192+t195-t79*(t191-t17*t46*(t178+t184-param1*t2*t8*u6)*2.5E1)*2.251799813685248E15)*9.860761315262648E-36;
		fh_udg[42*ng+i] = 0.0;
		fh_udg[43*ng+i] = 0.0;
		fh_udg[44*ng+i] = 0.0;
		fh_udg[45*ng+i] = 0.0;
		fh_udg[46*ng+i] = 0.0;
		fh_udg[47*ng+i] = 0.0;
		fh_udg[48*ng+i] = param8;

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
		double t29 = 1.0/u1;
		double t30 = t29*u2*1.0E2;
		double t31 = exp(t30);
		double t32 = t31+1.0;
		double t33 = log(t32);
		double t34 = u5*u5;
		double t35 = u6*u6;
		double t36 = t34+t35;
		double t37 = t29*t36;
		double t39 = param1-1.0;
		double t40 = u1*u4*2.0;
		double t41 = t34*u1;
		double t42 = t35*u1;
		double t43 = u2*u2;
		double t44 = u3*u3;
		double t45 = -t40+t41+t42+t43+t44;
		double t46 = 1.0/(u1*u1);
		double t47 = param1*t39*t45*t46*(1.0/2.0);
		double t38 = t37-t47;
		double t48 = sqrt(2.0);
		double t49 = t38*t38;
		double t50 = 1.0/(u1*u1*u1);
		double t51 = exp(-t30);
		double t52 = t51+1.0;
		double t53 = log(t52);
		double t54 = t29*u3*1.0E2;
		double t55 = exp(t54);
		double t56 = t55+1.0;
		double t57 = log(t56);
		double t58 = param1*t34*t39*t45*t50*2.0;
		double t59 = t49+t58;
		double t60 = sqrt(t59);
		double t61 = t37-t47+t60;
		double t62 = sqrt(t61);
		double t63 = param1*t35*t39*t45*t50*2.0;
		double t64 = t49+t63;
		double t65 = sqrt(t64);
		double t66 = t37-t47+t65;
		double t67 = sqrt(t66);
		double t71 = t33*4.503599627370496E15;
		double t72 = 1.0/3.141592653589793;
		double t73 = t53*1.0E6;
		double t74 = t33*1.0E6;
		double t75 = exp(-t54);
		double t76 = t75+1.0;
		double t77 = log(t76);
		double t78 = t57*1.0E6;
		double t79 = t48*t62*5.0E7;
		double t80 = t48*t67*5.0E7;
		double t81 = t48*t62*5.0E1;
		double t82 = t48*t67*5.0E1;
		double t83 = t33+t53-t57-t77+t81-t82;
		double t84 = t48*t62*2.251799813685248E17;
		double t68 = t71+t84+log(exp(t29*u2*-1.0E2)+1.0)*4.503599627370496E15+t83*(t72*atan(t73+t74-t78+t79-t80-log(exp(t29*u3*-1.0E2)+1.0)*1.0E6)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		double t69 = nl2*t7*uh2;
		double t70 = nl1*t7*uh2;
		double t85 = t53*4.503599627370496E15+t71+t84+t83*(t72*atan(t73+t74-t77*1.0E6-t78+t79-t80)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		fh_uh[0*ng+i] = -param2;
		fh_uh[1*ng+i] = nl1*t7*(t5*-2.0+t15+t16+t17-param1*uh4*2.0)*(-1.0/2.0)+nl1*t2*(t3*-3.0-t4+t10+t11+t12+t13+t14-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)-nl2*t2*t9-nl2*t7*uh5*uh6;
		fh_uh[2*ng+i] = nl2*t7*(t6*-2.0+t15+t16+t17-param1*uh4*2.0)*(-1.0/2.0)+nl2*t2*(-t3-t4*3.0+t10+t11+t12+t13+t14-t6*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)-nl1*t2*t9-nl1*t7*uh5*uh6;
		fh_uh[3*ng+i] = nl1*t18*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t6*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)+nl2*t18*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t6*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t6*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)-nl1*t2*(t5*uh2*-2.0+param1*t5*uh2+param1*t6*uh2-param1*uh2*uh4*2.0+uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t2*(t6*uh3*-2.0+param1*t5*uh3+param1*t6*uh3-param1*uh3*uh4*2.0+uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh_uh[4*ng+i] = nl2*t20;
		fh_uh[5*ng+i] = -nl1*t20;
		fh_uh[6*ng+i] = 0.0;
		fh_uh[7*ng+i] = nl1;
		fh_uh[8*ng+i] = -param3+t27+nl1*t7*(uh2*6.0-param1*uh2*2.0)*(1.0/2.0);
		fh_uh[9*ng+i] = t28+nl2*t7*(uh2*2.0-param1*uh2*2.0)*(1.0/2.0);
		fh_uh[10*ng+i] = nl1*t2*(t3*3.0+t4-t12-t13-t14+t21-param1*t3*3.0+t5*uh1*2.0)*(1.0/2.0)-nl2*t2*t24*(1.0/2.0);
		fh_uh[11*ng+i] = -nl2*t7*uh6;
		fh_uh[12*ng+i] = nl1*t7*uh6;
		fh_uh[13*ng+i] = 0.0;
		fh_uh[14*ng+i] = nl2;
		fh_uh[15*ng+i] = t69+nl1*t7*(uh3*2.0-param1*uh3*2.0)*(1.0/2.0);
		fh_uh[16*ng+i] = -param4+t70+nl2*t7*(uh3*6.0-param1*uh3*2.0)*(1.0/2.0);
		fh_uh[17*ng+i] = nl2*t2*(t3+t4*3.0-t11-t13-t14+t21-param1*t4*3.0+t6*uh1*2.0)*(1.0/2.0)-nl1*t2*t24*(1.0/2.0);
		fh_uh[18*ng+i] = nl2*t7*uh5;
		fh_uh[19*ng+i] = -nl1*t7*uh5;
		fh_uh[20*ng+i] = 0.0;
		fh_uh[21*ng+i] = 0.0;
		fh_uh[22*ng+i] = nl1*t7*t26*(-1.0/2.0);
		fh_uh[23*ng+i] = nl2*t7*t26*(-1.0/2.0);
		fh_uh[24*ng+i] = -param5+nl1*param1*t7*uh2+nl2*param1*t7*uh3;
		fh_uh[25*ng+i] = 0.0;
		fh_uh[26*ng+i] = 0.0;
		fh_uh[27*ng+i] = 0.0;
		fh_uh[28*ng+i] = 0.0;
		fh_uh[29*ng+i] = -nl2*uh6-nl1*param1*uh5;
		fh_uh[30*ng+i] = -nl1*uh6+nl2*t7*(uh1*uh5*4.0-param1*uh1*uh5*2.0)*(1.0/2.0);
		fh_uh[31*ng+i] = nl2*t2*(uh1*uh2*uh6*2.0-uh1*uh3*uh5*4.0+param1*uh1*uh3*uh5*2.0)*(-1.0/2.0)-nl1*t2*(uh1*uh3*uh6*2.0+param1*uh1*uh2*uh5*2.0)*(1.0/2.0);
		fh_uh[32*ng+i] = -param6+t27;
		fh_uh[33*ng+i] = -t28;
		fh_uh[34*ng+i] = nl1*(t68*t68)*4.930380657631324E-36;
		fh_uh[35*ng+i] = 0.0;
		fh_uh[36*ng+i] = -nl2*uh5+nl1*t7*(uh1*uh6*4.0-param1*uh1*uh6*2.0)*(1.0/2.0);
		fh_uh[37*ng+i] = -nl1*uh5-nl2*param1*uh6;
		fh_uh[38*ng+i] = nl1*t2*(uh1*uh2*uh6*-4.0+uh1*uh3*uh5*2.0+param1*uh1*uh2*uh6*2.0)*(-1.0/2.0)-nl2*t2*(uh1*uh2*uh5*2.0+param1*uh1*uh3*uh6*2.0)*(1.0/2.0);
		fh_uh[39*ng+i] = -t69;
		fh_uh[40*ng+i] = -param7+t70;
		fh_uh[41*ng+i] = nl2*(t85*t85)*4.930380657631324E-36;
		fh_uh[42*ng+i] = 0.0;
		fh_uh[43*ng+i] = 0.0;
		fh_uh[44*ng+i] = 0.0;
		fh_uh[45*ng+i] = 0.0;
		fh_uh[46*ng+i] = nl1;
		fh_uh[47*ng+i] = nl2;
		fh_uh[48*ng+i] = -param8;

	}
}

void fhatonly_mhd2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];
	double param7 = param[6];
	double param8 = param[7];

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
		double t17 = 1.0/u1;
		double t18 = t17*u2*1.0E2;
		double t19 = exp(t18);
		double t20 = t19+1.0;
		double t21 = log(t20);
		double t22 = u5*u5;
		double t23 = u6*u6;
		double t24 = t22+t23;
		double t25 = t17*t24;
		double t27 = param1-1.0;
		double t28 = u1*u4*2.0;
		double t29 = t22*u1;
		double t30 = t23*u1;
		double t31 = u2*u2;
		double t32 = u3*u3;
		double t33 = -t28+t29+t30+t31+t32;
		double t34 = 1.0/(u1*u1);
		double t35 = param1*t27*t33*t34*(1.0/2.0);
		double t26 = t25-t35;
		double t36 = sqrt(2.0);
		double t37 = t26*t26;
		double t38 = 1.0/(u1*u1*u1);
		double t39 = exp(-t18);
		double t40 = t39+1.0;
		double t41 = log(t40);
		double t42 = t17*u3*1.0E2;
		double t43 = exp(t42);
		double t44 = t43+1.0;
		double t45 = log(t44);
		double t46 = param1*t22*t27*t33*t38*2.0;
		double t47 = t37+t46;
		double t48 = sqrt(t47);
		double t49 = t25-t35+t48;
		double t50 = sqrt(t49);
		double t51 = param1*t23*t27*t33*t38*2.0;
		double t52 = t37+t51;
		double t53 = sqrt(t52);
		double t54 = t25-t35+t53;
		double t55 = sqrt(t54);
		double t57 = t21*4.503599627370496E15;
		double t58 = 1.0/3.141592653589793;
		double t59 = t41*1.0E6;
		double t60 = t21*1.0E6;
		double t61 = exp(-t42);
		double t62 = t61+1.0;
		double t63 = log(t62);
		double t64 = t45*1.0E6;
		double t65 = t36*t50*5.0E7;
		double t66 = t36*t55*5.0E7;
		double t67 = t36*t50*5.0E1;
		double t68 = t36*t55*5.0E1;
		double t69 = t21+t41-t45-t63+t67-t68;
		double t70 = t36*t50*2.251799813685248E17;
		double t56 = t57+t70+log(exp(t17*u2*-1.0E2)+1.0)*4.503599627370496E15+t69*(t58*atan(t59+t60-t64+t65-t66-log(exp(t17*u3*-1.0E2)+1.0)*1.0E6)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		double t71 = t41*4.503599627370496E15+t57+t70+t69*(t58*atan(t59+t60-t63*1.0E6-t64+t65-t66)*2.0-1.0)*2.251799813685248E15+1.433540275E9;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+param2*(u1-uh1);
		fh[1*ng+i] = param3*(u2-uh2)-nl1*t2*(t3*-3.0-t4+t8+t9+t10+t12+t13-t5*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl2*t2*t7;
		fh[2*ng+i] = param4*(u3-uh3)-nl2*t2*(-t3-t4*3.0+t8+t9+t10+t12+t13-t11*uh1*2.0-param1*uh1*uh4*2.0)*(1.0/2.0)+nl1*t2*t7;
		fh[3*ng+i] = param5*(u4-uh4)-nl1*t14*(-t3*uh2-t4*uh2+param1*t3*uh2+param1*t4*uh2-t5*uh1*uh2*2.0+param1*t5*uh1*uh2+param1*t11*uh1*uh2-param1*uh1*uh2*uh4*2.0+uh1*uh3*uh5*uh6*2.0)*(1.0/2.0)-nl2*t14*(-t3*uh3-t4*uh3+param1*t3*uh3+param1*t4*uh3-t11*uh1*uh3*2.0+param1*t5*uh1*uh3+param1*t11*uh1*uh3-param1*uh1*uh3*uh4*2.0+uh1*uh2*uh5*uh6*2.0)*(1.0/2.0);
		fh[4*ng+i] = -nl2*t16+nl1*uh7+param6*(u5-uh5);
		fh[5*ng+i] = nl1*t16+nl2*uh7+param7*(u6-uh6);
		fh[6*ng+i] = param8*(u7-uh7)+nl1*(t56*t56)*uh5*4.930380657631324E-36+nl2*(t71*t71)*uh6*4.930380657631324E-36;

	}
}

