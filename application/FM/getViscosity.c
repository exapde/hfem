
#define PI (3.141592653589793)

// Written by: C. Nguyen & P. Fernandez

void constantDynamicViscosity(double *mu, double *mu_udg, int ng, int nc, int computeJacobian)
{
    // mu_udg is not computed in this case since its value is always zero
    int i;

    for (i=0; i<ng; i++) {
        mu[i] = 1.0;
    }

//    if (computeJacobian == 1) {
//        int j;
//        for (j=0; j<nc; j++)
//            for (i=0; i<ng; i++)
//                mu_udg[j*ng+i] = 0.0;
//    }
}



void sutherlandLaw2d(double *mu, double *mu_udg, double *pg, double *udg, double *param, int ng, int nc, int ncu, int computeJacobian)
{
    // This function assumes that rho_inf = || v_inf || = 1 with the employed non-dimensionalization.

    double param1 = param[0];
    double param2 = param[1];
    double param3 = param[2];
    double param4 = param[3];
    double param5 = param[4];
    double param6 = param[5];

    // Sutherland's law parameters:
    double Ssutherland = 111;       // Sutherland's constant
    double Tinf = 300;              // Dimensional temperature (in the same units as Ssutherland)

    // Regularization constants:
    double b = 100.0;
    double c = 1.0/2.0 - atan(b)/PI;

    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];

        double u1 = udg[0*ng+i];
        double u2 = udg[1*ng+i];
        double u3 = udg[2*ng+i];
        double u4 = udg[3*ng+i];

        double t2 = 1.0/(u1*u1);
        double t3 = Ssutherland+Tinf;
        double t4 = 1.0/param1;
        double t5 = 1.0/(param5*param5);
        double t6 = 1.0/u1;
        double t7 = u2*u2;
        double t8 = t2*t7*(1.0/2.0);
        double t9 = u3*u3;
        double t10 = t2*t9*(1.0/2.0);
        double t11 = t8+t10;
        double t18 = t11*u1;
        double t12 = -t18+u4;
        double t13 = param1-1.0;
        double t14 = 1.0/t13;
        double t15 = Tinf*t4*t5*t6*t12*t14;
        double t16 = Ssutherland+t15;
        double t17 = 1.0/t16;
        double t19 = t4*t5*t6*t12*t14;
        double t20 = pow(t19,3.0/2.0);

        mu[i] = c+t3*t17*t20*(atan(b*t3*t17*t20)/3.141592653589793+1.0/2.0);
    }

    if (computeJacobian == 1) {

        for (int i = 0; i <ng; i++) {
            double x1 = pg[0*ng+i];
            double x2 = pg[1*ng+i];

            double u1 = udg[0*ng+i];
            double u2 = udg[1*ng+i];
            double u3 = udg[2*ng+i];
            double u4 = udg[3*ng+i];

            double t2 = 1.0/(u1*u1);
            double t3 = 1.0/param1;
            double t4 = 1.0/(param5*param5);
            double t5 = u2*u2;
            double t6 = t2*t5*(1.0/2.0);
            double t7 = u3*u3;
            double t8 = t2*t7*(1.0/2.0);
            double t9 = t6+t8;
            double t16 = t9*u1;
            double t10 = -t16+u4;
            double t11 = param1-1.0;
            double t12 = 1.0/t11;
            double t13 = 1.0/u1;
            double t14 = 1.0/(u1*u1*u1);
            double t15 = Ssutherland+Tinf;
            double t17 = Tinf*t3*t4*t10*t12*t13;
            double t18 = Ssutherland+t17;
            double t19 = 1.0/t18;
            double t20 = t3*t4*t10*t12*t13;
            double t21 = 1.0/3.141592653589793;
            double t22 = pow(t20,3.0/2.0);
            double t23 = b*t15*t19*t22;
            double t24 = atan(t23);
            double t25 = t21*t24;
            double t26 = t25+1.0/2.0;
            double t27 = t5*t14;
            double t28 = t7*t14;
            double t29 = t27+t28;
            double t33 = t29*u1;
            double t30 = t6+t8-t33;
            double t31 = t10*t10;
            double t32 = 1.0/(t18*t18);
            double t34 = Tinf*t3*t4*t12*t13*t30;
            double t35 = Tinf*t2*t3*t4*t10*t12;
            double t36 = t34+t35;
            double t37 = t2*t3*t4*t10*t12;
            double t38 = t3*t4*t12*t13*t30;
            double t39 = t37+t38;
            double t40 = sqrt(t20);
            double t41 = b*b;
            double t42 = 1.0/(param1*param1*param1);
            double t43 = 1.0/(param5*param5*param5*param5*param5*param5);
            double t44 = t15*t15;
            double t45 = 1.0/(t11*t11*t11);
            double t46 = t10*t14*t31*t32*t41*t42*t43*t44*t45;
            double t47 = t46+1.0;
            double t48 = 1.0/t47;

            mu_udg[0*ng+i] = t15*t22*t26*t32*t36-t15*t19*t26*t39*t40*(3.0/2.0)+t15*t19*t21*t22*t48*(b*t15*t22*t32*t36-b*t15*t19*t39*t40*(3.0/2.0));
            mu_udg[1*ng+i] = -t15*t19*t21*t22*t48*(b*t2*t3*t4*t12*t15*t19*t40*u2*(3.0/2.0)-Tinf*b*t2*t3*t4*t12*t15*t22*t32*u2)-t2*t3*t4*t12*t15*t19*t26*t40*u2*(3.0/2.0)+Tinf*t2*t3*t4*t12*t15*t22*t26*t32*u2;
            mu_udg[2*ng+i] = -t15*t19*t21*t22*t48*(b*t2*t3*t4*t12*t15*t19*t40*u3*(3.0/2.0)-Tinf*b*t2*t3*t4*t12*t15*t22*t32*u3)-t2*t3*t4*t12*t15*t19*t26*t40*u3*(3.0/2.0)+Tinf*t2*t3*t4*t12*t15*t22*t26*t32*u3;
            mu_udg[3*ng+i] = t15*t19*t21*t22*t48*(b*t3*t4*t12*t13*t15*t19*t40*(3.0/2.0)-Tinf*b*t3*t4*t12*t13*t15*t22*t32)+t3*t4*t12*t13*t15*t19*t26*t40*(3.0/2.0)-Tinf*t3*t4*t12*t13*t15*t22*t26*t32;
            for (int j=4; j<nc; j++) {
                mu_udg[j*ng+i] = 0.0;
            }
        }
    }
}



//void sutherlandLaw2d(double *mu, double *mu_udg, double *pg, double *udg, double *param, int ng, int nc, int ncu, int computeJacobian)
//{
//    // This function assumes that rho_inf = || v_inf || = 1 with the employed non-dimensionalization.
//
//    double param1 = param[0];
//    double param2 = param[1];
//    double param3 = param[2];
//    double param4 = param[3];
//    double param5 = param[4];
//    double param6 = param[5];
//
//    // Sutherland's law parameters:
//    double Ssutherland = 111;       // Sutherland's constant
//    double Tinf = 300;              // Dimensional temperature (in the same units as Ssutherland)
//
//    // Regularization constants:
//    double b = 100.0;
//    double c = 1.0/2.0 - atan(b)/PI;
//
//    for (int i = 0; i <ng; i++) {
//        double x1 = pg[0*ng+i];
//        double x2 = pg[1*ng+i];
//
//        double u1 = udg[0*ng+i];
//        double u2 = udg[1*ng+i];
//        double u3 = udg[2*ng+i];
//        double u4 = udg[3*ng+i];
//
//        double t2 = 1.0/u1;
//        double t3 = 1.0/(u1*u1);
//        double t4 = Ssutherland+Tinf;
//        double t5 = 1.0/param1;
//        double t6 = 1.0/(param5*param5);
//        double t7 = u2*u2;
//        double t8 = t3*t7*(1.0/2.0);
//        double t9 = u3*u3;
//        double t10 = t3*t9*(1.0/2.0);
//        double t11 = t8+t10;
//        double t18 = t11*u1;
//        double t12 = -t18+u4;
//        double t13 = param1-1.0;
//        double t14 = 1.0/t13;
//        double t15 = Tinf*t2*t5*t6*t12*t14;
//        double t16 = Ssutherland+t15;
//        double t17 = 1.0/t16;
//        double t19 = t2*t5*t6*t12*t14;
//        double t20 = pow(t19,3.0/2.0);
//
//        mu[i] = c+t2*t4*t17*t20*(atan(b*t2*t4*t17*t20)/3.141592653589793+1.0/2.0);
//    }
//
//    if (computeJacobian == 1) {
//
//        for (int i = 0; i <ng; i++) {
//            double x1 = pg[0*ng+i];
//            double x2 = pg[1*ng+i];
//
//            double u1 = udg[0*ng+i];
//            double u2 = udg[1*ng+i];
//            double u3 = udg[2*ng+i];
//            double u4 = udg[3*ng+i];
//
//            double t2 = 1.0/(u1*u1);
//            double t3 = 1.0/u1;
//            double t4 = Ssutherland+Tinf;
//            double t5 = 1.0/param1;
//            double t6 = 1.0/(param5*param5);
//            double t7 = u2*u2;
//            double t8 = t2*t7*(1.0/2.0);
//            double t9 = u3*u3;
//            double t10 = t2*t9*(1.0/2.0);
//            double t11 = t8+t10;
//            double t18 = t11*u1;
//            double t12 = -t18+u4;
//            double t13 = param1-1.0;
//            double t14 = 1.0/t13;
//            double t15 = Tinf*t3*t5*t6*t12*t14;
//            double t16 = Ssutherland+t15;
//            double t17 = 1.0/t16;
//            double t19 = t3*t5*t6*t12*t14;
//            double t20 = pow(t19,3.0/2.0);
//            double t21 = 1.0/3.141592653589793;
//            double t22 = b*t3*t4*t17*t20;
//            double t23 = atan(t22);
//            double t24 = t21*t23;
//            double t25 = t24+1.0/2.0;
//            double t26 = 1.0/(u1*u1*u1);
//            double t27 = t7*t26;
//            double t28 = t9*t26;
//            double t29 = t27+t28;
//            double t33 = t29*u1;
//            double t30 = t8+t10-t33;
//            double t31 = t12*t12;
//            double t32 = 1.0/(t16*t16);
//            double t34 = Tinf*t3*t5*t6*t14*t30;
//            double t35 = Tinf*t2*t5*t6*t12*t14;
//            double t36 = t34+t35;
//            double t37 = t2*t5*t6*t12*t14;
//            double t38 = t3*t5*t6*t14*t30;
//            double t39 = t37+t38;
//            double t40 = sqrt(t19);
//            double t41 = b*b;
//            double t42 = 1.0/(param1*param1*param1);
//            double t43 = 1.0/(param5*param5*param5*param5*param5*param5);
//            double t44 = 1.0/(u1*u1*u1*u1*u1);
//            double t45 = t4*t4;
//            double t46 = 1.0/(t13*t13*t13);
//            double t47 = t12*t31*t32*t41*t42*t43*t44*t45*t46;
//            double t48 = t47+1.0;
//            double t49 = 1.0/t48;
//
//            mu_udg[0*ng+i] = -t2*t4*t17*t20*t25+t3*t4*t20*t25*t32*t36-t3*t4*t17*t25*t39*t40*(3.0/2.0)-t3*t4*t17*t20*t21*t49*(b*t2*t4*t17*t20-b*t3*t4*t20*t32*t36+b*t3*t4*t17*t39*t40*(3.0/2.0));
//            mu_udg[1*ng+i] = -t3*t4*t17*t20*t21*t49*(b*t4*t5*t6*t14*t17*t26*t40*u2*(3.0/2.0)-Tinf*b*t4*t5*t6*t14*t20*t26*t32*u2)-t4*t5*t6*t14*t17*t25*t26*t40*u2*(3.0/2.0)+Tinf*t4*t5*t6*t14*t20*t25*t26*t32*u2;
//            mu_udg[2*ng+i] = -t3*t4*t17*t20*t21*t49*(b*t4*t5*t6*t14*t17*t26*t40*u3*(3.0/2.0)-Tinf*b*t4*t5*t6*t14*t20*t26*t32*u3)-t4*t5*t6*t14*t17*t25*t26*t40*u3*(3.0/2.0)+Tinf*t4*t5*t6*t14*t20*t25*t26*t32*u3;
//            mu_udg[3*ng+i] = t3*t4*t17*t20*t21*t49*(b*t2*t4*t5*t6*t14*t17*t40*(3.0/2.0)-Tinf*b*t2*t4*t5*t6*t14*t20*t32)+t2*t4*t5*t6*t14*t17*t25*t40*(3.0/2.0)-Tinf*t2*t4*t5*t6*t14*t20*t25*t32;
//
//            for (int j=4; j<nc; j++) {
//                mu_udg[j*ng+i] = 0.0;
//            }
//        }
//    }
//}



void sutherlandLaw3d(double *mu, double *mu_udg, double *pg, double *udg, double *param, int ng, int nc, int ncu, int computeJacobian)
{
    // This function assumes that rho_inf = || v_inf || = 1 with the employed non-dimensionalization.

    double param1 = param[0];
    double param2 = param[1];
    double param3 = param[2];
    double param4 = param[3];
    double param5 = param[4];
    double param6 = param[5];

    // Sutherland's law parameters:
    double Ssutherland = 111;       // Sutherland's constant
    double Tinf = 300;              // Dimensional temperature (in the same units as Ssutherland)

    // Regularization constants:
    double b = 100.0;
    double c = 1.0/2.0 - atan(b)/PI;

    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];
        double x3 = pg[2*ng+i];

        double u1 = udg[0*ng+i];
        double u2 = udg[1*ng+i];
        double u3 = udg[2*ng+i];
        double u4 = udg[3*ng+i];
        double u5 = udg[4*ng+i];

        double t2 = 1.0/(u1*u1);
        double t3 = Ssutherland+Tinf;
        double t4 = 1.0/param1;
        double t5 = 1.0/(param5*param5);
        double t6 = 1.0/u1;
        double t7 = u2*u2;
        double t8 = t2*t7*(1.0/2.0);
        double t9 = u3*u3;
        double t10 = t2*t9*(1.0/2.0);
        double t11 = u4*u4;
        double t12 = t2*t11*(1.0/2.0);
        double t13 = t8+t10+t12;
        double t20 = t13*u1;
        double t14 = -t20+u5;
        double t15 = param1-1.0;
        double t16 = 1.0/t15;
        double t17 = Tinf*t4*t5*t6*t14*t16;
        double t18 = Ssutherland+t17;
        double t19 = 1.0/t18;
        double t21 = t4*t5*t6*t14*t16;
        double t22 = pow(t21,3.0/2.0);

        mu[i] = c+t3*t19*t22*(atan(b*t3*t19*t22)/3.141592653589793+1.0/2.0);
    }

    if (computeJacobian == 1) {

        for (int i = 0; i <ng; i++) {
            double x1 = pg[0*ng+i];
            double x2 = pg[1*ng+i];
            double x3 = pg[2*ng+i];

            double u1 = udg[0*ng+i];
            double u2 = udg[1*ng+i];
            double u3 = udg[2*ng+i];
            double u4 = udg[3*ng+i];
            double u5 = udg[4*ng+i];

            double t2 = 1.0/(u1*u1);
            double t3 = Ssutherland+Tinf;
            double t4 = 1.0/param1;
            double t5 = 1.0/(param5*param5);
            double t6 = 1.0/u1;
            double t7 = u2*u2;
            double t8 = t2*t7*(1.0/2.0);
            double t9 = u3*u3;
            double t10 = t2*t9*(1.0/2.0);
            double t11 = u4*u4;
            double t12 = t2*t11*(1.0/2.0);
            double t13 = t8+t10+t12;
            double t19 = t13*u1;
            double t14 = -t19+u5;
            double t15 = param1-1.0;
            double t16 = 1.0/t15;
            double t17 = Tinf*t4*t5*t6*t14*t16;
            double t18 = Ssutherland+t17;
            double t20 = 1.0/(u1*u1*u1);
            double t21 = t4*t5*t6*t14*t16;
            double t22 = pow(t21,3.0/2.0);
            double t23 = 1.0/t18;
            double t24 = t7*t20;
            double t25 = t9*t20;
            double t26 = t11*t20;
            double t27 = t24+t25+t26;
            double t35 = t27*u1;
            double t28 = t8+t10+t12-t35;
            double t29 = 1.0/3.141592653589793;
            double t30 = b*t3*t22*t23;
            double t31 = atan(t30);
            double t32 = t29*t31;
            double t33 = t32+1.0/2.0;
            double t34 = t2*t4*t5*t14*t16;
            double t36 = t4*t5*t6*t16*t28;
            double t37 = t34+t36;
            double t38 = sqrt(t21);
            double t39 = 1.0/(t18*t18);
            double t40 = Tinf*t4*t5*t6*t16*t28;
            double t41 = Tinf*t2*t4*t5*t14*t16;
            double t42 = t40+t41;
            double t43 = t14*t14;
            double t44 = b*b;
            double t45 = 1.0/(param1*param1*param1);
            double t46 = 1.0/(param5*param5*param5*param5*param5*param5);
            double t47 = t3*t3;
            double t48 = 1.0/(t15*t15*t15);
            double t49 = t14*t20*t39*t43*t44*t45*t46*t47*t48;
            double t50 = t49+1.0;
            double t51 = 1.0/t50;

            mu_udg[0*ng+i] = t3*t23*t33*t37*t38*(-3.0/2.0)+t3*t22*t33*t39*t42-t3*t22*t23*t29*t51*(b*t3*t23*t37*t38*(3.0/2.0)-b*t3*t22*t39*t42);
            mu_udg[1*ng+i] = -t3*t22*t23*t29*t51*(b*t2*t3*t4*t5*t16*t23*t38*u2*(3.0/2.0)-Tinf*b*t2*t3*t4*t5*t16*t22*t39*u2)-t2*t3*t4*t5*t16*t23*t33*t38*u2*(3.0/2.0)+Tinf*t2*t3*t4*t5*t16*t22*t33*t39*u2;
            mu_udg[2*ng+i] = -t3*t22*t23*t29*t51*(b*t2*t3*t4*t5*t16*t23*t38*u3*(3.0/2.0)-Tinf*b*t2*t3*t4*t5*t16*t22*t39*u3)-t2*t3*t4*t5*t16*t23*t33*t38*u3*(3.0/2.0)+Tinf*t2*t3*t4*t5*t16*t22*t33*t39*u3;
            mu_udg[3*ng+i] = -t3*t22*t23*t29*t51*(b*t2*t3*t4*t5*t16*t23*t38*u4*(3.0/2.0)-Tinf*b*t2*t3*t4*t5*t16*t22*t39*u4)-t2*t3*t4*t5*t16*t23*t33*t38*u4*(3.0/2.0)+Tinf*t2*t3*t4*t5*t16*t22*t33*t39*u4;
            mu_udg[4*ng+i] = t3*t22*t23*t29*t51*(b*t3*t4*t5*t6*t16*t23*t38*(3.0/2.0)-Tinf*b*t3*t4*t5*t6*t16*t22*t39)+t3*t4*t5*t6*t16*t23*t33*t38*(3.0/2.0)-Tinf*t3*t4*t5*t6*t16*t22*t33*t39;

            for (int j=5; j<nc; j++) {
                mu_udg[j*ng+i] = 0.0;
            }
        }
    }
}



//void sutherlandLaw3d(double *mu, double *mu_udg, double *pg, double *udg, double *param, int ng, int nc, int ncu, int computeJacobian)
//{
//    // This function assumes that rho_inf = || v_inf || = 1 with the employed non-dimensionalization.
//
//    double param1 = param[0];
//    double param2 = param[1];
//    double param3 = param[2];
//    double param4 = param[3];
//    double param5 = param[4];
//    double param6 = param[5];
//
//    // Sutherland's law parameters:
//    double Ssutherland = 111;       // Sutherland's constant
//    double Tinf = 300;              // Dimensional temperature (in the same units as Ssutherland)
//
//    // Regularization constants:
//    double b = 100.0;
//    double c = 1.0/2.0 - atan(b)/PI;
//
//    for (int i = 0; i <ng; i++) {
//        double x1 = pg[0*ng+i];
//        double x2 = pg[1*ng+i];
//        double x3 = pg[2*ng+i];
//
//        double u1 = udg[0*ng+i];
//        double u2 = udg[1*ng+i];
//        double u3 = udg[2*ng+i];
//        double u4 = udg[3*ng+i];
//        double u5 = udg[4*ng+i];
//
//        double t2 = 1.0/u1;
//        double t3 = 1.0/(u1*u1);
//        double t4 = Ssutherland+Tinf;
//        double t5 = 1.0/param1;
//        double t6 = 1.0/(param5*param5);
//        double t7 = u2*u2;
//        double t8 = t3*t7*(1.0/2.0);
//        double t9 = u3*u3;
//        double t10 = t3*t9*(1.0/2.0);
//        double t11 = u4*u4;
//        double t12 = t3*t11*(1.0/2.0);
//        double t13 = t8+t10+t12;
//        double t20 = t13*u1;
//        double t14 = -t20+u5;
//        double t15 = param1-1.0;
//        double t16 = 1.0/t15;
//        double t17 = Tinf*t2*t5*t6*t14*t16;
//        double t18 = Ssutherland+t17;
//        double t19 = 1.0/t18;
//        double t21 = t2*t5*t6*t14*t16;
//        double t22 = pow(t21,3.0/2.0);
//
//        mu[i] = c+t2*t4*t19*t22*(atan(b*t2*t4*t19*t22)/3.141592653589793+1.0/2.0);
//    }
//
//    if (computeJacobian == 1) {
//
//        for (int i = 0; i <ng; i++) {
//            double x1 = pg[0*ng+i];
//            double x2 = pg[1*ng+i];
//            double x3 = pg[2*ng+i];
//
//            double u1 = udg[0*ng+i];
//            double u2 = udg[1*ng+i];
//            double u3 = udg[2*ng+i];
//            double u4 = udg[3*ng+i];
//            double u5 = udg[4*ng+i];
//
//            double t2 = 1.0/(u1*u1);
//            double t3 = 1.0/u1;
//            double t4 = Ssutherland+Tinf;
//            double t5 = 1.0/param1;
//            double t6 = 1.0/(param5*param5);
//            double t7 = u2*u2;
//            double t8 = t2*t7*(1.0/2.0);
//            double t9 = u3*u3;
//            double t10 = t2*t9*(1.0/2.0);
//            double t11 = u4*u4;
//            double t12 = t2*t11*(1.0/2.0);
//            double t13 = t8+t10+t12;
//            double t20 = t13*u1;
//            double t14 = -t20+u5;
//            double t15 = param1-1.0;
//            double t16 = 1.0/t15;
//            double t17 = Tinf*t3*t5*t6*t14*t16;
//            double t18 = Ssutherland+t17;
//            double t19 = 1.0/t18;
//            double t21 = t3*t5*t6*t14*t16;
//            double t22 = pow(t21,3.0/2.0);
//            double t23 = 1.0/(u1*u1*u1);
//            double t24 = 1.0/3.141592653589793;
//            double t25 = b*t3*t4*t19*t22;
//            double t26 = atan(t25);
//            double t27 = t24*t26;
//            double t28 = t27+1.0/2.0;
//            double t29 = t7*t23;
//            double t30 = t9*t23;
//            double t31 = t11*t23;
//            double t32 = t29+t30+t31;
//            double t36 = t32*u1;
//            double t33 = t8+t10+t12-t36;
//            double t34 = t14*t14;
//            double t35 = 1.0/(t18*t18);
//            double t37 = Tinf*t3*t5*t6*t16*t33;
//            double t38 = Tinf*t2*t5*t6*t14*t16;
//            double t39 = t37+t38;
//            double t40 = t2*t5*t6*t14*t16;
//            double t41 = t3*t5*t6*t16*t33;
//            double t42 = t40+t41;
//            double t43 = sqrt(t21);
//            double t44 = b*b;
//            double t45 = 1.0/(param1*param1*param1);
//            double t46 = 1.0/(param5*param5*param5*param5*param5*param5);
//            double t47 = 1.0/(u1*u1*u1*u1*u1);
//            double t48 = t4*t4;
//            double t49 = 1.0/(t15*t15*t15);
//            double t50 = t14*t34*t35*t44*t45*t46*t47*t48*t49;
//            double t51 = t50+1.0;
//            double t52 = 1.0/t51;
//
//            mu_udg[0*ng+i] = -t2*t4*t19*t22*t28+t3*t4*t22*t28*t35*t39-t3*t4*t19*t28*t42*t43*(3.0/2.0)-t3*t4*t19*t22*t24*t52*(b*t2*t4*t19*t22-b*t3*t4*t22*t35*t39+b*t3*t4*t19*t42*t43*(3.0/2.0));
//            mu_udg[1*ng+i] = -t3*t4*t19*t22*t24*t52*(b*t4*t5*t6*t16*t19*t23*t43*u2*(3.0/2.0)-Tinf*b*t4*t5*t6*t16*t22*t23*t35*u2)-t4*t5*t6*t16*t19*t23*t28*t43*u2*(3.0/2.0)+Tinf*t4*t5*t6*t16*t22*t23*t28*t35*u2;
//            mu_udg[2*ng+i] = -t3*t4*t19*t22*t24*t52*(b*t4*t5*t6*t16*t19*t23*t43*u3*(3.0/2.0)-Tinf*b*t4*t5*t6*t16*t22*t23*t35*u3)-t4*t5*t6*t16*t19*t23*t28*t43*u3*(3.0/2.0)+Tinf*t4*t5*t6*t16*t22*t23*t28*t35*u3;
//            mu_udg[3*ng+i] = -t3*t4*t19*t22*t24*t52*(b*t4*t5*t6*t16*t19*t23*t43*u4*(3.0/2.0)-Tinf*b*t4*t5*t6*t16*t22*t23*t35*u4)-t4*t5*t6*t16*t19*t23*t28*t43*u4*(3.0/2.0)+Tinf*t4*t5*t6*t16*t22*t23*t28*t35*u4;
//            mu_udg[4*ng+i] = t3*t4*t19*t22*t24*t52*(b*t2*t4*t5*t6*t16*t19*t43*(3.0/2.0)-Tinf*b*t2*t4*t5*t6*t16*t22*t35)+t2*t4*t5*t6*t16*t19*t28*t43*(3.0/2.0)-Tinf*t2*t4*t5*t6*t16*t22*t28*t35;
//
//            for (int j=5; j<nc; j++) {
//                mu_udg[j*ng+i] = 0.0;
//            }
//        }
//    }
//}



void getViscosity(double *mu, double *mu_udg, double *pg, double *udg, double *param, Int viscosityModel, int ng, int nc, int ncu, int nd, int computeJacobian)
{

    if (viscosityModel == 0)            // Constant dynamic viscosity
        constantDynamicViscosity(mu, mu_udg, ng, nc, computeJacobian);
    else if (viscosityModel == 1) {     // Sutherland's law
        if (nd == 2)
            sutherlandLaw2d(mu, mu_udg, pg, udg, param, ng, nc, ncu, computeJacobian);
        else if (nd == 3)
            sutherlandLaw3d(mu, mu_udg, pg, udg, param, ng, nc, ncu, computeJacobian);
    }
    else {
        printf("Viscosity model viscosityModel = %d not implemented\n", viscosityModel);
        printf("Execution will be terminated\n");
        exit(-1);
    }

}
