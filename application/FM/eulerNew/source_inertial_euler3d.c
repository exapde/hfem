
// Written by C. Nguyen and P. Fernandez

void source_inertial_euler3d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
// Source term for inertial forces in a rotating frame when solving in relative velocity.
		
	double RO_x = 0.0;
	double RO_y = 0.0;
	double RO_z = 0.0;

	double W_x = 0.0;
	double W_y = 0.0;
	double W_z = 0.0;

	double dWdt_x = 0.0;
	double dWdt_y = 0.0;
	double dWdt_z = 0.0;

	error("Need to determine angular velocity, angular acceleration and location of center of rotation.\n");

	for (int i = 0; i <ng; i++) {
		double OX = pg[0*ng+i];
		double OY = pg[1*ng+i];
		double OZ = pg[2*ng+i];

        double r = udg[0*ng+i];
        double ru = udg[1*ng+i];
        double rv = udg[2*ng+i];
        double rw = udg[3*ng+i];

		double t2 = OX+RO_x;
		double t3 = OY+RO_y;
		double t4 = OZ+RO_z;
		double t5 = W_x*t3;
		double t6 = t5-W_y*t2;
		double t7 = W_x*t4;
		double t8 = t7-W_z*t2;
		double t9 = W_y*t4;
		double t10 = t9-W_z*t3;

		s[1*ng+i] += W_z*rv*2.0-W_y*rw*2.0-r*(W_y*t6+W_z*t8)-r*(dWdt_y*t4-dWdt_z*t3);
		s[2*ng+i] += W_z*ru*-2.0+W_x*rw*2.0+r*(W_x*t6-W_z*t10)+r*(dWdt_x*t4-dWdt_z*t2);
		s[3*ng+i] += W_y*ru*2.0-W_x*rv*2.0+r*(W_x*t8+W_y*t10)-r*(dWdt_x*t3-dWdt_y*t2);
	}

    if (computeJacobian == 1) {
        for (int i = 0; i <ng; i++) {
            double OX = pg[0*ng+i];
            double OY = pg[1*ng+i];
            double OZ = pg[2*ng+i];

            double t2 = OX+RO_x;
            double t3 = OY+RO_y;
            double t4 = OZ+RO_z;
            double t5 = W_x*t3;
            double t6 = W_x*t4;
            double t7 = W_y*t4;
            double t8 = W_y*2.0;

            s_udg[1*ng+i] += -dWdt_y*t4+dWdt_z*t3-W_y*(t5-W_y*t2)-W_z*(t6-W_z*t2);
            s_udg[2*ng+i] += dWdt_x*t4-dWdt_z*t2+W_x*(t5-W_y*t2)-W_z*(t7-W_z*t3);
            s_udg[3*ng+i] += -dWdt_x*t3+dWdt_y*t2+W_x*(t6-W_z*t2)+W_y*(t7-W_z*t3);
            s_udg[7*ng+i] += W_z*-2.0;
            s_udg[8*ng+i] += t8;
            s_udg[11*ng+i] += W_z*2.0;
            s_udg[13*ng+i] += W_x*-2.0;
            s_udg[16*ng+i] += -t8;
            s_udg[17*ng+i] += W_x*2.0;
        }
    }
}
