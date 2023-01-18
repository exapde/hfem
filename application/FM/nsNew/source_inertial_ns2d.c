
// Written by C. Nguyen and P. Fernandez

void source_inertial_ns2d(double *s, double *s_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    // Source term for inertial forces in a rotating frame when solving in relative velocity.
		
	double RO_x = 0.0;
	double RO_y = 0.0;
    double W_z = 0.0;
    double dWdt_z = 0.0;
    
    //////////////
//     if (time < 1.0)
//         W_z = 0.0;
//     else if (time < 2.0) {
//         W_z = 0.2*(time-1.0);
//         dWdt_z = 0.2;
//     }
//     else
//         W_z = 0.2;
//     
//     if (time < 1.0) {
//         W_z = 0.1 + 0.1*time;
//         dWdt_z = 0.1;
//     }
//     else {
//         W_z = 0.2;
//     }
    //////////////
    
    error("Need to determine angular velocity and location of center of rotation.\n");

	for (int i = 0; i <ng; i++) {
		double OX = pg[0*ng+i];
		double OY = pg[1*ng+i];
        
		double r = udg[0*ng+i];
        double ru = udg[1*ng+i];
        double rv = udg[2*ng+i];
        
		s[1*ng+i] += r*W_z*W_z*(OX+RO_x) + 2.0*W_z*rv + r*dWdt_z*(OY+RO_y);
		s[2*ng+i] += r*W_z*W_z*(OY+RO_y) - 2.0*W_z*ru - r*dWdt_z*(OX+RO_x);
	}

    if (computeJacobian == 1) {
        for (int i = 0; i <ng; i++) {
            double OX = pg[0*ng+i];
            double OY = pg[1*ng+i];

            s_udg[0*ncu*ng+1*ng+i] +=  W_z*W_z*(OX+RO_x) + dWdt_z*(OY+RO_y);
            s_udg[0*ncu*ng+2*ng+i] +=  W_z*W_z*(OY+RO_y) - dWdt_z*(OX+RO_x);
            s_udg[1*ncu*ng+2*ng+i] += -2.0*W_z;
            s_udg[2*ncu*ng+1*ng+i] +=  2.0*W_z;
        }
    }
}
