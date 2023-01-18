
// Written by: C. Nguyen & P. Fernandez

void source_HIT_ns3d(double *s, double *s_udg, double *pg, double *udg, appstruct &app, double *param,
				 double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    // Forcing term for HIT with dynamic constant B * (k0/k). From 
    
	int i, j, k;
    double B = TBD, k0 = TBD, k;

	for (i = 0; i < ncu*ng; i++)
		s[i] = 0.0;
    
    for (j = 0; j < ng; j++) {
        k = (udg[1*ng+j]*udg[1*ng+j]+udg[2*ng+j]*udg[2*ng+j]+udg[3*ng+j]*udg[3*ng+j]) / (2*udg[0*ng+j]);
        s[1*ng+j] = B * (k0 / k) * udg[1*ng+j];
        s[2*ng+j] = B * (k0 / k) * udg[2*ng+j];
        s[3*ng+j] = B * (k0 / k) * udg[3*ng+j];
    }
    
	if (computeJacobian == 1) {
		int len = ng*ncu*nc;
        for (i = 0; i < len; i++)
            s_udg[i] = 0.0;
        
        for (j = 0; j < ng; j++) {
            k = (udg[1*ng+j]*udg[1*ng+j]+udg[2*ng+j]*udg[2*ng+j]+udg[3*ng+j]*udg[3*ng+j]) / (2*udg[0*ng+j]);
            dk_d0 = - (udg[1*ng+j]*udg[1*ng+j]+udg[2*ng+j]*udg[2*ng+j]+udg[3*ng+j]*udg[3*ng+j]) / (2*udg[0*ng+j]*udg[0*ng+j]);
            dk_d1 = udg[1*ng+j] / udg[0*ng+j];
            dk_d2 = udg[2*ng+j] / udg[0*ng+j];
            dk_d3 = udg[3*ng+j] / udg[0*ng+j];
            
            s_udg[0*ng*ncu+1*ng+j] = - B * dk_d0 * (k0 / (k*k)) * udg[1*ng+j];
            s_udg[0*ng*ncu+2*ng+j] = - B * dk_d0 * (k0 / (k*k)) * udg[2*ng+j];
            s_udg[0*ng*ncu+3*ng+j] = - B * dk_d0 * (k0 / (k*k)) * udg[3*ng+j];
            
            s_udg[1*ng*ncu+1*ng+j] = - B * dk_d1 * (k0 / (k*k)) * udg[1*ng+j] + B * (k0 / k);
            s_udg[1*ng*ncu+2*ng+j] = - B * dk_d1 * (k0 / (k*k)) * udg[2*ng+j];
            s_udg[1*ng*ncu+3*ng+j] = - B * dk_d1 * (k0 / (k*k)) * udg[3*ng+j];
            
            s_udg[2*ng*ncu+1*ng+j] = - B * dk_d2 * (k0 / (k*k)) * udg[1*ng+j];
            s_udg[2*ng*ncu+2*ng+j] = - B * dk_d2 * (k0 / (k*k)) * udg[2*ng+j] + B * (k0 / k);
            s_udg[2*ng*ncu+3*ng+j] = - B * dk_d2 * (k0 / (k*k)) * udg[3*ng+j];
            
            s_udg[3*ng*ncu+1*ng+j] = - B * dk_d3 * (k0 / (k*k)) * udg[1*ng+j];
            s_udg[3*ng*ncu+2*ng+j] = - B * dk_d3 * (k0 / (k*k)) * udg[2*ng+j];
            s_udg[3*ng*ncu+3*ng+j] = - B * dk_d3 * (k0 / (k*k)) * udg[3*ng+j] + B * (k0 / k);
        }
	}
}
