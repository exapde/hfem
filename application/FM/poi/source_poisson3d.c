void source_poisson3d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
    double pi = 3.14159265358979323846264;
	for (i = 0; i <ng*ncu; i++)		
        s[i] = 3*pi*pi*sin(pi*pg[i])*sin(pi*pg[ng+i])*sin(pi*pg[2*ng+i]);
	for (i = 0; i <ng*ncu*nc; i++)
		s_udg[i] = 0.0;
}

