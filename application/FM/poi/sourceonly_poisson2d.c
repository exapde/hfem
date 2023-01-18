void sourceonly_poisson2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
    double pi = 3.14159265358979323846264;
	for (i = 0; i <ng*ncu; i++)
		s[i] = 2*pi*pi*sin(pi*pg[i])*sin(pi*pg[ng+i]);
}

