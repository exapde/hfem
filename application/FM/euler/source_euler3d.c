void source_euler3d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
	for (i = 0; i <ng*ncu; i++)
		s[i] = 0.0;
	for (i = 0; i <ng*ncu*nc; i++)
		s_udg[i] = 0.0;
}

void sourceonly_euler3d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
	for (i = 0; i <ng*ncu; i++)
		s[i] = 0.0;
}

