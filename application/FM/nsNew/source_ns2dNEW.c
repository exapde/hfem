
// Written by: C. Nguyen & P. Fernandez

void source_ns2dNEW(double *s, double *s_udg, double *pg, double *udg, appstruct &app, double *param,
				 double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
	int i;

	for (i = 0; i <ng*ncu; i++)
		s[i] = 0.0;

	if (computeJacobian == 1) {
		int len = ng*ncu*nc;
		for (i = 0; i < len; i++)
			s_udg[i] = 0.0;
	}
}
