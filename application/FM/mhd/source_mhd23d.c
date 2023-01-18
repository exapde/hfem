void source_mhd23d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
	double alpha1 = param[2];
	double alpha2 = sqrt(0.18*alpha1);
	double alpha  = (alpha1*alpha1)/(alpha2*alpha2);
	for (i = 0; i <ng*ncu; i++){
		s[i] = 0.0;
		if(i>=ng*(ncu-1)){
		  s[i] = -udg[i]*alpha;
		}
	}
	for (i = 0; i <ng*ncu*nc; i++){
		s_udg[i] = 0.0;
		if(i>=ng*(ncu-1)*(nc-1)){
		  s_udg[i] = -alpha;
		}
	}
}

void sourceonly_mhd23d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	int i;
	for (i = 0; i <ng*ncu; i++)
		s[i] = 0.0;
}

