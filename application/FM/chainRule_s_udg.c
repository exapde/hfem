
// Written by: C. Nguyen & P. Fernandez

void chainRule_s_udg(double *s_udg, double *s_mu, double *mu_udg, int ng, int nc, int ncu, int nd)
{
    int i, j, k;
    int sz1 = ng;
    int sz2 = sz1 * ncu;
    
    for (i=0; i<nc; i++)
        for (j=0; j<ncu; j++)
            for (k=0; k<ng; k++)
                s_udg[k+j*sz1+i*sz2] += s_mu[k+j*sz1] * mu_udg[k+i*sz1];
    // s_udg: ng / ncu / nc
    // s_mu: ng / ncu
    // mu_udg: ng / nc
    
    // Alternative implementation (requires permuting s_udg, s_mu and mu_udg)
//     Int i;
//     Int oneInt = 1, is_udg = ncu*nc;
//     char chn = 'N';
//     double one = 1.0; 
//     
//     for (i=0; i<ng; i++) {
//         DGEMM(&chn, &chn, &ncu, &nc, &oneInt, &one, &s_mu[i*ncu], &ncu, &mu_udg[i*nc], 
//                         &oneInt, &one, &s_udg[i*is_udg], &ncu);
//     }
//     // s_udg: ncu / nc / ng
//     // s_mu: ncu / ng
//     // mu_udg: nc / ng
}
