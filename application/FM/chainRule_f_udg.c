
// Written by: C. Nguyen & P. Fernandez

void chainRule_f_udg(double *f_udg, double *f_mu, double *mu_udg, int ng, int nc, int ncu, int nd)
{
    int i, j, k, l;
    int sz1 = ng;
    int sz2 = sz1 * ncu;
    int sz3 = sz2 * nd;
    
    for (i=0; i<nc; i++)
        for (j=0; j<nd; j++)
            for (k=0; k<ncu; k++)
                for (l=0; l<ng; l++)
                    f_udg[l+k*sz1+j*sz2+i*sz3] += f_mu[l+k*sz1+j*sz2] * mu_udg[l+i*sz1];
    // f_udg: ng / ncu / nd / nc
    // f_mu: ng / ncu / nd
    // mu_udg: ng / nc
    
    // Alternative implementation (requires permuting f_udg, f_mu and mu_udg)
//     Int i;
//     Int na = ncu*nd, oneInt = 1, if_udg = ncu*nd*nc;
//     char chn = 'N';
//     double one = 1.0; 
//     
//     for (i=0; i<ng; i++) {
//         DGEMM(&chn, &chn, &na, &nc, &oneInt, &one, &f_mu[i*na], &na, &mu_udg[i*nc], 
//                         &oneInt, &one, &f_udg[i*if_udg], &na);
//     }
//     // f_udg: ncu / nd / nc / ng
//     // f_mu: ncu / nd / ng
//     // mu_udg: nc / ng
}
