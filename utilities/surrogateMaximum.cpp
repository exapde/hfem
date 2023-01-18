
// Written by: C. Nguyen & P. Fernandez

void surrogateMaximum(double * surrogMax, double * DsurrogMax_D_a1, double * DsurrogMax_D_a2, double a1, double a2, double alpha)
{
    /* This function computes a smooth surrogate for the maximum function, as well as its derivatives w.r.t. the two inputs */

    double actualMax;

    if (a1 > a2) {
        actualMax = a1;
    }
    else {
        actualMax = a2;
    }

    *surrogMax = actualMax + (1/alpha) * log(exp(alpha*(a1-actualMax)) + exp(alpha*(a2-actualMax)));
    *DsurrogMax_D_a1 = 1/(1+exp(alpha*(a2-a1)));
    *DsurrogMax_D_a2 = 1/(1+exp(alpha*(a1-a2)));

}
