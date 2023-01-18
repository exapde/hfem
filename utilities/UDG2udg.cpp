#ifndef __UDG2UDG
#define __UDG2UDG

// Written by: C. Nguyen & P. Fernandez

void UDG2udg(double * udg, double * UDG, double * pg, appstruct &app, int numPoints, int ncu, int ncq, int nc, int nd)
{
    /* This function converts the conservative variables from the undeformed domain to the deformed (physical) domain when the ALE formulation is used */
    /* This function should be called only if ALEflag == 2 or ALEflag == 3 */

    int i, j, k, m, n, len, sz2, sz3;
    double ggi;

    Int ALEflag = app.ALEflag;
    Int flag_q = app.flag_q;

    if (ALEflag != 2 && ALEflag != 3) {
        printf("UDG2udg function was called with invalid value of ALEflag = %d.\n",ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get Gg */
    double * Gg = &pg[3 * nd * numPoints];

    /* Get grad_g1 */
    double * grad_g1;
    if (ALEflag == 3){
        grad_g1 = &pg[(3 * nd + nd * nd + 2) * numPoints];
    }
    else {
        grad_g1 = new double[numPoints * nd];
        len = numPoints*nd;
        for (i = 0; i < len; i++) {
            grad_g1[i] = 0.0;
        }
    }


    /* Compute Ginv */
    double *Ginv = new double[numPoints * nd * nd];
    double * Gg11, * Gg21, * Gg31, * Gg12, * Gg22, * Gg32, * Gg13, * Gg23, * Gg33;
    double * Ginv11, * Ginv21, * Ginv31, * Ginv12, * Ginv22, * Ginv32, * Ginv13, * Ginv23, * Ginv33;
    if (nd == 2)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg12 = &Gg[2*numPoints];
        Gg22 = &Gg[3*numPoints];
        Ginv11 = &Ginv[0];
        Ginv21 = &Ginv[1*numPoints];
        Ginv12 = &Ginv[2*numPoints];
        Ginv22 = &Ginv[3*numPoints];

        for (i=0; i<numPoints; i++) {
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;
            Ginv11[i] =   Gg22[i] / ggi;
            Ginv12[i] = - Gg12[i] / ggi;
            Ginv21[i] = - Gg21[i] / ggi;
            Ginv22[i] =   Gg11[i] / ggi;
        }
    }
    else if (nd == 3)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg31 = &Gg[2*numPoints];
        Gg12 = &Gg[3*numPoints];
        Gg22 = &Gg[4*numPoints];
        Gg32 = &Gg[5*numPoints];
        Gg13 = &Gg[6*numPoints];
        Gg23 = &Gg[7*numPoints];
        Gg33 = &Gg[8*numPoints];

        Ginv11 = &Ginv[0];
        Ginv21 = &Ginv[1*numPoints];
        Ginv31 = &Ginv[2*numPoints];
        Ginv12 = &Ginv[3*numPoints];
        Ginv22 = &Ginv[4*numPoints];
        Ginv32 = &Ginv[5*numPoints];
        Ginv13 = &Ginv[6*numPoints];
        Ginv23 = &Ginv[7*numPoints];
        Ginv33 = &Ginv[8*numPoints];

        for (i=0; i<numPoints; i++) {
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;
            Ginv11[i] = (Gg22[i] * Gg33[i] - Gg23[i] * Gg32[i]) / ggi;
            Ginv12[i] = (Gg23[i] * Gg31[i] - Gg21[i] * Gg33[i]) / ggi;
            Ginv13[i] = (Gg21[i] * Gg32[i] - Gg22[i] * Gg31[i]) / ggi;
            Ginv21[i] = (Gg13[i] * Gg32[i] - Gg12[i] * Gg33[i]) / ggi;
            Ginv22[i] = (Gg11[i] * Gg33[i] - Gg13[i] * Gg31[i]) / ggi;
            Ginv23[i] = (Gg12[i] * Gg31[i] - Gg11[i] * Gg32[i]) / ggi;
            Ginv31[i] = (Gg12[i] * Gg23[i] - Gg13[i] * Gg22[i]) / ggi;
            Ginv32[i] = (Gg13[i] * Gg21[i] - Gg11[i] * Gg23[i]) / ggi;
            Ginv33[i] = (Gg11[i] * Gg22[i] - Gg12[i] * Gg21[i]) / ggi;
        }
    }


    /* Compute u */
    for (j=0; j < ncu; j++) {
        for (i = 0; i < numPoints; i++) {
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;
            udg[i + j * numPoints] = UDG[i + j * numPoints] / ggi;
        }
    }


    /* Compute q */
    double *q_tmp;
    if (flag_q == 1) {
        q_tmp = new double[numPoints * ncu * nd];

        sz2 = numPoints * ncu;
        for (k=0; k<nd; k++)
            for (j=0; j<ncu; j++)
                for (i=0; i<numPoints; i++) {
                    if (ALEflag == 3)
                        ggi = pg[(3 * nd + nd * nd) * numPoints + i];
                    else
                        ggi = 1.0;
                    q_tmp[i + j * numPoints + k * sz2] = UDG[i + j * numPoints + (1 + k) * sz2] / ggi +
                                                         UDG[i + j * numPoints] * grad_g1[i + k * numPoints];
                }

        sz3 = numPoints * nd;
        for (i=0; i<nd; i++)
            for (j=0; j<ncu; j++)
                for (k=0; k<numPoints; k++) {
                    udg[k + j * numPoints + (1+i) * sz2] = q_tmp[k + j * numPoints + 0 * sz2] * Ginv[k + 0 * numPoints + i * sz3];
                    for (m=1; m<nd; m++)
                        udg[k + j * numPoints + (1+i) * sz2] += q_tmp[k + j * numPoints + m * sz2]*Ginv[k + m * numPoints + i * sz3];
                }
    }

    delete[] Ginv;
    if (ALEflag == 2) {
        delete[] grad_g1;
    }
    if (flag_q == 1) {
        delete[] q_tmp;
    }
}

#endif
