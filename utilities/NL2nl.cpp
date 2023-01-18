#ifndef __NL2NL
#define __NL2NL

// Written by: C. Nguyen & P. Fernandez

void NL2nl(double * nl, double * NL, double * pg, appstruct &app, int numPoints, int nd)
{
    int i, j, k, sz2 = numPoints * nd;
    Int ALEflag = app.ALEflag;
    double nlNorm, ggi;

    if (ALEflag != 2 && ALEflag != 3) {
        printf("Invalid value of ALEflag = %d in NL2nl function.\n", ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get Gg */
    double * Gg = &pg[3 * numPoints * nd];

    /* Compute GinvT */
    double * GinvT = new double[numPoints * nd * nd];
    double * Gg11, * Gg21, * Gg31, * Gg12, * Gg22, * Gg32, * Gg13, * Gg23, * Gg33;
    double * GinvT11, * GinvT21, * GinvT31, * GinvT12, * GinvT22, * GinvT32, * GinvT13, * GinvT23, * GinvT33;
    if (nd == 2)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[1*numPoints];
        Gg12 = &Gg[2*numPoints];
        Gg22 = &Gg[3*numPoints];
        GinvT11 = &GinvT[0];
        GinvT21 = &GinvT[1*numPoints];
        GinvT12 = &GinvT[2*numPoints];
        GinvT22 = &GinvT[3*numPoints];

        for (i=0; i<numPoints; i++) {
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;
            GinvT11[i] =   Gg22[i] / ggi;
            GinvT12[i] = - Gg21[i] / ggi;
            GinvT21[i] = - Gg12[i] / ggi;
            GinvT22[i] =   Gg11[i] / ggi;
        }
    }
    else if (nd == 3)
    {
        Gg11 = &Gg[0];
        Gg21 = &Gg[numPoints];
        Gg31 = &Gg[2*numPoints];
        Gg12 = &Gg[3*numPoints];
        Gg22 = &Gg[4*numPoints];
        Gg32 = &Gg[5*numPoints];
        Gg13 = &Gg[6*numPoints];
        Gg23 = &Gg[7*numPoints];
        Gg33 = &Gg[8*numPoints];

        GinvT11 = &GinvT[0];
        GinvT21 = &GinvT[1*numPoints];
        GinvT31 = &GinvT[2*numPoints];
        GinvT12 = &GinvT[3*numPoints];
        GinvT22 = &GinvT[4*numPoints];
        GinvT32 = &GinvT[5*numPoints];
        GinvT13 = &GinvT[6*numPoints];
        GinvT23 = &GinvT[7*numPoints];
        GinvT33 = &GinvT[8*numPoints];

        for (i=0; i<numPoints; i++) {
            if (ALEflag == 3)
                ggi = pg[(3 * nd + nd * nd) * numPoints + i];
            else
                ggi = 1.0;
            GinvT11[i] = (Gg22[i] * Gg33[i] - Gg23[i] * Gg32[i]) / ggi;
            GinvT12[i] = (Gg13[i] * Gg32[i] - Gg12[i] * Gg33[i]) / ggi;
            GinvT13[i] = (Gg12[i] * Gg23[i] - Gg13[i] * Gg22[i]) / ggi;
            GinvT21[i] = (Gg23[i] * Gg31[i] - Gg21[i] * Gg33[i]) / ggi;
            GinvT22[i] = (Gg11[i] * Gg33[i] - Gg13[i] * Gg31[i]) / ggi;
            GinvT23[i] = (Gg13[i] * Gg21[i] - Gg11[i] * Gg23[i]) / ggi;
            GinvT31[i] = (Gg21[i] * Gg32[i] - Gg22[i] * Gg31[i]) / ggi;
            GinvT32[i] = (Gg12[i] * Gg31[i] - Gg11[i] * Gg32[i]) / ggi;
            GinvT33[i] = (Gg11[i] * Gg22[i] - Gg12[i] * Gg21[i]) / ggi;
        }
    }


    // Compute nl
    for (i=0; i<numPoints; i++) {
        if (ALEflag == 2)
            ggi = 1.0;
        else if (ALEflag == 3)
            ggi = pg[(3 * nd + nd * nd) * numPoints + i];

        for (j = 0; j < nd; j++) {
            nl[j * numPoints + i] = ggi * GinvT[0 * sz2 + j * numPoints + i] * NL[0 * numPoints + i];
            for (k=1; k<nd; k++)
                nl[j * numPoints + i] += ggi * GinvT[k * sz2 + j * numPoints + i] * NL[k * numPoints + i];
        }

        nlNorm = 0.0;
        for (j = 0; j < nd; j++)
            nlNorm += nl[j * numPoints + i] * nl[j * numPoints + i];
        nlNorm = sqrt(nlNorm);

        for (j = 0; j < nd; j++)
            nl[j * numPoints + i] = nl[j * numPoints + i] / nlNorm;
    }


    // Deallocate dynamic memory
    delete[] GinvT;
}

#endif
