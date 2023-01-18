#ifndef __UH2UH
#define __UH2UH

// Written by: C. Nguyen & P. Fernandez

void UH2uh(double * uh, double * UH, double * pg, appstruct &app, int numPoints, int nch, int nd)
{
    /* This function should be called only if ALEflag == 3 */

    int i, j;
    Int ALEflag = app.ALEflag;

    if (ALEflag != 3) {
        printf("UH2uh function was called with invalid value of ALEflag = %d.\n",ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get gg */
    double * gg = &pg[(3 * nd + nd * nd) * numPoints];

    /* Compute u */
    for (i=0; i<nch; i++)
        for (j=0; j<numPoints; j++)
            uh[j+i*numPoints] = UH[j+i*numPoints] / gg[j];

}

#endif
