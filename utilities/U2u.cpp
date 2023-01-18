#ifndef __U2U
#define __U2U

// Written by: C. Nguyen & P. Fernandez

void U2u(double * udg, double * UDG, double * pg, appstruct &app, int numPoints, int ncu, int nd)
{
    /* This function should be called only if ALEflag == 3 */

    int i, j;
    Int ALEflag = app.ALEflag;

    if (ALEflag != 3) {
        printf("U2u function was called with invalid value of ALEflag = %d.\n",ALEflag);
        printf("Execution will be terminated.\n");
        exit(-1);
    }

    /* Get gg */
    double * gg = &pg[(3 * nd + nd * nd) * numPoints];

    /* Compute u */
    for (i=0; i<ncu; i++)
        for (j=0; j<numPoints; j++)
            udg[j+i*numPoints] = UDG[j+i*numPoints] / gg[j];

}

#endif
