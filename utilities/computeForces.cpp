#ifndef __COMPUTEFORCES
#define __COMPUTEFORCES

// Written by: C. Nguyen & P. Fernandez

void computeForces(double *forces, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims, Int bouID)
{
    Int j, k, n, ie, id, ig, is;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int nfe = ndims[2];
    Int npv = ndims[9];
    Int npf = ndims[10];
    Int nqfR = master.nqfR;
    Int nc = ndims[19];
    Int ncu = ndims[20];
    Int nch = ndims[23];
    Int ndf = npf*nfe;
    
    Int nein;
    if (sys.nproc > 1)
        nein = sys.elempartpts[0] + sys.elempartpts[1];
    else
        nein = ndims[5];
    
    Int *bf;
    double *gwfcR = &master.gwfcR[0];

    for (id = 0; id < nd; id++)
        forces[id] = 0.0;

    for (ie = 0; ie < nein; ie++) {
        double *udg = &sol.UDG[ie*npv*nc];
        double *udg_ref = &sol.UDGi[ie*npv*nc];
        double *av = &sol.avField_p1CG[ie*npv];
//         double *av = &sol.avField[ie*npv];
        double *uh = &temp.uh[0];
        double *uh_ref = &temp.uh_ref[0];
        double *pn = &mesh.dgnodes[ie*npv*ncd];
        for (j=0; j<ndf; j++) {
            n = mesh.elcon[ie*ndf+j];
            for (k=0; k<nch; k++) {
                uh[k*ndf+j] = sol.UH[n*nch+k];
                uh_ref[k*ndf+j] = sol.UHi[n*nch+k];
            }
        }
//         facegeom(mesh, master, temp, ndims, ie);
        facegeom(&temp.pfR[0], &temp.JfR[0], &temp.XxfR[0], &temp.jacfR[0], &temp.nlR[0], &temp.pf[0], &master.shapftR[0], &ndims[0], master.nqfR);
//         faceflux(mesh, master, app, sol, temp, ndims, ie);
        faceflux(&temp.fhR[0], &temp.fhR_u[0], &temp.fhR_uh[0], &temp.ufR[0], &temp.uf_refR[0], 
                 &temp.uhR[0], &temp.uh_refR[0], &temp.avfR[0], &pn[0], &udg[0], &udg_ref[0], 
                 &uh[0], &uh_ref[0], &av[0], &temp.pfR[0], &temp.nlR[0], &master.shapftR[0], 
                 mesh, master, app, temp, &ndims[0], master.nqfR, 0);
        
        bf = &mesh.bf[ie*nfe];
        
        for (is = 0; is < nfe; is++) {
                if (bf[is] == -bouID) {
                    
                    for (ig = 0; ig < nqfR; ig++)
                        for (id = 0; id < nd; id++) {
                            forces[id] += gwfcR[ig] * temp.jacfR[is*nqfR+ig] * temp.fhR[(id+1)*nfe*nqfR+is*nqfR+ig];
                        }
                }
            }
    }
}

void computeMomentum(double *momentum, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims)
{
    Int ie, id, ig;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int npv = ndims[9];
    Int nqvR = master.nqvR;
    Int nc = ndims[19];
    Int ncu = ndims[20];
    
    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    Int nein;
    if (sys.nproc > 1)
        nein = sys.elempartpts[0] + sys.elempartpts[1];
    else
        nein = ndims[5];
    
    double *gwvlR = &master.gwvlR[0];
    double *shapvtR = &master.shapvtR[0];
    double *udgR = &temp.udgp[0];
    double *udg;

    for (id = 0; id < nd; id++)
        momentum[id] = 0.0;
    
    for (ie = 0; ie < nein; ie++) {
        udg = &sol.UDG[ie*npv*nc];

//         elementgeom(mesh, master, temp, ndims, ie);
        elementgeom(&temp.pR[0], &temp.J_R[0], &temp.XxR[0], &temp.jacR[0], &mesh.dgnodes[ie*npv*ncd], &master.shapvtR[0], &ndims[0], master.nqvR);
        
        /* Compute solution at Gauss points */
        DGEMM(&chn, &chn, &nqvR, &nc, &npv, &one, shapvtR, &nqvR, udg,
                        &npv, &zero, udgR, &nqvR);
        
        for (ig = 0; ig < nqvR; ig++)
            for (id = 0; id < nd; id++)
                momentum[id] += gwvlR[ig] * temp.jacR[ig] * udgR[(id+1)*nqvR+ig];
    }
}

void computeMeasure(double *measure, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims)
{
    Int ie, ig;
    Int ncd = ndims[1];
    Int npv = ndims[9];
    Int nqvR = master.nqvR;
    
    Int nein;
    if (sys.nproc > 1)
        nein = sys.elempartpts[0] + sys.elempartpts[1];
    else
        nein = ndims[5];
    
    double *gwvlR = &master.gwvlR[0];
    
    *measure = 0.0;
    
    for (ie = 0; ie < nein; ie++) {
//         elementgeom(mesh, master, temp, ndims, ie);
        elementgeom(&temp.pR[0], &temp.J_R[0], &temp.XxR[0], &temp.jacR[0], &mesh.dgnodes[ie*npv*ncd], &master.shapvtR[0], &ndims[0], nqvR);
        
        for (ig = 0; ig < nqvR; ig++)
            *measure += gwvlR[ig] * temp.jacR[ig];
    }
}

void computeElementMeasure(sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims)
{
    Int ie, ig;
    Int nd = ndims[0];
    Int ncd = ndims[1];
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nqvR = master.nqvR;
    
    double *gwvlR = &master.gwvlR[0];
    
    if (mesh.elemMeasure.size() != ne)
        mesh.elemMeasure.resize(ne);
    if (mesh.hAvg.size() != ne)
        mesh.hAvg.resize(ne);
    
    for (ie = 0; ie < ne; ie++) {
//         elementgeom(mesh, master, temp, ndims, ie);
        elementgeom(&temp.pR[0], &temp.J_R[0], &temp.XxR[0], &temp.jacR[0], &mesh.dgnodes[ie*npv*ncd], &master.shapvtR[0], &ndims[0], master.nqvR);
        
        mesh.elemMeasure[ie] = 0.0;
        for (ig = 0; ig < nqvR; ig++)
            mesh.elemMeasure[ie] += gwvlR[ig] * temp.jacR[ig];
        
        mesh.hAvg[ie] = pow(mesh.elemMeasure[ie], 1.0 / ((double) nd));
    }
}

#ifdef  HAVE_MPI
void computeForcesMPI(double *forces, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims, Int bouID)
{
    Int nd = ndims[0];
    double *localForcesContribution = new double[nd];
    
    computeForces(&localForcesContribution[0], sys, elem, mesh, master, sol, app, temp, ndims, bouID);
    
    for (Int id = 0; id < nd; id++) {
        forces[id] = 0.0;
        MPI_Allreduce(&localForcesContribution[id], &forces[id], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    delete[] localForcesContribution;
}

void computeMomentumMPI(double *momentum, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims)
{
    Int nd = ndims[0];
    double *localMomentumContribution = new double[nd];
    
    computeMomentum(&localMomentumContribution[0], sys, elem, mesh, master, sol, app, temp, ndims);
    
    for (Int id = 0; id < nd; id++) {
        momentum[id] = 0.0;
        MPI_Allreduce(&localMomentumContribution[id], &momentum[id], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    delete[] localMomentumContribution;
}

void computeMeasureMPI(double *measure, sysstruct &sys, elemstruct &elem, 
        meshstruct &mesh, masterstruct &master, solstruct &sol, appstruct &app, 
        tempstruct &temp, Int* ndims)
{
    double localMeasureContribution[0];
    *measure = 0.0;
    
    computeMeasure(&localMeasureContribution[0], sys, elem, mesh, master, sol, app, temp, ndims);
    
    MPI_Allreduce(&localMeasureContribution[0], &measure[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
#endif

#endif
