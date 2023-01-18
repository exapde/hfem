#ifndef __SOLUTIONUPDATE
#define __SOLUTIONUPDATE

// Written by: C. Nguyen & P. Fernandez

void getU(solstruct &sol, sysstruct &sys, meshstruct &mesh, appstruct &app, tempstruct &temp, Int* ndims,
        double* udg, double* duh_Vector, double alpha, Int ie)
{
    /* TODO: Use a different strategy to update U with quasi-Newton */
    Int nfe = ndims[2];
    Int npf = ndims[10];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];
    Int nch = ndims[23];
    Int ncf = ndims[75];

    Int ndf = npf*nfe;
    Int na = npv*ncu;
    Int i, nb, lenDinvF;

    Int isEDGelement = mesh.isEDGelement[ie];

    if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0)) {
        nb = ndf*nch;
    }
    else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {
        nb = ncf*nch;
    }
    if (app.hybrid == 0 || app.hybrid == 2) {
        lenDinvF = npv*ncu*ndf*nch;
    }
    else if (app.hybrid == 1) {
        lenDinvF = npv*ncu*ncf*nch;
    }
    
    Int inc = 1, n;
    char chn = 'N', charu = 'U';
    double one = 1.0, minusalpha = -alpha;

    double *DinvF = &sol.DinvF[ie*lenDinvF];
    double *DinvRu = &sol.DinvRu[ie*na];
    double *duh = &temp.uh[0];
    
    Int e = mesh.elementtype[ie];
    Int *cg2dg = &mesh.cg2dg[e][0];

    if (app.hybrid == 0 || (app.hybrid == 2 && isEDGelement == 0)) {
        for (int j=0; j<ndf; j++) {
            n = mesh.elcon[ie*ndf+j];
            for (int k=0; k<nch; k++)
                duh[j*nch+k] = duh_Vector[n*nch+k];
        }
    }
    else if (app.hybrid == 1 || (app.hybrid == 2 && isEDGelement == 1)) {
        for (int j=0; j<ncf; j++) {
            n = mesh.elcon[ie*ndf+cg2dg[j]];
            for (int k=0; k<nch; k++)
                duh[j*nch+k] = duh_Vector[n*nch+k];
        }
    }

    /* U <- U - alpha * (inv(D)*F) * DUH */
    DGEMV(&chn, &na, &nb, &minusalpha, &DinvF[0], &na, &duh[0], &inc, &one, &udg[0], &inc);

    /* U <- U + alpha * inv(D)* Ru */
    DAXPY(&na, &alpha, &DinvRu[0], &inc, &udg[0], &inc);
}

void updateSol(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, double* udg, double* uh, double* duh, double alpha)
{
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];

    /* update UH */
    for (Int i = 0; i< sol.UH.size(); i++)
        uh[i] += alpha*duh[i];
    
    /* update UDG */
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for
        for (Int ie=0; ie<ne; ie++) {
            /* compute U */
            getU(sol, sys, mesh, app, temps[this_thread], &ndims[0], &udg[ie*npv*nc], &duh[0], alpha, ie);

            /* Compute Q */
            if (app.flag_q==1)
                getQ(elems[this_thread], mesh, master, app, sol, temps[this_thread], 
                     &udg[ie*npv*nc], &udg[ie*npv*nc+npv*ncu], &uh[0], &sol.SH[ie*npv*nc+npv*ncu], ie);
        }
// // //     }
}

#ifdef  HAVE_MPI
void updateSolMPI(sysstruct &sys, elemstruct* elems, meshstruct &mesh, masterstruct &master,
        solstruct &sol, appstruct &app, tempstruct* temps, Int* ndims, double* udg, double* uh, double* duh, double alpha)
{
    // TODO: Send the vectors as soon as they are computed, instead of in this function

    Int i, n, ri, neighbor, inc = 1;
    Int bsz = sys.blkSize;
    Int ne = ndims[5];
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];

    /* initialize request counter */
    Int request_counter = 0;
    
    /* copy some portion of duh to buffsend */
    for (i=0; i<sys.nentsend; i++) {
        ri = sys.entsend[i];
        DCOPY(&bsz, &duh[bsz*ri], &inc, &sys.buffsend[bsz*i], &inc);
    }

    /* non-blocking send */
    Int nsend, psend = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nsend = sys.entsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&sys.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    Int nrecv, precv = 0;
    for (n=0; n<sys.nnbsd; n++) {
        neighbor = sys.nbsd[n];
        nrecv = sys.entrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&sys.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &sys.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    /* update interior dofs of UH while sending and receiving data */
    Int indof = sys.BJ_nrows*bsz;
    for (i = 0; i< indof; i++)
        uh[i] += alpha*duh[i];

    /* update interior dofs of UDG while sending and receiving data */
    Int nein = sys.elempartpts[0];
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        int this_thread = 0;
// // //         #pragma omp for
        for (Int ie=0; ie<nein; ie++) {
            /* compute U */
            getU(sol, sys, mesh, app, temps[this_thread], &ndims[0], &udg[ie*npv*nc], &duh[0], alpha, ie);

            /* Compute Q */
            if (app.flag_q==1)
                getQ(elems[this_thread], mesh, master, app, sol, temps[this_thread],
                     &udg[ie*npv*nc], &udg[ie*npv*nc+npv*ncu], &uh[0], &sol.SH[ie*npv*nc+npv*ncu], ie);
        }
// // //     }
    
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, sys.requests, sys.statuses);

    /* copy buffrecv to duh */
    for (i=0; i<sys.nentrecv; i++) {
        ri = sys.entrecv[i];
        DCOPY(&bsz, &sys.buffrecv[bsz*i], &inc, &duh[bsz*ri], &inc);
    }
    
    /* update remaining dofs of UH */
    for (i=indof; i<sol.UH.size(); i++)
        uh[i] += alpha*duh[i];

    /* update remaining dofs of UDG */
// // //     #pragma omp parallel num_threads(sys.noThreads)
// // //     {
// // //         int this_thread = omp_get_thread_num();
        this_thread = 0;
// // //         #pragma omp for
        for (Int ie=nein; ie<ne; ie++) {
            /* compute U */
            getU(sol, sys, mesh, app, temps[this_thread], &ndims[0], &udg[ie*npv*nc], &duh[0], alpha, ie);

            /* Compute Q */
            if (app.flag_q==1)
                getQ(elems[this_thread], mesh, master, app, sol, temps[this_thread],
                     &udg[ie*npv*nc], &udg[ie*npv*nc+npv*ncu], &uh[0], &sol.SH[ie*npv*nc+npv*ncu], ie);
        }
// // //     }
}
#endif

#endif
