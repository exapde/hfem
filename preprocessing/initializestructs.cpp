#ifndef __INITIALIZESTRUCTS
#define __INITIALIZESTRUCTS

void cleardmdstruct(dmdstruct &dmd) 
{        
    dmd.nbsd.clear();    // neighboring cpus
    dmd.intelem.clear(); // nonoverlapping global elements
    dmd.intent.clear();  // nonoverlapping global entities
    dmd.elempart.clear(); // overlapping global elements
    dmd.elempartpts.clear(); // classifiers of overlapping global elements: (interior, interface, exterior)
    dmd.entpart.clear();  // overlapping global entities
    dmd.entpartpts.clear(); // classifiers of overlapping global entities: (interior, interface, exterior)
    dmd.elemrecv.clear();  // local elements received from neighboring cpus
    dmd.elemrecvpts.clear(); // classifiers of local elements received from neighboring cpus: (# elements from ncpu1, ...) 
    dmd.elemsend.clear();  // local elements sent to neighboring cpus
    dmd.elemsendpts.clear(); // classifiers of local elements sent to neighboring cpus: (# elements to ncpu1, ...)    
    dmd.entrecv.clear();    // local entities received from neighboring cpus
    dmd.entrecvpts.clear(); // classifiers of local entities received from neighboring cpus: (# entities from ncpu1, ...)          
    dmd.entsend.clear();  // local entities sent to neighboring cpus
    dmd.entsendpts.clear(); // classifiers of local entities sent to neighboring cpus: (# elements to ncpu1, ...)    
    dmd.vecrecv.clear();  // local vectors received from neighboring cpus to perform matrix-vector product
    dmd.vecrecvpts.clear(); // classifiers of local vectors received from neighboring cpus: (# vectors from ncpu1, ...)          
    dmd.vecsend.clear();    // local vectors sent to neighboring cpus to perform matrix-vector product
    dmd.vecsendpts.clear(); // classifiers of local vectors sent to neighboring cpus: (# vectors to ncpu1, ...)             
    dmd.matrecv.clear();    // local matrices received from neighboring cpus to construct the preconditioner
    dmd.matrecvpts.clear(); // classifiers of local matrices received from neighboring cpus: (# matrices from ncpu1, ...)                
    dmd.matsend.clear();   // local matrices sent to neighboring cpus to construct the preconditioner
    dmd.matsendpts.clear();  // classifiers of local matrices sent to neighboring cpus: (# matrices to ncpu1, ...)                                         
    dmd.rowent2elem.clear(); // global entity-to-element connectivities
    dmd.colent2elem.clear(); // global entity-to-element connectivities
    dmd.rowent2ent.clear();  // global entity-to-entity connectivities
    dmd.colent2ent.clear();  // global entity-to-entity connectivities
    dmd.bcrs_rowent2elem.clear(); // local entity-to-element connectivities
    dmd.bcrs_colent2elem.clear(); // local entity-to-element connectivities
    dmd.bcrs_rowent2ent.clear();  // local entity-to-entity connectivities
    dmd.bcrs_colent2ent.clear();  // local entity-to-entity connectivities    
    dmd.elemmap.clear();  // element reordering 
    dmd.entmap.clear();  // entity reordering 
    dmd.ent2ind.clear();  // global-to-local entity mapping  
    dmd.elcon.clear();    // local element-to-entity connectivities
    dmd.t2f.clear();      // local elemeent-to-face connectivities
    //dmd.t2t;      // local elemeent-to-element  connectivities
}

void initializeStructs(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, elemstruct &elems, sysstruct &sys, tempstruct &temps)
{
    Int nd = app.nd;
    Int ncd = app.ncd;
    Int nc = app.nc;
    Int ncu = app.ncu;
    Int ncq = app.ncq;
    Int ncp = app.ncp;
    Int nch = app.nch;
    Int nco = app.nco;    
    Int ne = mesh.ne;
    Int nf = mesh.nf;
    Int nv = mesh.nv;
    Int ndh = mesh.ndofuh;
    Int nfe = mesh.nfemax;
    Int nve = mesh.nvemax;
    Int nvf = mesh.nvfmax;
    Int npe = mesh.npemax;
    Int npf = mesh.npfmax;
    Int nme = mesh.nmemax;
    Int nmf = mesh.nmfmax;
    Int nge = mesh.ngemax;
    Int ngeR = mesh.ngeRmax;   
    Int ngf = mesh.ngfmax;    
    Int ngfR = mesh.ngfRmax;    
    Int ndf = mesh.ndfmax;
    Int ncf = mesh.ncfmax;
    Int cbsr_nrows = sys.numEntities;
    Int cbsr_nblks = sys.numBlocks;
    Int maxBlocksPerRow = sys.maxBlocksPerRow;
    Int blkSize = sys.blkSize;
    Int BJ_nrows = sys.BJ_nrows;    
    Int nproc = sys.nproc;    
    Int nentrecv = sys.nentrecv;
    Int nentsend = sys.nentsend;    
    Int nelemrecv = sys.nelemrecv;
    Int nelemsend = sys.nelemsend;    
    Int nvecrecv = sys.nvecrecv;
    Int nvecsend = sys.nvecsend;                
    Int nmatrecv = sys.nmatrecv;
    Int nmatsend = sys.nmatsend;  
    Int nnbsd     = sys.nnbsd;
    
    Int i, j, k, l;

    /* initialize app struct */
    app.tdep = app.flag[0];             /* flag to determine unsteady or steady problems */
    app.wave = app.flag[1];             /* flag to determine wave or non-wave problems */
    app.alag = app.flag[2];             /* flag to determine augmented lagrangian */
    app.adjoint = app.flag[3];          /* flag to determine adjoint or primal problems */
    app.linearProblem = app.flag[4];    /* flag to determine linear or nonlinear problems */
    app.flag_q = app.flag[5];           /* flag to determine if the Q equation is also discretized  */
    app.flag_p = app.flag[6];           /* flag to determine if the P equation is also discretized  */
    app.flag_g = app.flag[7];           /* flag to determine if the GCL equation is also discretized  */
    app.debugmode = app.flag[8];
    app.overlappinglevel = app.flag[9];
    app.reuseOrdering = app.flag[10];
    app.reuseJacobian = app.flag[11];
    app.reuseResidual = app.flag[12];        
    app.jacobianStep = app.flag[13];  /* reuse Jacobian for every jacovianStep time steps */
    app.orderingStep = app.flag[14];  /* reuse Ordering for every jacovianStep time steps */    
    app.quasiNewton = app.flag[15];    
            
    app.fc_u = app.factor[0];           /* factor when discretizing the time derivative of the U equation */
    app.fc_q = app.factor[1];           /* factor when discretizing the time derivative of the Q equation */
    app.fc_p = app.factor[2];           /* factor when discretizing the time derivative of the P equation */
    app.time = 0.0;           /* initialize time */
    if (app.factor.size()>3)
        app.time = app.factor[3];   
    else
        app.time = 0;
    if (app.factor.size()>4)
        app.dtfc = app.factor[4];       /* needed only for augmented Lagrangian */
    else
        app.dtfc = 0;
    if (app.factor.size()>5)
        app.alpha = app.factor[5];      /* needed only for augmented Lagrangian */
    else
        app.alpha = 0;
    
    app.hybrid = app.problem[0];        /* 0: HDG; 1: EDG; 2: IEDG, HEDG */
    app.appname = app.problem[1];       /* 0: Euler; 1: Compressible Navier-Stokes; etc. */
    app.temporalscheme = app.problem[2];
    app.dirkOrder = app.problem[3];    /* temporal accuracy order */
    app.BDFsteps = app.problem[3];
    app.dirkStage = app.problem[4];    /* DIRK stages */
    app.convStabMethod = app.problem[5];  // Flag for convective stabilization tensor. 0: Constant tau, 1: Lax-Friedrichs; 2: Roe.
    app.diffStabMethod = app.problem[6];  // Flag for diffusive stabilization tensor. 0: No diffusive stabilization.
    app.rotatingFrame = app.problem[7];   // Flag for rotating frame. Options: 0: Velocities are respect to a non-rotating frame. 1: Velocities are respect to a rotating frame.
    app.viscosityModel = app.problem[8];  // Flag for viscosity model. 0: Constant kinematic viscosity; 1: Sutherland's law    
    app.SGSmodel = app.problem[9];        // Flag for sub-grid scale (SGS) model. Only available for 3D solver.
                                        //  0: No SGS model. 1: Static Smagorinsky/Yoshizawa/Knight model. 
                                        //  2: Static WALE/Yoshizawa/Knight model. 3: Static Vreman/Yoshizawa/Knight model.
                                        //  4: Dynamic Smagorinsky/Yoshizawa/Knight model.        
    app.ALEflag = app.problem[10];                    // Flag for Arbitrary Lagrangian-Eulerian (ALE) formulation. 0: No ALE; 1: Translation; 2: Translation + rotation; 3: Translation + rotation + deformation
    app.AVflag = app.problem[11];                     // Flag for artificial viscosity. 0: No artificial viscosity; 1: Homogeneous artificial viscosity (C. Nguyen's formulation); 2: Hypersonic homogeneous artificial viscosity (C. Nguyen's formulation)
                                        //                                3: Isotropic artificial viscosity (D. Moro's formulation). 4: Latest version of the model (taking the best of all previous models)
                                        //                                8: Density smoothness sensor (Per's approach)    
    app.linearSolver = app.problem[12];  /* 0: direct solver; 1: gmres; etc. */
    
    app.rampFactor = 1.0;               // Ramp factor for artificial viscosity flux
    
// app.flag   = [app.tdep app.wave app.alag app.adjoint app.linearproblem app.flg_q app.flg_p app.flg_g... 
//               app.debugmode app.overlappinglevel app.preconditionertype app.preconditionerside...
//               app.reuseOrdering app.reuseJacobian app.reuseResidual app.jacobianStep app.orderingStep...
//               app.quasiNewton app.quasiNewtonAccuracy app.orthogMethod app.reorderMethod...
//               app.schurImplementation app.matvecImplementation app.precSolveImplementation...
//               app.precPrecision app.matvecPrecision app.orthogPrecision app.adaptiveGMREStol app.flag];
             
    if (app.tdep == 1) {
        app.DIRKnotConverged.resize(app.dt.size()*app.dirkStage);
        for (i = 0; i < app.DIRKnotConverged.size(); i++)
            app.DIRKnotConverged[i] = 0;
        app.BDFnotConverged.resize(app.dt.size());
        for (i = 0; i < app.BDFnotConverged.size(); i++)
            app.BDFnotConverged[i] = 0;
    }        
        
    
    /////////////////// SOL STRUCT ///////////////////    
    Int n1, nein;
    if (nproc == 1)
        nein = ne;
    else if (nproc > 1)
        nein = sys.elempartpts[0]+sys.elempartpts[1];
    if (app.hybrid == 0 || app.hybrid == 2)
        n1 = nch*npf*nfe*npe*ncu*ne;        // This is more than necessary for IEDG
    else if (app.hybrid == 1)
        n1 = nch*ncf*npe*ncu*ne;
    
    sol.UDG2Write.resize(npe*ncu);       // Only required if writeQflag = 1.
    sol.SH.resize(npe*nc*ne);
    for (i=0; i<sol.SH.size(); i++)
        sol.SH[i] = 0.0;
    sol.DinvRu.resize(npe*ncu*ne);
    sol.DinvF.resize(n1);
    

    
    /////////////////// SYS STRUCT ///////////////////            
    sys.NewtonTol = app.solversparam[0];         // Newton tolerance for non-linear system
    sys.tol = app.solversparam[1];               // GMRES tolerance for linear system    
    sys.NewtonMaxiter = app.problem[13];          // Maximum number of Newton/quasi-Newton iterations
    sys.trueNewtonMaxiter = app.problem[13];         // Maximum number of Newton iterations
    sys.maxiter = 1000;//app.problem[14];             // Maximum number of GMRES iterations
    sys.restart = 200;//app.problem[15];              // Parameter k in GMRES(k)
    sys.preconditioner = app.flag[16];         // 0: Restricted additive Schwarz (RAS), 1: Subdomain-wise block Jacobi (BJ), 2: Entity-wise block Jacobi (BJ)
    sys.preconditionerSide = app.flag[17];     // -1: No preconditioner, 0: Left preconditioner. 1: Right preconditioner
    sys.quasiNewtonAccuracy = app.flag[18];        // Format in which Dinv and K are stored for quasi-Newton. 0: Single precision; 1: Double precision
    sys.orthogMethod = app.flag[19];           // Orthogonalization method in GMRES. 0: MGS, 1: ICGS, 2: IMGS (only MGS and ICGS are available in MPI code), 3: ICGS without convergence guarantee
    sys.reorderMethod = app.flag[20];          // 0: No reordering. 1: Approximate (inexact) MDF. 2: Exact MDF. 3: Approximate (inexact) MDF with constraints (only valid for parallel BJ preconditioner). 4: Exact MDF with constraints (only valid for parallel BJ preconditioner)    
    
    sys.print = 1;                      // 0: No info is printed. 1: Summary of CPU-time breakdown is printed. 2: Detailed CPU-time breakdown is printed. 3: GMRES residual is also printed.
    sys.schurImplementation = app.flag[21];        // Implementation flag for Schur complement. Options: 0 and 1
    sys.matvecImplementation = app.flag[22];       // Implementation flag for matrix-vector product. Options: 0, 1, 2 (finite differences for Gateaux derivative)
    sys.precSolveImplementation = app.flag[23];    // Implementation flag for BILU(0) solve. Only valid for RAS and subdomain-wise BJ. Options: 0 and 1
    sys.precPrecision = app.flag[24];              // Precision for preconditioner solve. Options: 0: Single precision. 1: Double precision
    sys.matvecPrecision = app.flag[25];            // Precision for matrix-vector product. Options: 0: Single precision. 1: Double precision
    sys.orthogPrecision = app.flag[26];            // Precision for orthogonalization. Options: 0: Single precision. 1: Double precision
    sys.adaptiveGMREStol = app.flag[27];           // Adaptive GMRES tolerance for nonlinear problems. Options: 1: Adaptive method. 0: Non-adaptive method. Note the adaptive strategy is not used for linear PDEs regardless the value of this flag.    
    
    sys.linearSolvesClocks = 0;         // Time spent in lear solves since the Jacobian matrix was computed the last time (only applies for quasi-Newton)
    sys.computeGlobalEnt2entWeight = 0; // Only applies for MPI code. If set to 1, the Frobenius norm of the Schur matrix blocks are computed, stored in a file and the execution is terminated.
    sys.linearProblem = app.linearProblem;  // Required for adaptive GMRES tolerance (in nonlinear problems) based on rNorm vs. NewtonTol    
    if (sys.adaptiveGMREStol == 1 && sys.linearProblem == 0)
        sys.adaptGMREStol = 1;      // Only applies for linear problems with adaptiveGMREStol == 1. Has two meanings:
                                    // Before executing GMRES routine: If adapting GMRES tolerance is allowed.
                                    // After executing GMRES routine: If GMRES tolerance was adapted in the previous linear solve
    else
        sys.adaptGMREStol = 0;
    sys.robustMode = 0;                
        
    //////////////////////////////////////////////


    // For quasi-Newton:
    if (app.quasiNewton == 1) {
        sol.ipivD.resize(npe*ncu*ne);
        if (sys.quasiNewtonAccuracy == 0) {
            sol.DinvFloat.resize(npe*ncu*npe*ncu*ne);
            sol.Kfloat.resize(n1);
        }
        else if (sys.quasiNewtonAccuracy == 1) {
            sol.Dinv.resize(npe*ncu*npe*ncu*ne);
            sol.K.resize(n1);
        }
    }
    
    // For AV model No. 6:
    if (app.AVflag == 6 || app.AVflag == 8 || app.AVflag == 9) {
        sol.avField_DG.resize(npe*ne);
        sol.avField_p1CG.resize(npe*ne);
        sol.avField.resize(npe*ne);
    }
    else if (app.AVflag == 10) {
        sol.avField_DG.resize(3*npe*ne);
        sol.avField_p1CG.resize(3*npe*ne);
        sol.avField.resize(3*npe*ne);
    }
    
//     if (app.AVflag == 8) {
//         Int projporder = porder - 1;
//         master.projLowP.resize(npe*npe);
//         getProjectionMatrix(&ndims[0], projporder, &master.projLowP[0]);
//     }
    
//     // For accelerated residual solver:
//     if (app.linearProblem == 2) {        
//         sol.Un.resize(npe*ncu*ne+nch*ndh);
//         sol.Um.resize(npe*ncu*ne+nch*ndh);
//         sol.Vn.resize(npe*ncu*ne+nch*ndh);
//         sol.R.resize(npe*ncu*ne+nch*ndh);
//     }    
    
    // Vectors for matrix-free matrix-vector product:
    if (sys.matvecImplementation == 2) {
        sol.UDG4FD.resize(npe*nc*ne);
        sol.UH4FD.resize(nch*ndh);
    }
    
//     // Average solution:
//     sol.UDG_avg.resize(npe*nc*ne);
//     sol.UH_avg.resize(nch*ndh);
    
//     // Initial guess of Minimal Residual algorithm:
//     sol.UDG_initMR.resize(npe*nc*ne);
//     sol.UH_initMR.resize(nch*ndh);
//     
    
    // Check inputs:
    if (nproc == 1 && (sys.preconditioner != -1 && sys.preconditioner != 0))
        error("Preconditioner not valid for serial version of the code.\n");
    
    if ((sys.precPrecision == 0) && (sys.preconditioner == 0 || sys.preconditioner == 1) && (sys.reorderMethod == 3 || sys.reorderMethod == 4))
        error("MDF with constraints not implemented for BILU0 in mixed-precision.\n");
    
    if ((sys.nproc == 1 || sys.preconditioner != 0) && (sys.reorderMethod == 3 || sys.reorderMethod == 4))
        error("MDF with constraints only implemented for parallel RAS preconditioner.\n");
    
    if (app.dirkStage != 1 && app.dirkStage != 2 && app.dirkStage != 3)
        error("DIRK stage not implemented yet %d.\n");
    
    if (app.adjoint != 0)
        error("Adjoint solver not implemented yet.\n");
    
    if (nproc == 1 && sys.matvecPrecision == 0)
        error("Single precision matrix-vector product not implemented in serial version of the code.\n");
    
    if (nproc == 1 && sys.orthogPrecision == 0)
        error("Single precision orthogonalization not implemented in serial version of the code.\n");
    
    if (nproc > 1 && sys.orthogMethod == 2)
        error("IMGS orthogonlization not implemented in MPI code yet.\n");
    
    if (sys.orthogMethod != 0 && sys.orthogMethod != 1 && sys.orthogMethod != 2 && sys.orthogMethod != 3)
        error("Invalid orthogonalization method.\n");
    
    if (sys.preconditioner != -1 && sys.preconditioner != 0 && sys.preconditioner != 1 && sys.preconditioner != 2)
        error("Invalid preconditioner.\n");
    
    if (sys.preconditionerSide != 0 && sys.preconditionerSide != 1)
        error("Invalid preconditioner side.\n");
    
    if (app.quasiNewton != 0 && app.quasiNewton != 1)
        error("Quasi-Newton flag has invalid value.\n");
    
    if (app.rotatingFrame != 0 && app.rotatingFrame != 1)
        error("Rotating frame flag has invalid value.\n");
    
    if (app.convStabMethod != 0 && app.convStabMethod != 1 && app.convStabMethod != 2 && app.convStabMethod != 3 && app.convStabMethod != 4)
        error("Flag for stabilization of convective operator has invalid value.\n");

    if (app.diffStabMethod != 0 && app.diffStabMethod != 1 && app.diffStabMethod != 2)
        error("Flag for stabilization of difussive operator has invalid value.\n");
    
    if (app.adjoint != 0 && app.adjoint != 1)
        error("Adjoint flag has invalid value.\n");
    
            
    // Get length for elem and temp structures:
    Int elemtempLen = cbsr_nblks + cbsr_nrows + nge*ncd + nge*nd*nd + nge + nge*nd*nd + npf*nfe*ncd + 
                ngf*nfe*ncd + ngf*nfe*nd*(nd-1) + ngf*nfe*nd + ngf*nfe + ngf*nfe + nge*nc + nge*nc + nge*ncu*nd + 
                nge*ncu*nd*nc + nge*ncu + nge*ncu*nc + nge*ncu*(nd+1)*nc + ndf*nch + ndf*nch + 
                ndf*nc + ndf*nc + ngf*nfe*nc + ngf*nfe*nc + ngf*nfe*nch + ngf*nfe*nch + 3*nge + 3*npf*nfe + 
                3*ngf*nfe + ngf*nfe*nch + ngf*nfe*nch*nc + ngf*nfe*nch*nch + ndf*nch + npf*ndf*nch*nc + 
                npe*npe*nch*nc + npe*ncu*ndf*nch + npf*npf*nfe*nch*nch + npf*npf*nfe + ngf*nch + 
                ngf*nch*nc + ngf*nch*nch + ngf*ncd + ngf*nc + ngf*nch + ngf*nd + ngf*nc + ngf*nch + 3*ngf + 
                npe*ncu + npe*npe*ncu*nc + npe*ncu*ndf*nch + nch*ndf + nch*ndf + 
                nch*ndf*npe*nc + nch*ndf*nch*ndf + npe * npe * ncu * ncu + npe * npe * ncu * ncu + 
                max(npe * npe, ncu * ncu) + max(npe * ndf * ncu, ncu * ncu * npf) + max(npe,ncu);
    
	if (app.flag_q == 1) {
        elemtempLen += npe*ncu*npe*ncu + npe*ncu*ncu*ndf + ncu*ndf*npe*ncu + 
            ncu*ndf*ncu*ndf + npe*ncq + npe*npe + npe*npe*nd + npe*ndf*nd;
    }
    
//     // To account for the new fields in temp structure:
//     elemtempLen += ncd*max(max(nqvR,nqvJ),nqvQ) + nc*max(max(nqvR,nqvJ),nqvQ) + nc*max(max(nqvR,nqvJ),nqvQ) + 3*max(max(nqvR,nqvJ),nqvQ) + 3*nfe*npf + nfe*nqfR*nch + 
//                 nfe*nqfR*nch*nc + nfe*nqfR*nch*nch + nfe*nqfR*nch + nfe*nqfR*nc + nfe*nqfR*ncd + nfe*nqfR*nd + nfe*nqfR + nfe*nqfR*nd + nfe*nqfR*nch + nfe*nqfR*nc + 
//                 3*nfe*nqfR + nfe*nqfJ*nch + nfe*nqfJ*nch*nc + nfe*nqfJ*nch*nch + nfe*nqfJ*nch + nfe*nqfJ*nc + nfe*nqfJ*ncd + nfe*nqfJ*nd + nfe*nqfJ + 
//                 nfe*nqfJ*nd + nfe*nqfJ*nch + nfe*nqfJ*nc + 3*nfe*nqfJ + nfe*nqfR*nch + nfe*nqfJ*nch*nc + nfe*nqfJ*nch*nch + nqfR*nch + nqfR*nch*nc + nqfR*nch*nch + nqfR*nc + 
//                 nqfR*nch + nqfR*nc + nqfR*nch + nqfR*ncd + nfe*nqfR*nd*nd + nfe*nqfR*nd*nd + nqfR*nd + 3*nqfR + nqfJ*nch + nqfJ*nch*nc + nqfJ*nch*nch + nqfJ*nc + nqfJ*nch + 
//                 nqfJ*nc + nqfJ*nch + nqfJ*ncd + nfe*nqfJ*nd*nd + nfe*nqfJ*nd*nd + nqfJ*nd + 3*nqfJ + nqfJ*nfe*npe*nd + nfe*nqfJ*ncu*npe*ncu + nqfJ*nfe*ndf*nd + 
//                 ndf*ncu*ncu*nfe*nqfJ + nqvJ*npe*nd + nqvJ*ndf*nd + ncu*nqvJ*(nd+1)*ncu*npe + nqvJ*(nd+1)*ncu*ncu*ndf + nqvQ*ncd + nqvQ*nd*nd + nqvQ + nqvQ*nd*nd + 
//                 nfe*nqfQ*ncd + nfe*nqfQ*nd*nd + nfe*nqfQ*nd*nd + nfe*nqfQ + nfe*nqfQ*nd + nfe*nqfQ*nd + nqvJ*ncu +
//                 nqvJ*ncu*nd + nqvJ*ncu*nc + nqvJ*ncu*nd*nc + nqvJ + nqvJ*nd*nd + nqvJ*ncd + nqvJ*nd*nd + nqvR*ncu + nqvR*ncu*nd + nqvR*ncu*nc + 
//                 nqvR*ncu*nd*nc + nqvR + nqvR*nd*nd + nqvR*ncd + nqvR*nd*nd;
    
    /* Initialize sys struct */
    
    // Note: The following sizes are more than required for BJ preconditioner
    sys.oent2ent.resize(cbsr_nblks);            // Only required for implementation No. 3 of BILU0
    sys.oent2entStart.resize(cbsr_nrows+1);     // Only required for implementation No. 3 of BILU0
    sys.LUoent2ent.resize(cbsr_nblks);          // Only required for implementation No. 1 of BILU0
    sys.Loent2entStart.resize(cbsr_nrows+1);    // Only required for implementation No. 1 of BILU0
    sys.Uoent2entStart.resize(cbsr_nrows+1);    // Only required for implementation No. 1 of BILU0

    if (nproc>1) {
        sys.buffsend.resize(max(nentsend*blkSize,nelemsend*nc*npe));
        sys.buffrecv.resize(max(nentrecv*blkSize,nelemrecv*nc*npe));
        if (sys.preconditioner == 0) {
            sys.buffsendmat.resize(blkSize*blkSize*sys.nmatsend);
            sys.buffrecvmat.resize(blkSize*blkSize*sys.nmatrecv);
        }
    }

#ifdef  HAVE_MPI
     if (nproc>1) {
        sys.requests = (MPI_Request *) malloc( 2*sys.nnbsd * sizeof(MPI_Request) );
        sys.statuses = (MPI_Status *) malloc( 2*sys.nnbsd * sizeof(MPI_Status) );
     }
#endif
    
    // Allocate memory for Jacobian matrix:
    sys.Hg.resize(blkSize*blkSize*cbsr_nblks);
    sys.Mg.resize(blkSize*blkSize*cbsr_nblks);        
    
//     // Allocate memory for preconditioner:
//     if (sys.preconditioner == 0)
//         sys.Mg.resize(blkSize*blkSize*cbsr_nblks);
//     else if (sys.preconditioner == 1 || sys.preconditioner == 2) {      
//         if (nproc == 1)
//             sys.Mg.resize(blkSize*blkSize*cbsr_nblks);
//         else if (nproc > 1) {
//             sys.Mg.resize(blkSize*blkSize*BJ_nblks);
//         }
//     }
    
    // Allocate memory for RHS:
    sys.Rg.resize(blkSize*cbsr_nrows);          // This version of RHS is modified during GMRES
    sys.Rg_0.resize(blkSize*cbsr_nrows);        // This version of RHS is NOT modified during GMRES
    if (sys.matvecImplementation == 2) {
        sys.Rg4FD.resize(blkSize*cbsr_nrows);
        sys.r4FD.resize(blkSize*cbsr_nrows);
    }

//     if (nproc > 1) {
//         sys.ent2entWeight.resize( BJ_nblks + BK_nblks );
//     }
    sys.ent2entWeight.resize(cbsr_nblks);

    sys.Mx.resize(blkSize*cbsr_nrows);
    sys.Mv = &sys.Mx[0];                    // Shares memory with sys.Mx
    
    sys.rev.resize(sys.maxiter+1);
    sys.x.resize(blkSize*cbsr_nrows);
    sys.xDenseRow.resize(sys.noThreads);
    if (sys.matvecImplementation == 1) {
        for (i = 0; i < sys.noThreads; i++)
            sys.xDenseRow[i].resize(blkSize*maxBlocksPerRow);
    }
    //sys.xDense.resize(blkSize*cbsr_nblks);
    sys.r.resize(blkSize*cbsr_nrows);
    sys.rtmp.resize(blkSize);
    sys.rtmps.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.rtmps[i].resize(sys.rtmp.size());
    if (nproc > 1)
        sys.v.resize(max(blkSize*BJ_nrows*sys.restart + max(blkSize*cbsr_nrows,blkSize*BJ_nrows),2*sys.noThreads*elemtempLen));           // Shares memory with elem and temp structures
    else
        sys.v.resize(max(blkSize*cbsr_nrows*(sys.restart+1),2*sys.noThreads*elemtempLen));                          // Shares memory with elem and temp structures
    if (sys.orthogMethod == 1 || sys.orthogMethod == 2 || sys.orthogMethod == 3) {
        sys.s.resize(sys.restart+1);        // +1 is required since we compute qNorm by appending q to Q
        if (nproc > 1)
            sys.stmp.resize(sys.s.size());
        sys.stmps.resize(sys.noThreads);
        for (i = 0; i < sys.noThreads; i++)
            sys.stmps[i].resize(sys.s.size());
    }
    if (sys.orthogPrecision == 0) {
        if (nproc > 1) {
            error("This is wrong.\n");
            sys.v_sp.resize(blkSize*BJ_nrows*sys.restart + max(blkSize*cbsr_nrows,blkSize*BJ_nrows));
        }
        else {
            error("This is wrong\n");
            sys.v_sp.resize(blkSize*cbsr_nrows*(sys.restart+1));
        }
        if (sys.orthogMethod == 1 || sys.orthogMethod == 2 || sys.orthogMethod == 3) {
            sys.s_sp.resize(sys.restart+1);        // +1 is required since we compute qNorm by appending q to Q
            if (nproc > 1)
                sys.stmp_sp.resize(sys.s_sp.size());
            sys.stmps_sp.resize(sys.noThreads);
            for (i = 0; i < sys.noThreads; i++)
                sys.stmps_sp[i].resize(sys.s_sp.size());
        }
    }
    sys.y.resize(sys.restart);
    sys.hy.resize(sys.restart+1);
    sys.e1.resize(sys.restart+1);
    for(i=0; i<sys.restart+1; i++)
        sys.e1[i] = 0.0;
    sys.hh.resize((sys.restart+1)*sys.restart);
    sys.hhls.resize((sys.restart+1)*sys.restart);
    sys.work.resize(blkSize);
    sys.works.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.works[i].resize(sys.work.size());
    sys.ipiv.resize(blkSize);
    sys.ipivs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.ipivs[i].resize(sys.ipiv.size());
    sys.ordered2unordered.resize(cbsr_nrows);
    sys.unordered2ordered.resize(cbsr_nrows);
    sys.workls.resize(2*(sys.restart+1));
    for (i=0; i<cbsr_nrows; i++) {
        sys.unordered2ordered[i] = i;
        sys.ordered2unordered[i] = i;
    }
    
    // Vectors for mixed-precision algorithm arrays:
    if (sys.precPrecision == 0) {
        if (sys.preconditioner == 0 || sys.preconditioner == 1 || sys.preconditioner == 2)
            sys.r_sp.resize(cbsr_nrows*blkSize);                    // This is more memory than necessary for preconditioners 1 and 2
        if (sys.preconditioner == 0 || sys.preconditioner == 1)
            sys.rDenseRow_sp.resize(blkSize*maxBlocksPerRow);       // Only used in applyBILU0_sp
        sys.Mg_sp.resize(sys.Mg.size());
        sys.rtmp_sp.resize(blkSize);
        sys.rtmps_sp.resize(sys.noThreads);
        for (i = 0; i < sys.noThreads; i++)
            sys.rtmps_sp[i].resize(sys.rtmp_sp.size());
    }
    if (sys.matvecPrecision == 0) {
        sys.Hg_sp.resize(sys.Hg.size());
        if (sys.matvecImplementation == 1) {
            sys.xDenseRow_sp.resize(sys.noThreads);
            for (i = 0; i < sys.noThreads; i++)
                sys.xDenseRow_sp[i].resize(sys.xDenseRow[i].size());
        }
        sys.x_sp.resize(blkSize*cbsr_nrows);
        if (nproc > 1)
            sys.v_sp.resize(blkSize*BJ_nrows);
        else
            sys.v_sp.resize(blkSize*cbsr_nrows);
        //if ((nproc > 1) && (sys.preconditioner == -1 || sys.preconditioner == 1 || sys.preconditioner == 2))
        sys.Kg_sp.resize(sys.Kg.size());
    }
    
    // Vectors for approximate MDF ordering:
    sys.HiiInv.resize(blkSize*blkSize);
    sys.HiiInvs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.HiiInvs[i].resize(sys.HiiInv.size());
    sys.HiiInvHij.resize(blkSize*blkSize);
    sys.HiiInvHijs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.HiiInvHijs[i].resize(sys.HiiInvHij.size());
    
    // Vectors for exact MDF ordering:
    sys.HkkInv.resize(blkSize*blkSize);
    sys.HkkInvs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.HkkInvs[i].resize(sys.HkkInv.size());
    sys.HikHkkInv.resize(blkSize*blkSize);
    sys.HikHkkInvs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.HikHkkInvs[i].resize(sys.HikHkkInv.size());
    sys.HikHkkInvHkj.resize(blkSize*blkSize);
    sys.HikHkkInvHkjs.resize(sys.noThreads);
    for (i = 0; i < sys.noThreads; i++)
        sys.HikHkkInvHkjs[i].resize(sys.HikHkkInvHkj.size());
    
    Int start = 0;
    sys.C = &sys.v[start];
    start += cbsr_nblks;
    sys.w = &sys.v[start];
    start += cbsr_nrows;
    
    /* Initialize temp struct */    
    temps.pg = &sys.v[start];
    start += nge*ncd;
    temps.Jg = &sys.v[start];
    start += nge*nd*nd;
    temps.jacg = &sys.v[start];
    start += nge;
    temps.Xxg = &sys.v[start];
    start += nge*nd*nd;
    temps.pf = &sys.v[start];
    start += npf*nfe*ncd;
    temps.pgf = &sys.v[start];
    start += ngf*nfe*ncd;
    temps.Jgf = &sys.v[start];
    start += ngf*nfe*nd*(nd-1);
    temps.nlgf = &sys.v[start];
    start += ngf*nfe*nd;
    temps.jacgf = &sys.v[start];
    start += ngf*nfe;
    temps.nlgjac = &sys.v[start];
    start += ngf*nfe;
    temps.udgg = &sys.v[start];
    start += nge*nc;
    temps.odgg = &sys.v[start];
    start += nge*nco;
//     temps.udgg_ref.resize(nge*nc);
    temps.f = &sys.v[start];
    start += nge*ncu*nd;
    temps.f_udg = &sys.v[start];
    start += nge*ncu*nd*nc;
    temps.s = &sys.v[start];
    start += nge*ncu;
    temps.s_udg = &sys.v[start];
    start += nge*ncu*nc;
    temps.wrl = &sys.v[start];
    start += nge*ncu*(nd+1)*nc;
    temps.wrk = &temps.wrl[0];            // Shares memory with temp.wrl
    temps.uh = &sys.v[start];
    start += ndf*nch;
    //temps.uh_ref = &sys.v[start];
    //start += ndf*nch;
//     temps.uh_ref.resize(ndf*nch);
    temps.uf = &sys.v[start];
    start += ndf*nc;
    temps.of = &sys.v[start];
    start += ndf*nco;
//     temps.uf_ref.resize(ndf*nc);
    temps.ugf = &sys.v[start];
    start += ngf*nfe*nc;
    temps.ogf = &sys.v[start];
    start += ngf*nfe*nco;
//     temps.ugf_ref.resize(ngf*nfe*nc);
    temps.uhg = &sys.v[start];
    start += ngf*nfe*nch;
    //temps.uhrefg = &sys.v[start];
    //start += ngf*nfe*nch;
//     temps.uhrefg.resize(ngf*nfe*nch);
    //temps.avg_p1CG = &sys.v[start];
    //start += 3*nge;
//     temps.avg_p1CG.resize(3*nge);
    //temps.avf_p1CG = &sys.v[start];
    //start += 3*npf*nfe;
//     temps.avf_p1CG.resize(3*npf*nfe);
    //temps.avfg_p1CG = &sys.v[start];
    //start += 3*ngf*nfe;
//     temps.avfg_p1CG.resize(3*ngf*nfe);
    temps.fh = &sys.v[start];
    start += ngf*nfe*nch;
    temps.fh_u = &sys.v[start];
    start += ngf*nfe*nch*nc;
    temps.fh_uh = &sys.v[start];
    start += ngf*nfe*nch*nch;
    temps.Rutmp = &sys.v[start];
    start += ndf*nch;
    temps.BDtmp = &sys.v[start];
    start += npf*ndf*nch*nc;
    temps.BDt = &sys.v[start];
    start += npe*npe*nch*nc;
    temps.Ft = &sys.v[start];
    start += npe*ncu*ndf*nch;
    temps.Ftmp = &sys.v[start];
    start += npf*npf*nfe*nch*nch;
    temps.Etmp = &sys.v[start];
    start += npf*npf*nfe;
    temps.fb = &sys.v[start];
    start += ngf*nch;
    temps.fb_u = &sys.v[start];
    start += ngf*nch*nc;
    temps.fb_uh = &sys.v[start];
    start += ngf*nch*nch;
    temps.pft = &sys.v[start];
    start += ngf*ncd;
    temps.uft = &sys.v[start];
    start += ngf*nc;
    temps.uht = &sys.v[start];
    start += ngf*nch;
    temps.nlt = &sys.v[start];
    start += ngf*nd;
    //temps.uft_ref = &sys.v[start];
    //start += ngf*nc;
//     temps.uft_ref.resize(ngf*nc);
    //temps.uht_ref = &sys.v[start];
    //start += ngf*nch;
//     temps.uht_ref.resize(ngf*nch);
    //temps.avft_p1CG = &sys.v[start];
    //start += 3*ngf;
//     temps.avft_p1CG.resize(3*ngf);

    temps.ipiv.resize(npe*ncu);
    temps.ind.resize(npe);
    temps.jnd.resize(npe);

    if (app.flag_q == 1) {
        temps.BMiC = &sys.v[start];
        start += npe*ncu*npe*ncu;
        temps.BMiE = &sys.v[start];
        start += npe*ncu*ncu*ndf;
        temps.GMiC = &sys.v[start];
        start += ncu*ndf*npe*ncu;
        temps.GMiE = &sys.v[start];
        start += ncu*ndf*ncu*ndf;
    }

    temps.shapMiC.resize(nge*npe*nd); 
    temps.shapMiE.resize(nge*ndf*nd); 
    temps.wrlshapMiC.resize(nge*(nd+1)*ncu*npe*ncu);        
    temps.wrlshapMiE.resize(nge*(nd+1)*ncu*ncu*ndf);

    temps.MiCf.resize(npf*nfe*npe*nd); 
    temps.shapMiCf.resize(ngf*nfe*npe*nd); 
    temps.fhushapMiCf.resize(ngf*nfe*npe*ncu*ncu);
    temps.Df.resize(npf*nfe*npe*ncu*ncu);    
    temps.MiEf.resize(npf*nfe*ndf*nd); 
    temps.shapMiEf.resize(ngf*nfe*ndf*nd); 
    temps.fhushapMiEf.resize(ngf*nfe*ncu*ncu*ndf);
    temps.Ff.resize(npf*nfe*ncu*ncu*ndf);

    if (app.hybrid == 1 || app.hybrid == 2) {       // For EDG, IEDG and HEDG only
        temps.K_tmp.resize(ncu*npe*ncf*nch);
        temps.F_tmp.resize(ncf*nch*ncu*npe);
        temps.H_tmp.resize(ncf*nch*ncf*nch);
        temps.Rh_tmp.resize(ncf*nch);
    }

    // Initialize elem struct
    /*  [M -C  E] [dq] = [Rq]
        [B  D  F] [du] = [Ru]
        [G  K  H] [dh] = [Rh] */
    if (app.flag_q == 1) {
        elems.Rq = &sys.v[start];
        start += npe*ncq;
        elems.M = &sys.v[start];
        start += npe*npe;
        elems.C = &sys.v[start];
        start += npe*npe*nd;
        elems.E = &sys.v[start];
        start += npe*ndf*nd;
    }
    elems.Ru = &sys.v[start];
    start += npe*ncu;
    elems.BD = &sys.v[start];
    start += npe*npe*ncu*nc;
    elems.F = &sys.v[start];
    start += npe*ncu*ndf*nch;
    elems.Rh = &sys.v[start];
    start += nch*ndf;
    elems.Rhonly = &sys.v[start];
    start += nch*ndf;
    elems.GK = &sys.v[start];
    start += nch*ndf*npe*nc;
    elems.H = &sys.v[start];
    start += nch*ndf*nch*ndf;
    elems.D_tmp = &sys.v[start];
    start += npe * npe * ncu * ncu;
    elems.D_LU = &sys.v[start];
    start += npe * npe * ncu * ncu;
    elems.D_LU_tmp = &sys.v[start];
    start += max(npe * npe, ncu * ncu);
    elems.DiF_ij = &sys.v[start];
    start += max(npe * ndf * ncu, ncu * ncu * npf);
    elems.Ru_i = &sys.v[start];
    start += max(npe,ncu);
    elems.pivDii.resize(max(npe,ncu));
    elems.workDii.resize(max(npe,ncu));
    elems.workD.resize(npe*ncu);
    elems.F_dense.resize(ncu*npf*ncu*npf*nfe);
    elems.F_conv.resize(npe*nch*ndf*ncu);
    elems.D_inv.resize(ncu*npe*ncu*npe);
    elems.D_inv_extCol.resize(ncu*npe*ncu*ndf);
    elems.Di_tmp.resize(ncu*npe*ncu*npe);
    elems.DiF_conv.resize(ncu*npe*nch*ndf);
    elems.DiF_extRow.resize(ncu*ndf*nch*ndf);
    elems.DiF_tmp.resize(ncu*npe*nch*ndf);
    elems.K_tmp.resize(nch*ndf*npe*ncu);
    elems.K.resize(nch*ndf*npe*ncu);
    for(j=0; j<nch*npf*nfe*npe*nc; j++)
        elems.GK[j] = 0.0;
    for(j=0; j<npe*ncu*nch*npf*nfe; j++)
        elems.F[j] = 0.0;
        
//         temps.pp = &sys.v[start];
//         start += ncd*max(max(nqvR,nqvJ),nqvQ);
//         temps.udgp = &sys.v[start];
//         start += nc*max(max(nqvR,nqvJ),nqvQ);
//         temps.udgp_ref = &sys.v[start];
//         start += nc*max(max(nqvR,nqvJ),nqvQ);
//         temps.avp = &sys.v[start];
//         start += 3*max(max(nqvR,nqvJ),nqvQ);
//         
//         temps.avf = &sys.v[start];
//         start += 3*nfe*npf;

//         temps.fhR = &sys.v[start];
//         start += nfe*nqfR*nch;
//         temps.fhR_u = &sys.v[start];     // NOT NECESSARY
//         start += nfe*nqfR*nch*nc;           // NOT NECESSARY
//         temps.fhR_uh = &sys.v[start];    // NOT NECESSARY
//         start += nfe*nqfR*nch*nch;          // NOT NECESSARY
//         temps.uhR = &sys.v[start];
//         start += nfe*nqfR*nch;
//         temps.ufR = &sys.v[start];
//         start += nfe*nqfR*nc;
//         temps.pfR = &sys.v[start];
//         start += nfe*nqfR*ncd;
//         temps.nlR = &sys.v[start];
//         start += nfe*nqfR*nd;
//         temps.jacfR = &sys.v[start];
//         start += nfe*nqfR;
//         temps.nljacR = &sys.v[start];        // NOT NECESSARY
//         start += nfe*nqfR*nd;
//         temps.uh_refR = &sys.v[start];
//         start += nfe*nqfR*nch;
//         temps.uf_refR = &sys.v[start];
//         start += nfe*nqfR*nc;
//         temps.avfR = &sys.v[start];
//         start += 3*nfe*nqfR;
// 
//         temps.fhJ = &sys.v[start];
//         start += nfe*nqfJ*nch;
//         temps.fhJ_u = &sys.v[start];
//         start += nfe*nqfJ*nch*nc;
//         temps.fhJ_uh = &sys.v[start];
//         start += nfe*nqfJ*nch*nch;
//         temps.uhJ = &sys.v[start];
//         start += nfe*nqfJ*nch;
//         temps.ufJ = &sys.v[start];
//         start += nfe*nqfJ*nc;
//         temps.pfJ = &sys.v[start];
//         start += nfe*nqfJ*ncd;
//         temps.nlJ = &sys.v[start];
//         start += nfe*nqfJ*nd;
//         temps.jacfJ = &sys.v[start];
//         start += nfe*nqfJ;
//         temps.nljacJ = &sys.v[start];
//         start += nfe*nqfJ*nd;
//         temps.uh_refJ = &sys.v[start];
//         start += nfe*nqfJ*nch;
//         temps.uf_refJ = &sys.v[start];
//         start += nfe*nqfJ*nc;
//         temps.avfJ = &sys.v[start];
//         start += 3*nfe*nqfJ;
//         
//         temps.fhjacR = &sys.v[start];
//         start += nfe*nqfR*nch;
//         temps.fh_ujacJ = &sys.v[start];
//         start += nfe*nqfJ*nch*nc;
//         temps.fh_uhjacJ = &sys.v[start];
//         start += nfe*nqfJ*nch*nch;
// 
//         temps.fbR = &sys.v[start];
//         start += nqfR*nch;
//         temps.fbR_u = &sys.v[start];       // NOT NECESSARY
//         start += nqfR*nch*nc;                 // NOT NECESSARY
//         temps.fbR_uh = &sys.v[start];      // NOT NECESSARY
//         start += nqfR*nch*nch;                // NOT NECESSARY
//         temps.uftR = &sys.v[start];
//         start += nqfR*nc;
//         temps.uhtR = &sys.v[start];
//         start += nqfR*nch;
//         temps.uft_refR = &sys.v[start];
//         start += nqfR*nc;
//         temps.uht_refR = &sys.v[start];
//         start += nqfR*nch;
//         temps.pftR = &sys.v[start];
//         start += nqfR*ncd;
//         temps.JfR = &sys.v[start];
//         start += nfe*nqfR*nd*nd;
//         temps.XxfR = &sys.v[start];
//         start += nfe*nqfR*nd*nd;
//         temps.nltR = &sys.v[start];
//         start += nqfR*nd;
//         temps.avftR = &sys.v[start];
//         start += 3*nqfR;
// 
//         temps.fbJ = &sys.v[start];
//         start += nqfJ*nch;
//         temps.fbJ_u = &sys.v[start];
//         start += nqfJ*nch*nc;
//         temps.fbJ_uh = &sys.v[start];
//         start += nqfJ*nch*nch;
//         temps.uftJ = &sys.v[start];
//         start += nqfJ*nc;
//         temps.uhtJ = &sys.v[start];
//         start += nqfJ*nch;
//         temps.uft_refJ = &sys.v[start];
//         start += nqfJ*nc;
//         temps.uht_refJ = &sys.v[start];
//         start += nqfJ*nch;
//         temps.pftJ = &sys.v[start];
//         start += nqfJ*ncd;
//         temps.JfJ = &sys.v[start];
//         start += nfe*nqfJ*nd*nd;
//         temps.XxfJ = &sys.v[start];
//         start += nfe*nqfJ*nd*nd;
//         temps.nltJ = &sys.v[start];
//         start += nqfJ*nd;
//         temps.avftJ = &sys.v[start];
//         start += 3*nqfJ;
//         
//         temps.shapMiCfJ = &sys.v[start];
//         start += nqfJ*nfe*npe*nd;
//         temps.fh_ushapMiCfJ = &sys.v[start];
//         start += nfe*nqfJ*ncu*npe*ncu;
//         temps.shapMiEfJ = &sys.v[start];
//         start += nqfJ*nfe*ndf*nd;
//         temps.fh_ushapMiEfJ = &sys.v[start];
//         start += ndf*ncu*ncu*nfe*nqfJ;
//         
//         temps.shapMiC_J = &sys.v[start];
//         start += nqvJ*npe*nd;
//         temps.shapMiE_J = &sys.v[start];
//         start += nqvJ*ndf*nd;
//         temps.wrlshapMiC_J = &sys.v[start];
//         start += ncu*nqvJ*(nd+1)*ncu*npe;
//         temps.wrlshapMiE_J = &sys.v[start];
//         start += nqvJ*(nd+1)*ncu*ncu*ndf;
//         
//         temps.pvQ = &sys.v[start];
//         start += nqvQ*ncd;
//         temps.Jv_Q = &sys.v[start];
//         start += nqvQ*nd*nd;
//         temps.jacvQ = &sys.v[start];
//         start += nqvQ;
//         temps.XxvQ = &sys.v[start];
//         start += nqvQ*nd*nd;
// 
//         temps.pfQ = &sys.v[start];
//         start += nfe*nqfQ*ncd;
//         temps.Jf_Q = &sys.v[start];
//         start += nfe*nqfQ*nd*nd;
//         temps.XxfQ = &sys.v[start];
//         start += nfe*nqfQ*nd*nd;
//         temps.jacfQ = &sys.v[start];
//         start += nfe*nqfQ;
//         temps.nlfQ = &sys.v[start];
//         start += nfe*nqfQ*nd;
//         temps.nlfjacQ = &sys.v[start];
//         start += nfe*nqfQ*nd;
//         
//         temps.sJ = &sys.v[start];
//         start += nqvJ*ncu;
//         temps.fJ = &sys.v[start];
//         start += nqvJ*ncu*nd;
//         temps.sJ_udg = &sys.v[start];
//         start += nqvJ*ncu*nc;
//         temps.fJ_udg = &sys.v[start];
//         start += nqvJ*ncu*nd*nc;
//         temps.jacJ = &sys.v[start];
//         start += nqvJ;
//         temps.XxJ = &sys.v[start];
//         start += nqvJ*nd*nd;
//         temps.pJ = &sys.v[start];
//         start += nqvJ*ncd;
//         temps.J_J = &sys.v[start];
//         start += nqvJ*nd*nd;
//         
//         temps.sR = &sys.v[start];
//         start += nqvR*ncu;
//         temps.fR = &sys.v[start];
//         start += nqvR*ncu*nd;
//         temps.sR_udg = &sys.v[start];      // NOT NECESSARY
//         start += nqvR*ncu*nc;
//         temps.fR_udg = &sys.v[start];      // NOT NECESSARY
//         start += nqvR*ncu*nd*nc;
//         temps.jacR = &sys.v[start];
//         start += nqvR;
//         temps.XxR = &sys.v[start];
//         start += nqvR*nd*nd;
//         temps.pR = &sys.v[start];
//         start += nqvR*ncd;
//         temps.J_R = &sys.v[start];
//         start += nqvR*nd*nd;
    
    
// // // //     // Compute element measure:
// // // //     wrfwerv; //THIS DOES NOT WORK WITH OPEN MP... THINK HOW TO FIX IT (MOST LIKELY CHANGING computeElementMeasure)
// // // //     computeElementMeasure(sys, elem, mesh, master, sol, app, temp, ndims);
}

#endif
