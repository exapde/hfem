#ifndef __DIGASOPRE_H
#define __DIGASOPRE_H

#include "wrapper.h"

//appstruct
struct appstruct {
    Int my_rank;
    Int nproc;
    Int nfile;
    Int nd;
    Int ncd;
    Int nc;
    Int ncu;
    Int nch;
    Int ncq;
    Int ncp;
    Int nco;    
    vector<Int> ndims;
    vector<Int> meshsize;
    vector<Int> solsize;
    vector<Int> porder;
    vector<Int> elemtype;
    vector<Int> nodetype;
    vector<Int> pgauss;
    vector<Int> pgaussR;
    vector<Int> quadtype;
    vector<Int> bcm;
    vector<double> bcs;
    vector<Int> bcd;
    vector<double> bcv;
    vector<double> dt;      
    vector<Int> flag;
    vector<double> factor;
    vector<Int> problem;
    vector<double> physicsparam;
    vector<double> solversparam;
    
    string filein;          // Name of binary file with input data
    string fileout;         // Name of binary file to write the solution
    
    double time;            /* current time */
    double fc_u;            /* factor when discretizing the time derivative of the U equation. Allow scalar field for local time stepping in steady problems? */
    double fc_q;            /* factor when discretizing the time derivative of the Q equation. Allow scalar field for local time stepping in steady problems? */
    double fc_p;            /* factor when discretizing the time derivative of the P equation. Allow scalar field for local time stepping in steady problems? */
    double dtfc;            /* needed only for augmented Lagrangian */
    double alpha;           /* needed only for augmented Lagrangian */

    Int tdep;               /* flag to determine unsteady or steady problems */
    Int wave;               /* flag to determine wave or non-wave problems */
    Int alag;               /* flag to determine augmented lagrangian */
    Int ALEflag;            // Flag for Arbitrary Lagrangian-Eulerian (ALE) formulation. 0: No ALE; 1: Translation; 2: Translation + rotation; 3: Translation + rotation + deformation
    Int adjoint;            /* flag to determine adjoint or primal problems */
    Int linearProblem;      /* flag to determine linear or nonlinear problems */
    Int flag_q;             /* flag to determine if the Q equation is also discretized  */
    Int flag_p;             /* flag to determine if the P equation is also discretized  */
    Int flag_g;             /* flag to determine if the GCL equation is also discretized  */
    Int debugmode;
    
    Int overlappinglevel;
    Int quasiNewton;        // flag to determine quasi Newton or full Newton
    Int reuseJacobian;      // flag to determine whether or not the Jacobian matrix needs to be computed
    Int reuseResidual;      // flag to determine whether or not the residual vector needs to be computed
    Int reusePreconditioner;// flag to determine whether or not the preconditioner needs to be computed
    Int reuseOrdering;      // flag to determine whether or not the DOF ordering needs to be computed

    Int temporalscheme;
    Int hybrid;             // Discretization scheme. 0: HDG; 1: EDG; 2: IEDG
    Int appname;            // 0: Compressible Euler; 1: Compressible Navier-Stokes; 2: Compressible regularized RANS-SA
    Int linearSolver;       // 0: direct solver; 1: gmres; etc. (future: CG, CGS, QMR)
    Int dirkStage;          // DIRK stages
    Int dirkOrder;          // DIRK order
    Int BDFsteps;           // Number of steps for BDF scheme
    Int jacobianStep;       // reuse Jacobian for every jacobianStep time steps
    Int orderingStep;       // reuse Ordering for every orderingStep time steps

    vector<Int> DIRKnotConverged;
    vector<Int> BDFnotConverged;
    
    Int AVflag;             // Flag for artificial viscosity. 0: No artificial viscosity; 1: Homogeneous artificial viscosity (C. Nguyen's formulation); 2: Hypersonic homogeneous artificial viscosity (C. Nguyen's formulation)
                                    // 3: Isotropic artificial viscosity (D. Moro's formulation). 4: Latest version of the model (taking the best of all previous models)
                                    // 8: Density smoothness sensor (Per's approach)
    double rampFactor;      // Ramp factor for artificial viscosity flux
    Int viscosityModel;     // Flag for viscosity model. 0: Constant dynamic viscosity; 1: Sutherland's law

    Int SGSmodel;           // Flag for sub-grid scale (SGS) model. 0: No SGS model. 1: Smagorinsky/Yoshizawa/Knight model. 
                            //                                      2: WALE/Yoshizawa/Knight model. 3: Vreman/Yoshizawa/Knight model. SGS model only available for 3D solver.   
    Int convStabMethod;     // Flag for convective stabilization tensor. 0: Constant tau, 1: Lax-Friedrichs; 2: Roe.
    Int diffStabMethod;     // Flag for diffusive stabilization tensor. 0: No diffusive stabilization.
    Int rotatingFrame;      // Flag for rotating frame. Options: 0: Velocities are respect to a non-rotating frame. 1: Velocities are respect to a rotating frame.
};


// element info struct
struct elementinfo {    
    Int dim;      // spatial dimension
    Int elemtype; // type     
    
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    Int nqf; // number of quad faces
    Int ntf; // number of tet faces
    vector< Int > facetype; // 0 triangular, 1 quadrilateral 
    vector< vector< Int > > face; // face connectivity     
    vector< vector< Int > > edge; // edge connectivity     
};

struct masterstruct {    
    Int dim;      // spatial dimension
    Int nd;      // spatial dimension
    Int elemtype; // % element type  
    Int nodetype; // node type  
    Int npe; // number of DG nodes for solution
    Int nme; // number of DG nodes for geometry
    Int npv; // number of DG nodes for solution
    Int nmv; // number of DG nodes for geometry
    Int nge; // number of Gauss nodes
    Int ngv; // number of Gauss nodes
    Int ngeR; // number of Gauss nodes
    Int ngvR; // number of Gauss nodes
    Int nfe; // number of faces    
    Int nle;  // number of edges 
    Int nve; // number of vertices 
    Int nplmax; // number number of DG nodes on each element vertice
    Int ndf;
    Int naf;
    Int nqf;
    Int nqfR;
    
    vector<Int> porder; // polynomial degrees along each direction
    vector<Int> pgauss;    
    vector<Int> pgaussR;    
    vector<Int> ndims;
    vector<Int> npf;
    vector<Int> nmf;
    vector<Int> ngf;
    vector<Int> ngfR;
    vector<Int> npfix;
    vector<Int> nmfix;
    vector<Int> ngfix;
    vector<Int> ngfRix;    
    vector<Int> nvf;
    vector<Int> npl;
    
    vector<Int> permnode;
    vector<Int> permedge;
    vector< vector<Int> > perm;
    vector< vector<Int> > permgeom;
    vector< vector<Int> > face;
    
    vector<double> plocvl;
    vector<Int> tlocvl;
    
    vector<double> gpvl;
    vector<double> gwvl;    
    vector<double> shapvl;
    vector<double> shapvt;
    vector<double> shapvg;    
    vector<double> shapvgdotshapvl;
    vector<double> gpvlR;
    vector<double> gwvlR;
    vector<double> shapvtR;
    vector<double> shapvgR;
    vector<double> shapmv;
    vector<double> shapmvR;
    
    vector< vector<double> > plocfc;
    vector< vector< Int > > tlocfc;  
    
    vector< vector<double> > gpfc;
    vector< vector<double> > gwfc;    
    vector< vector<double> > shapfc;
    vector< vector<double> > shapft;
    vector< vector<double> > shapfg;
    vector< vector<double> > shapfgdotshapfc;
    vector< vector<double> > gpfcR;
    vector< vector<double> > gwfcR;        
    vector< vector<double> > shapftR;
    vector< vector<double> > shapfgR;
    vector< vector<double> > shapmf;
    vector< vector<double> > shapmfR;
    
    vector<double> shapnv;
    vector<double> shapnvt;   
    vector<double> projLowP;
        
    // NEW FIELDS FOR NEW MATRIX ASSEMBLY:
    //Int pgaussR;
    Int pgaussJ;
    Int pgaussQ;
    
    Int quadTypeR;
    Int quadTypeJ;
    Int quadTypeQ;
    
    Int nqvR;
    Int nqvQ;
    Int nqvJ;
    vector<double> gpvlJ;
    vector<double> gpvlQ;    
    vector<double> gwvlJ;
    vector<double> gwvlQ;    
    vector<double> shapvlJ;
    vector<double> shapvlQ;    
    vector<double> shapvtJ;
    vector<double> shapvtQ;    
    vector<double> shapvgJ;
    vector<double> shapvgQ;
    vector<double> shapvgdotshapvlR;
    vector<double> shapvgdotshapvlJ;
    vector<double> shapvgdotshapvlQ;
        
    vector<Int> nqfQ;
    vector<Int> nqfJ;    
    vector< vector<double> > gpfcJ;
    vector< vector<double> > gpfcQ;    
    vector< vector<double> > gwfcJ;
    vector< vector<double> > gwfcQ;
    vector< vector<double> > shapfcJ;
    vector< vector<double> > shapfcQ;
    vector< vector<double> > shapftJ;
    vector< vector<double> > shapftQ;
    vector< vector<double> > shapfgJ;
    vector< vector<double> > shapfgQ;
    vector< vector<double> > shapfgdotshapfcR;
    vector< vector<double> > shapfgdotshapfcJ;
    vector< vector<double> > shapfgdotshapfcQ;
};
       
struct meshstruct {
    Int my_rank;
    Int dim;   // spatial dimension    
    Int ne; // number of elements
    Int nf; // number of faces
    Int nv; // number of vertices  
    Int ndofuh;
    Int nfemax; // maximum number of faces on elements    
    Int nlemax;  // maximum number of edges on elements
    Int nvemax; // maximum number of vertices on elements    
    Int npemax; // maximum number of nodes on elements
    Int nmemax; // maximum number of nodes on elements
    Int ngemax; // maximum number of nodes on elements
    Int ngeRmax; // maximum number of nodes on elements    
    Int nvfmax; // maximum number of vertices on faces
    Int npfmax; // maximum number of nodes on faces  
    Int ngfmax;
    Int ngfRmax;
    Int ndfmax;
    Int ncfmax;
    Int nmfmax; // maximum number of nodes on faces    
    Int nodetype; // 0 uniform, 1 optimal        
    Int hybrid;
    
    vector<Int> ndf;    
    vector<Int> ncf;    
    vector<Int> nfes;
    vector<Int> npes;
    vector< vector<Int> > npfixs;        
    vector< vector<Int> > npfs;        
    vector< vector<Int> > nmfs;        
    vector< vector<Int> > perm;        
    
    vector<Int> ndims;            
    vector<double> p; // coordinates of vertice
    vector<double> dgnodes; // coordinates of DG nodes   
    vector<Int> t;  // element-to-vertice connectivities
    vector<Int> dgnodesidx; // element index of dgnodes        
    vector<Int> elcon; // element-to-entity connectivities
    vector<Int> elconidx; // element index of elcon        
    //vector< vector<Int> > npfidx;  // face index of npf
    vector<Int> bf;  // numbering of boundary faces
    vector<Int> t2f; // element-to-face connectivities
    vector<Int> t2t; // element-to-element connectivities
    vector<Int> f;   // face-to-element connnectivities 
    vector<Int> extintelem;
    vector<Int> extintface;
    vector<Int> extintent;
    vector<Int> elem2cpu;
    vector<Int> face2cpu;
    vector<Int> ent2cpu;
    vector<Int> nbsd;
    vector< vector<Int> > cg2dg;
    vector< vector<Int> > dg2cg;             
            
    vector<Int> porder;    
    vector<Int> pgauss;
    vector<Int> quadtype;    
    //vector<elementstruct> elements; // types of elements in the mesh    
    //vector<masterstruct> masters; // types of masters in the mesh    
    vector<Int> elementtype;  // associated with the structure elements 
    vector<Int> physics;     // associated with the governing equations            
    vector<Int> isEDGface;    
    vector<Int> isEDGelement;    
    
    vector<double> elemMeasure; // Element measure: ne
    vector<double> hAvg;        // Characteristic element size: ne
    vector<double> M;           // Metric tensor at dgnodes
    vector<double> Minv;        // Inverse of metric tensor at dgnodes
};

struct dmdstruct {    
    Int my_rank; // index of the processor
    vector<Int> nbsd;    // neighboring cpus
    vector<Int> intelem; // nonoverlapping global elements
    vector<Int> intent;  // nonoverlapping global entities
    vector<Int> elempart; // overlapping global elements
    vector<Int> elempartpts; // classifiers of overlapping global elements: (interior, interface, exterior)
    vector<Int> entpart;  // overlapping global entities
    vector<Int> entpartpts; // classifiers of overlapping global entities: (interior, interface, exterior)
    vector<Int> elemrecv;  // local elements received from neighboring cpus
    vector<Int> elemrecvpts; // classifiers of local elements received from neighboring cpus: (# elements from ncpu1, ...) 
    vector<Int> elemsend;  // local elements sent to neighboring cpus
    vector<Int> elemsendpts; // classifiers of local elements sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> entrecv;    // local entities received from neighboring cpus
    vector<Int> entrecvpts; // classifiers of local entities received from neighboring cpus: (# entities from ncpu1, ...)          
    vector<Int> entsend;  // local entities sent to neighboring cpus
    vector<Int> entsendpts; // classifiers of local entities sent to neighboring cpus: (# elements to ncpu1, ...)    
    vector<Int> vecrecv;  // local vectors received from neighboring cpus to perform matrix-vector product
    vector<Int> vecrecvpts; // classifiers of local vectors received from neighboring cpus: (# vectors from ncpu1, ...)          
    vector<Int> vecsend;    // local vectors sent to neighboring cpus to perform matrix-vector product
    vector<Int> vecsendpts; // classifiers of local vectors sent to neighboring cpus: (# vectors to ncpu1, ...)             
    vector<Int> matrecv;    // local matrices received from neighboring cpus to construct the preconditioner
    vector<Int> matrecvpts; // classifiers of local matrices received from neighboring cpus: (# matrices from ncpu1, ...)                
    vector<Int> matsend;   // local matrices sent to neighboring cpus to construct the preconditioner
    vector<Int> matsendpts;  // classifiers of local matrices sent to neighboring cpus: (# matrices to ncpu1, ...)                                         
    vector<Int> rowent2elem; // global entity-to-element connectivities
    vector<Int> colent2elem; // global entity-to-element connectivities
    vector<Int> rowent2ent;  // global entity-to-entity connectivities
    vector<Int> colent2ent;  // global entity-to-entity connectivities
    vector<Int> bcrs_rowent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_colent2elem; // local entity-to-element connectivities
    vector<Int> bcrs_rowent2ent;  // local entity-to-entity connectivities
    vector<Int> bcrs_colent2ent;  // local entity-to-entity connectivities    
    vector<Int> elemmap;  // element reordering 
    vector<Int> entmap;  // entity reordering 
    vector<Int> ent2ind;  // global-to-local entity mapping  
    vector<Int> elcon;    // local element-to-entity connectivities
    vector<Int> t2f;      // local elemeent-to-face connectivities
    //vector<Int> t2t;      // local elemeent-to-element  connectivities
    Int  maxBlocksPerRow; // maximum number of entities per row
    Int  minBlocksPerRow; // minimum number of entities per row    
};

struct solstruct {
    vector<Int> ndims;
    vector<double> UDG; /* UDG = (U, Q, P) */
    vector<double> UH;  /* UH = UHAT */
    vector<double> PDG; /* Additional unknown for wave applications */
    vector<double> SH;  /* Source from the discretization of the time derivatives */
    vector<double> SP;  /* Additional source for wave applications */
    vector<double> ODG; /* other DG fields */
    
    vector<double> UDG4FD; /* Only used for matrix-free matrix-vector product */
    vector<double> UH4FD;  /* Only used for matrix-free matrix-vector product */
    
    vector<double> UDG_avg; // Time average UDG
    vector<double> UH_avg;  // Time average UH
    
    vector<double> Un; /* Un = (U, UHAT) */
    vector<double> Um; /* Um = (U, UHAT) */
    vector<double> Vn; /* Vn = (U, UHAT) */
    vector<double> R;  /* R = (RU, RUHAT) */
    
    vector<double> avField_DG;  // DG artificial viscosity field
    vector<double> avField_p1CG;  // p=1 CG artificial viscosity field
    vector<double> avField;
    
    vector<double> UDG2Write; // Auxiliary variable to write solution to a binary file. Only required if writeQflag = 1.
    
    // For Minimal Residual algorithm to compute the initial guess for Newton iteration
    vector<double> UDG_initMR;
    vector<double> UH_initMR;
    
    /* for time-dependent problems only */
    vector<double> UDGi; /* store UDG from the previous timestep */
    vector<double> UHi;  /* store UH from the previous timestep */
    vector<double> PDGi; /* store PDG from the previous timestep */
    
    /* used to recover DU */
    vector<double> DinvRu;  /* inverse(D)*Ru */
    vector<double> DinvF;   /* inverse(D)*F */

    /* Only needed for Quasi Newton. TODO: Wouldn't it be better to store F instead of DinvF for quasi-Newton? */
    vector<double> Dinv;    /* inv( dRu_bar/du ) */
    vector<float> DinvFloat;    /* inv( dRu_bar/du ) */
    vector<Int> ipivD;      /* pivots for LU factors of dRu_bar/du */
    vector<double> K;       /* dRh/du */
    vector<float> Kfloat;       /* dRh/du */

    double alpha;  /* stepsize for damped Newton */
    double rNorm;  /* the norm of the residual vector for both U and UH */
    
    vector<double> Cs;      // Constant for dynamic Smagorinsky SGS model: npv / ne
    
    /* call to write (UDG,UH) into a binary file */
    void writeSol2File(string filename, Int ne, Int n, Int* ndims, Int writeQflag);
    
    /* call to write (UDG_avg,UH_avg) into a binary file */
    void writeAvgSol2File(string filename_avg, Int ne, Int n, Int* ndims, Int writeQflag);
};

struct elemstruct {
/* storing the following system for one element
 [M -C  E] [dq] = [Rq]
 [B  D  F] [du] = [Ru]
 [G  K  H] [dh] = [Rh]
*/
    double*  M;
    double*  E;
    double*  C;
    double*  Rq;
    double*  BD;
    double*  F;
    double*  Ru;
    double*  GK;
    double*  H;
    double*  Rh;
    double*  Rhonly;
    vector<double> D1;
    vector<double> F1;
    vector<double> K1;
    vector<double> H1;
    vector<double> K;

    Int my_rank;
    
    vector<double> F_dense;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> F_conv;              ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_tmp;                     ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> K_tmp;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_LU;                      ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  D_LU_tmp;                  ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> D_inv;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> D_inv_extCol;        ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_conv;            ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_tmp;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> DiF_extRow;          ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  DiF_ij;                    ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> Di_tmp;              ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    double*  Ru_i;                      ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<Int> pivDii;                 ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> workDii;             ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
    vector<double> workD;               ////////// ONLY REQUIRED FOR ALTERNATIVE SCHUR IMPLEMENTATIONS
};


struct tempstruct {
    double* pg;
    double* Jg;
    double* jacg;
    double* Xxg;

    double* pgR;
    double* JgR;
    double* jacgR;
    double* XxgR;
    
    double* pf;
    double* pgf;
    double* Jgf;
    double* nlgf;
    double* jacgf;

    double* pfR;
    double* pgfR;
    double* JgfR;
    double* nlgfR;
    double* jacgfR;
    
    double* udgg;
    double* odgg;
    double* f;
    double* f_udg;
    double* s;
    double* s_udg;
    
    double* udggR;
    double* odggR;
    double* sR;
    double* fR;
    double* sR_udg;      // NOT NECESSARY
    double* fR_udg;      // NOT NECESSARY

    double* uh;    
    double* uf;
    double* of;
    double* uhg;
    double* ugf;
    double* ogf;
    double* fh;
    double* fh_u;
    double* fh_uh;

    double* uhgR;
    double* ugfR;
    double* ogfR;
    double* fhR;
    double* fhR_u;
    double* fhR_uh;    
    
    double* pft;
    double* uft;
    double* oft;
    double* uht;
    double* nlt;
    
    //here
    
    double* nlgjac;    
    double* udgg_ref;
    double* wrk;
    double* wrl;
    
    double* uh_ref;    
    double* uf_ref;    
    double* ugf_ref;    
    double* uhrefg;
    double* avg_p1CG;
    double* avf_p1CG;
    double* avfg_p1CG;
    double* Rutmp;
    double* BDtmp;
    double* BDt;    
    double* Ft;
    double* Ftmp;
    double* Etmp;
    double* fb;
    double* fb_u;
    double* fb_uh;    
    double* uft_ref;
    double* uht_ref;
    double* avft_p1CG;

    double* BMiE;
    double* BMiC;
    double* GMiE;
    double* GMiC;        
    
    vector<double> shapMiC; //(ngv*npv*nd); 
    vector<double> shapMiE; //(ngv*ndf*nd); 
    vector<double> wrlshapMiC; //(ngv*nd1*ncu*npv*ncu);        
    vector<double> wrlshapMiE; //(ngv*nd1*ncu*ncu*ndf);
    
    vector<double> MiCf; 
    vector<double> shapMiCf; 
    vector<double> fhushapMiCf;
    vector<double> Df;
    
    vector<double> MiEf; 
    vector<double> shapMiEf; 
    vector<double> fhushapMiEf;
    vector<double> Ff;
    
    vector<Int> ipiv;
    vector<Int> ind;
    vector<Int> jnd;

    vector<double> H_tmp;
    vector<double> K_tmp;
    vector<double> F_tmp;
    vector<double> Rh_tmp;
    
    // FIELDS FOR NEW MATRIX ASSEMBLY WITH DIFFERENT QUADRATURE RULES FOR JACOBIAN AND RESIDUAL:
    double* pp;
    double* udgp;
    double* udgp_ref;
    double* avp;
    
    double* avf;
    
    double* uhR;
    double* ufR;    
    double* nlR;
    double* jacfR;
    double* nljacR;//NOT NECESSARY
    double* uh_refR;
    double* uf_refR;
    double* avfR;
    
    double* fhJ;
    double* fhJ_u;
    double* fhJ_uh;
    double* uhJ;
    double* ufJ;
    double* pfJ;
    double* nlJ;
    double* jacfJ;
    double* nljacJ;
    double* uh_refJ;
    double* uf_refJ;
    double* avfJ;
    
    double* fhjacR;
    double* fh_ujacJ;
    double* fh_uhjacJ;
    
    double* fbR;
    double* fbR_u;
    double* fbR_uh;
    double* uftR;
    double* uhtR;
    double* uft_refR;
    double* uht_refR;
    double* pftR;
    double* JfR;
    double* XxfR;
    double* nltR;
    double* avftR;
    
    double* fbJ;
    double* fbJ_u;
    double* fbJ_uh;
    double* uftJ;
    double* uhtJ;
    double* uft_refJ;
    double* uht_refJ;
    double* pftJ;
    double* JfJ;
    double* XxfJ;
    double* nltJ;
    double* avftJ;
    
    double* shapMiCfJ;
    double* fh_ushapMiCfJ;
    double* shapMiEfJ;
    double* fh_ushapMiEfJ;
    
    double* pvQ;
    double* Jv_Q;
    double* jacvQ;
    double* XxvQ;
    
    double* pfQ;
    double* Jf_Q;
    double* XxfQ;
    double* jacfQ;
    double* nlfQ;
    double* nlfjacQ;
    
    double* sJ;
    double* fJ;
    double* sJ_udg;
    double* fJ_udg;
    double* jacJ;
    double* XxJ;
    double* pJ;
    double* J_J;
    
//     double* sR;
//     double* fR;
//     double* sR_udg;      // NOT NECESSARY
//     double* fR_udg;      // NOT NECESSARY
//     double* jacR;
//     double* XxR;
//     double* pR;
//     double* J_R;
    
    double* shapMiC_J;
    double* shapMiE_J;
    double* wrlshapMiC_J;
    double* wrlshapMiE_J;
};

struct sysstruct {
    /* TODO: Create an structure for the preconditioner */
    Int numEntities;        // Number of entities in the mesh (nfn in EDG or nf in HDG)
    Int numBlocks;          // Number of blocks in matrix of linear system
    Int blkSize;            // Block size (nch in EDG or nch*npf in HDG)
    Int maxBlocksPerRow;    // Maximum number of blocks per row in the matrix of the linear system
    Int BJ_nrows;           // Number of rows associated to entities in the processor
//     Int BJ_nblks;           // Number of blocks in rows associated to entities in the processor
//     Int BK_nrows;           // Number of rows associated to entities in the processor that contain columns not in the processor
//     Int BK_nblks;           // Number of blocks in rows associated to entities in the processor that contain columns not in the processor
    Int nproc;              // Number of processors
    Int nentpartpts;        
    Int nentrecv;
    Int nentsend;
    Int nelempartpts;        
    Int nelemrecv;
    Int nelemsend;
    Int nvecrecv;
    Int nvecsend;    
    Int nmatrecv;
    Int nmatsend;    
    Int nnbsd;
    Int my_rank;            // Index of the processor
    Int noThreads;         // Number of OpenMP threads
    Int computeGlobalEnt2entWeight;       // Only applies for MPI code. If set to 1, the Frobenius norm of the Schur matrix blocks are computed, stored in a file and the execution is terminated.
    
    vector<Int> ent2ent;
    vector<Int> blockStartJ;
    vector<Int> ent2entStart;
//    vector<Int> globalEnt2ent;      // Only used in (1) MPI code and (2) if sys.computeGlobalEnt2entWeight == 1
//    vector<Int> globalEnt2entStart; // Only used in (1) MPI code and (2) if sys.computeGlobalEnt2entWeight == 1
    vector<Int> oent2ent;           // Only required for implementation No. 2 of BILU0
    vector<Int> oent2entStart;      // Only required for implementation No. 2 of BILU0
    vector<Int> LUoent2ent;         // Only required for implementation No. 3 of BILU0
    vector<Int> Loent2entStart;     // Only required for implementation No. 3 of BILU0
    vector<Int> Uoent2entStart;     // Only required for implementation No. 3 of BILU0
//     vector<Int> BJ_rowpts;
//     vector<Int> BJ_colind;
//     vector<Int> BK_rowpts;
//     vector<Int> BK_colind;
//     vector<Int> BK_rowind;
    vector<double> ent2entWeight;       // Only from entities allocated to the processor (weights from entities in other processors are treated by those other processors)
//     vector<Int> ent2entWeightLen;       // Length of the vector to be received by processor No. 1 from each processor
    vector<Int> entpart;                // local-to-global entity mapping
    vector<Int> entpartpts;
    vector<Int> entrecv;
    vector<Int> entrecvpts;
    vector<Int> entsend;
    vector<Int> entsendpts;
    vector<Int> elempart;
    vector<Int> elempartpts;
    vector<Int> elemrecv;
    vector<Int> elemrecvpts;
    vector<Int> elemsend;
    vector<Int> elemsendpts;
    vector<Int> vecrecv;
    vector<Int> vecrecvpts;
    vector<Int> vecsend;
    vector<Int> vecsendpts;    
    vector<Int> matrecv;
    vector<Int> matrecvpts;
    vector<Int> matsend;
    vector<Int> matsendpts;
    vector<Int> nbsd;
    
// // // //     Int numRowsP;
// // // //     Int numBlocksP;
// // // //     Int numRowsJ;               // <- numEntities
// // // //     Int numBlocksJ;             // <- numBlocks
// // // //     Int nIntRows;               // <- BJ_nrows
// // // //     vector<Int> ent2entStartJ;  // <- ent2entStart
// // // //     vector<Int> ent2entJ;       // <- ent2ent
// // // //     vector<Int> ent2entStartP;
// // // //     vector<Int> ent2entP;
    
    // Variables for Minimal Residual algorithm for non-linear initial guess:
    vector<double> Ru_MR;
    vector<double> dRuda_MR;
    vector<double> a2_MR;
    vector<vector<double> > a2s_MR;
    vector<double> da_MR;
    vector<double> C_MR;
    vector<double> Clocal_MR;
    
    Int orthogMethod;               // Orthogonalization method in GMRES. 0: MGS, 1: ICGS, 2: IMGS (only MGS and ICGS are available in MPI code), 3: ICGS without convergence guarantee
    Int maxiter;                    // Maximum number of GMRES iterations
    Int restart;                    // Parameter k in GMRES(k)
    double tol;                     // Convergence tolerance for linear system
    Int preconditionerSide;         // 0: Left preconditioner. 1: Right preconditioner
    Int reorderMethod;              // 0: No reordering. 1: Approximate (inexact) MDF. 2: Exact MDF. 3: Approximate (inexact) MDF with constraints (only valid for parallel BJ preconditioner). 4: Exact MDF with constraints (only valid for parallel BJ preconditioner)
    Int preconditioner;             // -1: No preconditioner, 0: Restricted Additive Schwarz (RAS), 1: Subdomain-wise Block Jacobi (BJ), 2: Entity-wise block Jacobi
    Int NewtonMaxiter;              // Maximum number of Newton/quasi-Newton iterations
    Int trueNewtonMaxiter;          // Maximum number of true Newton iterations
    double NewtonTol;               // Newton tolerance for non-linear system

    double rNorm;                   /* the norm of the residual vector for both U and UH */

    Int print;                      // 0: no print, 1: print on screen
    Int schurImplementation;        // Implementation flag. 0 and 1 values are valid
    Int matvecImplementation;       // Implementation flag for matrix-vector product. Options: 0, 1, 2 (finite differences for Gateaux derivative)
    Int precSolveImplementation;    // Implementation flag for preconditioner solve. Options: 0 and 1
    Int precPrecision;              // Precision for preconditioner solve. Options: 0: Single precision. 1: Double precision
    Int matvecPrecision;            // Precision for matrix-vector product. Options: 0: Single precision. 1: Double precision
    Int orthogPrecision;            // Precision for orthogonalization. Options: 0: Single precision. 1: Double precision
    Int quasiNewtonAccuracy;        // Format in which Dinv and K are stored for quasi-Newton. 0: Single precision; 1: Double precision
    Int adaptiveGMREStol;           // Adaptive GMRES tolerance for nonlinear problems. Options: 1: Adaptive method. 0: Non-adaptive method. Note the adaptive strategy is not used for linear PDEs regardless the value of this flag.
    Int adaptGMREStol;              // Only applies for linear problems with adaptiveGMREStol == 1. Has two meanings:
                                    // Before executing GMRES routine: If adapting GMRES tolerance is allowed.
                                    // After executing GMRES routine: If GMRES tolerance was adapted in the previous linear solve
    Int robustMode;                 // Options: 0 (default): All flags to default value. 1: All flags to safety (or robust) mode
    
    Int linearProblem;      /* flag to determine linear or nonlinear problems */
    long linearSolvesClocks;         // Clocks spent in lear solves since the Jacobian matrix was computed the last time (only applies for quasi-Newton)
    long lastAssemblyAndPrecClocks;  // Clocks spent to perform the last matrix assembly and preconditioner computation (only applies for quasi-Newton)

    vector<double> Hg;              // Part of global matrix such that M_i ~= (Hg_i)^-1 (BILU0 inverse)
    vector<double> Kg;              // Part of global matrix not contained in Hg
    vector<double> Rg;              // Global RHS (this vector will be modified during GMRES)
    vector<double> Rg4FD;           // Global RHS (used for matrix-free matrix-vector product)
    vector<double> Rg_0;            // Global RHS (this vector will not be modified during GMRES)
    vector<double> x;               // Solution to linear system (Hg * x = Rg)
    vector<vector<double> > xDenseRow;       // Dense version of x for a particular row
//    vector<double> xDense;          // Dense version of x
    vector<double> Mg;              // Preconditioner matrix
    vector<double> r;               // Residual vector
    vector<double> r4FD;            // Residual vector (used for matrix-free matrix-vector product)
    vector<double> v;               // Krylov vectors (also have other uses)

    vector<double> s;               // Auxiliary variable for orthogonalization
    vector<double> stmp;            // Auxiliary variable for orthogonalization
    vector<vector<double> > stmps;            // Auxiliary variable for orthogonalization
    vector<double> Mx;              // Preconditioner matrix apply to the solution vector
    double* Mv;                     // Preconditioner matrix apply to the last vector in the Arnoldi basis

    // Single precision arrays for mixed-precision algorithms:
    vector<float> Hg_sp;
    vector<float> Kg_sp;
    vector<float> Mg_sp;
    vector<float> x_sp;
    vector<vector<float> > xDenseRow_sp;       // Dense version of x_sp for a particular row
    vector<float> v_sp;            // Krylov vectors in single precision
    vector<float> s_sp;               // Auxiliary variable for orthogonalization
    vector<float> stmp_sp;            // Auxiliary variable for orthogonalization
    vector<vector<float> > stmps_sp;            // Auxiliary variable for orthogonalization
    vector<float> r_sp;
    vector<float> rtmp_sp;
    vector<vector<float> > rtmps_sp;
    vector<float> rDenseRow_sp;

    vector<Int> ipiv;                   // Pivots for LU factors of diagonal blocks of BILU0 preconditioner
    vector<vector<Int> > ipivs;
    vector<Int> ordered2unordered;      // Mapping from ordered entities to unordered entities
    vector<Int> unordered2ordered;      // Mapping from unordered entities to ordered entities

    double* C;
    double* w;

    vector<double> e1;
    vector<double> rev;
    vector<double> y;
    vector<double> hy;
    vector<double> hh;
    vector<double> hhls;
    vector<double> rtmp;
    vector<vector<double> > rtmps;
    vector<double> work;
    vector<vector<double> > works;
    vector<double> workls;

    vector<double> buffsend;
    vector<double> buffrecv;
    vector<double> buffrecvEnt2entWeight;
    vector<double> buffsendmat;
    vector<double> buffrecvmat;
    
    // Auxiliary vectors for approximate MDF ordering:
    vector<double> HiiInv;
    vector<vector<double> > HiiInvs;
    vector<double> HiiInvHij;
    vector<vector<double> > HiiInvHijs;
    
    // Auxiliary vectors for exact MDF ordering:
    vector<double> HkkInv;
    vector<vector<double> > HkkInvs;
    vector<double> HikHkkInv;
    vector<vector<double> > HikHkkInvs;
    vector<double> HikHkkInvHkj;
    vector<vector<double> > HikHkkInvHkjs;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;

//    MPI_Request * requestsEnt2entWeight;
//    MPI_Status * statusesEnt2entWeight;
#endif
};

#endif
