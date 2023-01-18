#include "flux_poisson2d.c"
#include "flux_poisson3d.c"

// Written by: C. Nguyen & P. Fernandez

void flux_poisson(double *f, double *f_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            flux_poisson2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            flux_poisson3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fluxonly_poisson2d.c"
#include "fluxonly_poisson3d.c"

void fluxonly_poisson(double *f, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fluxonly_poisson2d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fluxonly_poisson3d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


#include "source_poisson2d.c"
#include "source_poisson3d.c"

void source_poisson(double *s, double *s_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            source_poisson2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            source_poisson3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "sourceonly_poisson2d.c"
#include "sourceonly_poisson3d.c"

void sourceonly_poisson(double *s, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            sourceonly_poisson2d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            sourceonly_poisson3d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fhat_poisson2d.c"
#include "fhat_poisson3d.c"

void fhat_poisson(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhat_poisson2d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhat_poisson3d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fhatonly_poisson2d.c"
#include "fhatonly_poisson3d.c"

void fhatonly_poisson(double *fh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhatonly_poisson2d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhatonly_poisson3d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fbou_poisson2d.c"
#include "fbou_poisson3d.c"

void fbou_poisson(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbou_poisson2d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbou_poisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void fbouonly_poisson(double *fh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbouonly_poisson2d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbouonly_poisson3d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

