#include "flux_appname2d.c"
#include "flux_appname3d.c"

// Written by: C. Nguyen & P. Fernandez

void flux_appname(double *f, double *f_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            flux_appname2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            flux_appname3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void fluxonly_appname(double *f, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fluxonly_appname2d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fluxonly_appname3d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


#include "source_appname2d.c"
#include "source_appname3d.c"

void source_appname(double *s, double *s_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            source_appname2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            source_appname3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void sourceonly_appname(double *s, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            sourceonly_appname2d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            sourceonly_appname3d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fhat_appname2d.c"
#include "fhat_appname3d.c"

void fhat_appname(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhat_appname2d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhat_appname3d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

void fhatonly_appname(double *fh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhatonly_appname2d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhatonly_appname3d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fbou_appname2d.c"
#include "fbou_appname3d.c"

void fbou_appname(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbou_appname2d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbou_appname3d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


void fbouonly_appname(double *fh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbouonly_appname2d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbouonly_appname3d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

