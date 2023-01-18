#include "../getan.c"
#include "flux_ns2d.c"
#include "flux_ns3d.c"

void flux_ns(double *f, double *f_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            flux_ns2d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            flux_ns3d(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void fluxonly_ns(double *f, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fluxonly_ns2d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fluxonly_ns3d(f, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "source_ns2d.c"
#include "source_ns3d.c"

void source_ns(double *s, double *s_udg, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            source_ns2d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            source_ns3d(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void sourceonly_ns(double *s, double *pg, double *udg, double *param, 
        double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            sourceonly_ns2d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            sourceonly_ns3d(s, pg, udg, param, time, ng, nc, ncu, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fhat_ns2d.c"
#include "fhat_ns3d.c"

void fhat_ns(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhat_ns2d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhat_ns3d(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

void fhatonly_ns(double *fh, double *pg, double *udg, 
        double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fhatonly_ns2d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        case 3:            
            fhatonly_ns3d(fh, pg, udg, uh, nl, param, time, ng, nc, ncu, nd, ncd);        
            break;
        default:
            exit(-1);
            break;
    }
}

#include "fbou_ns2d.c"
#include "fbou_ns3d.c"

void fbou_ns(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbou_ns2d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbou_ns3d(fh, fh_u, fh_uh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}

void fbouonly_ns(double *fh, double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nch, int nc, int nd, int ncd)
{
    switch (nd) {
        case 2:    
            fbouonly_ns2d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;        
        case 3:            
            fbouonly_ns3d(fh, pg, udg, uhg, nl,
                      ui, param, time, ib, ng, nch, nc, nd, ncd);
            break;
        default:
            exit(-1);
            break;
    }
}


