void fbou_poisson3d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{            
    
    int    i;                   
    double kappa   = param[0];
    double tau = param[1];   
    
    if (ib==1) { /* Dirichlet */
        for (i=0; i<ng*nch; i++) 
            fh[i] = tau*(ui[0]-uhg[i]);

        for (i=0; i<ng*nch*nc; i++)             
            fh_u[i] = 0.0;
        
        for (i=0; i<ng*nch*nch; i++)             
            fh_uh[i] = -tau;                
    }                                     
    else if (ib==2) { /* Neumman */
        fhat_poisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = fh[i] + ui[0];
    }                                     
    else if (ib==3) { /* Neumman */
        fhat_poisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = fh[i] + nl[i];
    }                                     
    else if (ib==4) { /* Neumman */
        fhat_poisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = fh[i] + nl[ng+i];
    }             
    else if (ib==5) { /* Neumman */
        fhat_poisson3d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = fh[i] + nl[2*ng+i];
    }                                     
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}


// void fbouonly_poisson3d(double *fh, double *pg, double *udg, double *uhg, double *nl,
//           double *ui, double *param, double time, int ib,
//           int ng, int nc, int nch, int nd, int ncd)
// {        
//     
//     int    i;                   
//     double kappa   = param[0];
//     double tau = param[1];    
//     
//     if (ib==1) { /* Dirichlet */
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = tau*(ui[0]-uhg[i]);           
//     }                                     
//     else if (ib==2) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + ui[0];
//     }                                     
//     else if (ib==3) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[i];
//     }                                     
//     else if (ib==4) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[ng+i];
//     }                                     
//     else if (ib==5) { /* Neumman */
//         fhatonly_poisson(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd); 
//         
//         for (i=0; i<ng*nch; i++) 
//             fh[i] = fh[i] + nl[2*ng+i];
//     }                                     
//     else {                        
//         printf("This error is in %s on line %d\n",__FILE__, __LINE__);
//         printf("Boundary condition %d is not implemented yet.", ib);            
//         exit(-1);                                    
//     }                
// }
// 
// 
