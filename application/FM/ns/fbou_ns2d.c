double finduplus(double yp, double u0)
{
    double kappa = 0.41, B = 5.0; 
    double C = exp(-kappa*B);
    double uk = kappa*u0;
    double y0 = u0 + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);

    double u1, ua, up, y1, ya;
    double tol = 1e-6;
    
//     cout<<"finduplus"<<endl;
//     cout<<yp<<"  "<<u0<<endl;
//     cout<<y0<<"  "<<uk<<endl;
    
    if (abs(yp-y0)<tol)
        up = u0;
    else {
        // find lower (u0) and upper (u1) bounds
        if (y0<yp)
        {
            u1 = (2.0*u0 > 1.0) ? 2*u0 : 1.0;
            uk = kappa*u1;
            y1 = u1 + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0); 
            while (y1<yp)
            {
                u1 = 2*u1;
                uk = kappa*u1;
                y1 = u1 + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);
            }
        }
        else
        {
            u1 = u0;
            y1 = y0;
            u0 = u1/2;
            uk = kappa*u0;
            y0 = u0 + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);            
            while (y0>yp)
            {                
                u0 = u0/2;
                uk = kappa*u0;
                y0 = u0 + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);
            }
        }

//         cout<<"here"<<endl;
//         cout<<y0<<"  "<<yp<<"  "<<y1<<endl;
//         cout<<u0<<"  "<<u1<<endl;
        
        // bisection
        ua = 0.5*(u0+u1);        
        uk = kappa*ua;
        ya = ua + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);
        while (abs(yp-ya)>tol)
        {            
            if (ya<yp)
                u0 = ua;
            else
                u1 = ua;
            ua = 0.5*(u0+u1);
            uk = kappa*ua;
            ya = ua + C*(exp(uk) - 1.0 - uk - 0.5*uk*uk - uk*uk*uk/6.0);
        }
        up = ua;
    }        
    
    return up;
}

void fbou_ns2d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{        
    
    double gam, epslm, gam1, Minf, M2, Re;    
    double nx, ny, tm;        
    int    i, j, k, m, n, nm, nk, nn;
                   
    gam   = param[0];
    epslm = param[1];
    Re    = param[2];
    Minf  = param[4];
    gam1  = gam-1.0;
    M2    = Minf*Minf;
    
    if (ib==1) { /* freestream / far-field */

        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *anm  = new double [ng*nch*nch*nch];
        double *Anm  = new double [ng*nch*nch*nch];                    
    
        getan(an,anm,uhg,nl,param,0,ng,nch,nd);
        getan(An,Anm,uhg,nl,param,1,ng,nch,nd);                                
                
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;

        for (i=0; i<ng*nch*nc; i++)             
            fh_u[i] = 0.0;
        
        for (i=0; i<ng*nch*nch; i++)             
            fh_uh[i] = 0.0;
        
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(ui[k] - uhg[nk]);

                    fh_u[nm] = an[nm] + An[nm];      
                    fh_uh[nm] = -2.0*An[nm];  
                    for (n=0; n<nch; n++) {
                        nk = k*nch*nch*ng+n*nch*ng+j*ng+i;
                        nn = n*ng+i;
                        fh_uh[nm] = fh_uh[nm] + (anm[nk] + Anm[nk])*(udg[nn] - uhg[nn]) - 
                                 (anm[nk] - Anm[nk])*(ui[n] - uhg[nn]);  
                    }
                }
                
        delete[] an; delete[] An; delete[] anm; delete[] Anm;     
        
    }                                     
    else if (ib==2) { /* adiabatic wall */

        fhat_ns(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);                          

        for (i=0; i<ng; i++) {

            fh[0*ng+i] =  udg[0*ng+i]-uhg[0*ng+i];
            fh[1*ng+i] = -uhg[1*ng+i];
            fh[2*ng+i] = -uhg[2*ng+i];

            fh_uh[0*ng+i] = -1.0;
            fh_uh[1*ng+i] = 0.0;
            fh_uh[2*ng+i] = 0.0; 
            fh_uh[4*ng+i] = 0.0; 
            fh_uh[5*ng+i] = -1.0;
            fh_uh[6*ng+i] = 0.0;
            fh_uh[8*ng+i] = 0.0;
            fh_uh[9*ng+i] = 0.0;
            fh_uh[10*ng+i] = -1.0;
            fh_uh[12*ng+i] = 0.0;
            fh_uh[13*ng+i] = 0.0;
            fh_uh[14*ng+i] = 0.0;

            for (j = 0; j <nch-1; j++)
                for (k = 0; k <nc; k++)
                    fh_u[k*nch*ng+j*ng+i] = 0.0;
            fh_u[0*ng+i] = 1.0;                                    
        }            
    }
    else if (ib==3) { /* isothermal wall */

        fhat_ns(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);

        for (i=0; i<ng; i++) { 
            fh[1*ng+i] = uhg[1*ng+i];
            fh[2*ng+i] = uhg[2*ng+i];
            fh[3*ng+i] = gam*gam1*M2*uhg[3*ng+i]/uhg[0*ng+i] - ui[3];

            fh_uh[1*ng+i] = 0.0;
            fh_uh[2*ng+i] = 0.0;
            fh_uh[3*ng+i] = -gam*gam1*M2*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
            fh_uh[5*ng+i] = 1.0;
            fh_uh[6*ng+i] = 0.0;
            fh_uh[7*ng+i] = 0.0;
            fh_uh[9*ng+i] = 0.0;
            fh_uh[10*ng+i] = 1.0;
            fh_uh[11*ng+i] = 0.0;
            fh_uh[13*ng+i] = 0.0;
            fh_uh[14*ng+i] = 0.0;
            fh_uh[15*ng+i] = gam*gam1*M2/uhg[0*ng+i];

            for (j = 0; j <nch-1; j++)
                for (k = 0; k <nc; k++)
                    fh_u[k*nch*ng+j*ng+i] = 0.0;
        }                         
    }
    else if (ib==4) { /* pressure far-field */
        
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *anm  = new double [ng*nch*nch*nch];
        double *Anm  = new double [ng*nch*nch*nch];
        double *uinf  = new double [ng*nch];
        double *uinfu  = new double [ng*nch*nch]; 
        
        double pinf = (ui[3]-0.5)*(gam-1);
        
        getan(an,anm,uhg,nl,param,0,ng,nch,nd);
        getan(An,Anm,uhg,nl,param,1,ng,nch,nd);
        
        for (i=0; i<ng*(nch-1); i++) 
            uinf[i] = udg[i];        
        for (i=0; i<ng; i++) {
            n = (nch-1)*ng+i;                            
            uinf[n] = pinf/(gam-1) + 0.5*udg[ng+i]*udg[ng+i]/udg[i] + 0.5*udg[2*ng+i]*udg[2*ng+i]/udg[i];
        }        
        
        int sz2 = ng * nch;
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++)
                for (i=0; i<ng; i++) {
                    nm = k*sz2+j*ng+i;
                    uinfu[nm] = (j==k) ? 1.0 : 0.0;
                }
        
        for (i = 0; i < ng; i++) {
            nm = 0 * sz2 + (1 + nd) * ng + i;
            uinfu[nm] = 0.0;
            for (j = 0; j < nd; j++)
                uinfu[nm] -= 0.5 * udg[(j + 1) * ng + i] * udg[(j + 1) * ng + i] / (udg[0*ng+i] * udg[0*ng+i]);

            for (j = 0; j < nd; j++) {
                nm = (j + 1) * sz2 + (1 + nd) * ng + i;
                uinfu[nm] = udg[(j + 1) * ng + i] / udg[0*ng+i];
            }

            nm = (1 + nd) * sz2 + (1 + nd) * ng + i;
            uinfu[nm] = 0.0;
        }
        
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;

        for (i=0; i<ng*nch*nc; i++)             
            fh_u[i] = 0.0;
        
        for (i=0; i<ng*nch*nch; i++)             
            fh_uh[i] = 0.0;

        int na, nb;
        for (k=0; k<nch; k++)        
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(uinf[nk] - uhg[nk]);

                    fh_u[nm] = an[nm] + An[nm];      
                    fh_uh[nm] = -2.0*An[nm];  
                    for (n=0; n<nch; n++) {
                        nk = k*nch*nch*ng+n*nch*ng+j*ng+i;
                        nn = n*ng+i;
                        na = n*nch*ng+j*ng+i;
                        nb = k*nch*ng+n*ng+i;
                        fh_u[nm]  = fh_u[nm] - (an[na] - An[na])*uinfu[nb];
                        fh_uh[nm] = fh_uh[nm] + (anm[nk] + Anm[nk])*(udg[nn] - uhg[nn]) - 
                                 (anm[nk] - Anm[nk])*(uinf[nn] - uhg[nn]);  
                    }
                }

                        
       delete[] an; delete[] An; delete[] anm; delete[] Anm; delete[] uinf; delete[] uinfu;         
    }                             
    else if (ib==5 | ib == 6) { /* freestream / far-field */

        double *ub  = new double [ng*nch];
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *anm  = new double [ng*nch*nch*nch];
        double *Anm  = new double [ng*nch*nch*nch];                    
    
        getan(an,anm,uhg,nl,param,0,ng,nch,nd);
        getan(An,Anm,uhg,nl,param,1,ng,nch,nd);                                
                
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;

        for (i=0; i<ng*nch*nc; i++)             
            fh_u[i] = 0.0;
        
        for (i=0; i<ng*nch*nch; i++)             
            fh_uh[i] = 0.0;        
        
        //double x = param[6];
        double x = (ib <= 5) ? 1.0 : 2.0;
        double mu = 1.0/Re; 
        double Rex = x*Re;        
        double Cf = 0.027/(pow(Rex,1.0/7.0));
        double tauw = Cf/2.0;
        double utau = sqrt(tauw);        
        double yplus, uplus=1.0;
        for (k=0; k<nch; k++)
            for (i=0; i<ng; i++)
            {                
                if (k == 1) {                             
                    yplus = pg[ng+i]*utau/mu;            
                    uplus = finduplus(yplus, uplus);
                    ub[ng+i] = uplus*utau;
                    if (ub[ng+i] > 1.0)
                       ub[ng+i] = 1.0;     
                    //cout<<ib<<"   "<<i<<"   "<<pg[ng+i]<<"   "<<ub[ng+i]<<endl;                    
                }
                else
                    ub[k*ng+i] = ui[k];
                cout<<ib<<"   "<<k<<"   "<<i<<"   "<<ub[k*ng+i]<<endl;                    
            }
        
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(ub[nk] - uhg[nk]);

                    fh_u[nm] = an[nm] + An[nm];      
                    fh_uh[nm] = -2.0*An[nm];  
                    for (n=0; n<nch; n++) {
                        nk = k*nch*nch*ng+n*nch*ng+j*ng+i;
                        nn = n*ng+i;
                        fh_uh[nm] = fh_uh[nm] + (anm[nk] + Anm[nk])*(udg[nn] - uhg[nn]) - 
                                 (anm[nk] - Anm[nk])*(ub[nk] - uhg[nn]);  
                    }
                }
                
        delete[] ub; delete[] an; delete[] An; delete[] anm; delete[] Anm;     
        
    }                                         
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}


void fbouonly_ns2d(double *fh, double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{        
    
    double gam, epslm, gam1, Minf, Re, M2;    
    double nx, ny, tm;        
    int    i, j, k, m, n, nm, nk, nn;
                   
    gam  = param[0];
    epslm= param[1];         
    Re   = param[2];
    Minf = param[4];    
    gam1 = gam-1.0;    
    M2   = Minf*Minf;            
    
    if (ib==1) { /* freestream */

        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        
        getanonly(an,uhg,nl,param,0,ng,nch,nd);
        getanonly(An,uhg,nl,param,1,ng,nch,nd);                        
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;
        
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(ui[k] - uhg[nk]);
                    }
                
        delete[] an; delete[] An; 
    }                                     
    else if (ib==2) { /* adiabatic */

        fhatonly_ns(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);                          

        for (i=0; i<ng; i++) {
            fh[0*ng+i] =  udg[0*ng+i]-uhg[0*ng+i];
            fh[1*ng+i] = -uhg[1*ng+i];
            fh[2*ng+i] = -uhg[2*ng+i]; 
        }            
    }
    else if (ib==3) { /* isothermal */

        fhatonly_ns(fh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);

        for (i=0; i<ng; i++) { 
            fh[1*ng+i] = uhg[1*ng+i];
            fh[2*ng+i] = uhg[2*ng+i];
            fh[3*ng+i] = gam*gam1*M2*uhg[3*ng+i]/uhg[0*ng+i] - ui[3];
        }                         
    }
    else if (ib==4) { /* pressure */
        
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *uinf  = new double [ng*nch];
        
        double pinf = (ui[3]-0.5)*(gam-1);
        
        getanonly(an,uhg,nl,param,0,ng,nch,nd);
        getanonly(An,uhg,nl,param,1,ng,nch,nd);                        
        
        for (i=0; i<ng*(nch-1); i++) 
            uinf[i] = udg[i];        
        for (i=0; i<ng; i++) {
            n = (nch-1)*ng+i;                            
            uinf[n] = pinf/(gam-1) + 0.5*udg[ng+i]*udg[ng+i]/udg[i] + 0.5*udg[2*ng+i]*udg[2*ng+i]/udg[i];
        }        
                                       
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;

        int na, nb;
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(uinf[nk] - uhg[nk]);
                    }                
                        
       delete[] an; delete[] An; delete[] uinf; 
    }                              
    else if (ib==5 | ib==6) { /* freestream / far-field */

        double *ub  = new double [ng*nch];
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
    
        getanonly(an,uhg,nl,param,0,ng,nch,nd);
        getanonly(An,uhg,nl,param,1,ng,nch,nd);                        
        
        for (i=0; i<ng*nch; i++) 
            fh[i] = 0.0;
                
        
        //double x = param[6];
        double x = (ib <= 5) ? 1.0 : 2.0;
        double mu = 1.0/Re; 
        double Rex = x*Re;        
        double Cf = 0.027/(pow(Rex,1.0/7.0));
        double tauw = Cf/2.0;
        double utau = sqrt(tauw);        
        double yplus, uplus=1.0;
        for (k=0; k<nch; k++)
            for (i=0; i<ng; i++)
            {                
                if (k == 1) {                    
                    yplus = pg[ng+i]*utau/mu;            
                    uplus = finduplus(yplus, uplus);
                    ub[ng+i] = uplus*utau;
                }
                else
                    ub[k*ng+i] = ui[k];
            }
                
        for (k=0; k<nch; k++)
            for (j=0; j<nch; j++) 
                 for (i=0; i<ng; i++) {            
                    nm = k*nch*ng+j*ng+i;
                    nk = k*ng+i;
                    fh[j*ng+i] = fh[j*ng+i] + (an[nm] + An[nm])*(udg[nk] - uhg[nk]) - 
                                 (an[nm] - An[nm])*(ub[nk] - uhg[nk]);
                    }
                
        delete[] ub; delete[] an; delete[] An;         
    }                                             
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}




