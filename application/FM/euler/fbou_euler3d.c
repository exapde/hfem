void fbou_euler3d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{        
    
    double gam, epslm, gam1, Minf, M2, Re;
    double nx, ny, nz, tm;
    int    i, j, k, m, n, nm, nk, nn;
                   
    gam  = param[0];
    epslm= param[1];                     
    Re   = param[2];
    Minf = param[4];    
    gam1 = gam-1.0;    
    M2   = Minf*Minf;   
    
    int sz2 = ng * nch;
    
    if (ib==1) { /* freestream */

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
    else if(ib==2) { /* inviscid wall */

        double ur, uru, urv, urw, urE, hr, hru, hrv, hrw, hrE, unu;
        double uinf1, uinf2, uinf3, uinf4, uinf5;

        for (i=0; i<ng; i++) {

            // normal vector
            nx  = nl[0*ng+i];              
            ny  = nl[1*ng+i]; 
            nz  = nl[2*ng+i]; 
            
            // UDG
            ur  = udg[0*ng+i];
            uru = udg[1*ng+i];
            urv = udg[2*ng+i];
            urw = udg[3*ng+i];
            urE = udg[4*ng+i];
            
            // UH
            hr  = uhg[0*ng+i];
            hru = uhg[1*ng+i];
            hrv = uhg[2*ng+i];
            hrw = uhg[3*ng+i];
            hrE = uhg[4*ng+i];
                              
            // make uinf
            unu   = uru*nx + urv*ny + urw*nz;
            uinf1 = ur;
            uinf2 = uru - unu*nx;
            uinf3 = urv - unu*ny;
            uinf4 = urw - unu*nz;
            uinf5 = urE;            
            
            // fh = uinf - uh;
            fh[0*ng+i] = uinf1 - hr;
            fh[1*ng+i] = uinf2 - hru;    
            fh[2*ng+i] = uinf3 - hrv;
            fh[3*ng+i] = uinf4 - hrw;
            fh[4*ng+i] = uinf5 - hrE;

            // fh_udg = uinfu;
            fh_u[0*ng+i] = 1.0;
            fh_u[1*ng+i] = 0.0;
            fh_u[2*ng+i] = 0.0;
            fh_u[3*ng+i] = 0.0;        
            fh_u[4*ng+i] = 0.0;
            
            fh_u[5*ng+i] = 0.0;
            fh_u[6*ng+i] = 1.0 - nx*nx;
            fh_u[7*ng+i] = 0.0 - nx*ny;
            fh_u[8*ng+i] = 0.0 - nx*nz;        
            fh_u[9*ng+i] = 0.0;
            
            fh_u[10*ng+i] = 0.0;
            fh_u[11*ng+i] = 0.0 - ny*nx;
            fh_u[12*ng+i] = 1.0 - ny*ny;
            fh_u[13*ng+i] = 0.0 - ny*nz;        
            fh_u[14*ng+i] = 0.0;
                        
            fh_u[15*ng+i] = 0.0;
            fh_u[16*ng+i] = 0.0 - nz*nx;
            fh_u[17*ng+i] = 0.0 - nz*ny;
            fh_u[18*ng+i] = 1.0 - nz*nz;        
            fh_u[19*ng+i] = 0.0;
                        
            fh_u[20*ng+i] = 0.0;            
            fh_u[21*ng+i] = 0.0;
            fh_u[22*ng+i] = 0.0;            
            fh_u[23*ng+i] = 0.0;
            fh_u[24*ng+i] = 1.0;

            //fh_uh 
            fh_uh[0*ng+i]  = -1.0;
            fh_uh[1*ng+i]  = 0.0;
            fh_uh[2*ng+i]  = 0.0;
            fh_uh[3*ng+i]  = 0.0;        
            fh_uh[4*ng+i]  = 0.0;
            
            fh_uh[5*ng+i]  = 0.0;
            fh_uh[6*ng+i]  = -1.0;
            fh_uh[7*ng+i]  = 0.0;
            fh_uh[8*ng+i]  = 0.0;        
            fh_uh[9*ng+i]  = 0.0;
            
            fh_uh[10*ng+i]  = 0.0;
            fh_uh[11*ng+i]  = 0.0;
            fh_uh[12*ng+i] = -1.0;
            fh_uh[13*ng+i] = 0.0;              
            fh_uh[14*ng+i] = 0.0;
                        
            fh_uh[15*ng+i]  = 0.0;
            fh_uh[16*ng+i]  = 0.0;
            fh_uh[17*ng+i] = 0.0;
            fh_uh[18*ng+i] = -1.0;              
            fh_uh[19*ng+i] = 0.0;
            
            fh_uh[20*ng+i] = 0.0;
            fh_uh[21*ng+i] = 0.0;
            fh_uh[22*ng+i] = 0.0;
            fh_uh[23*ng+i] = 0.0;
            fh_uh[24*ng+i] = -1.0;                                                
        }            
    }
    else if (ib==4) { /* pressure */
        
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *anm  = new double [ng*nch*nch*nch];
        double *Anm  = new double [ng*nch*nch*nch];
        double *uinf  = new double [ng*nch];
        double *uinfu  = new double [ng*nch*nch];
        
        double pinf = (ui[4]-0.5)*(gam-1);
        
        getan(an,anm,uhg,nl,param,0,ng,nch,nd);
        getan(An,Anm,uhg,nl,param,1,ng,nch,nd);                        
        
        for (i=0; i<ng*(nch-1); i++) 
            uinf[i] = udg[i];        
        for (i=0; i<ng; i++) {
            n = (nch-1)*ng+i;                            
            uinf[n] = pinf/(gam-1) + 0.5*udg[ng+i]*udg[ng+i]/udg[i] 
                    + 0.5*udg[2*ng+i]*udg[2*ng+i]/udg[i] + 0.5*udg[3*ng+i]*udg[3*ng+i]/udg[i];
        }        
        
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
   else if (ib==8) { // Subsonic turbomachinery inlet: Impose p0, T0, alpha, theta. Extrapolate R- = u - 2c/(gamma-1)
       // TODO: Code SA BC. Code ALE.

       double Vn, p0_target, T0_target, alpha_target, theta_target, aux1, aux2;

       p0_target = ui[0];
       T0_target = ui[1];
       alpha_target = ui[2];
       if (nd == 3) {
           theta_target = ui[3];
       }

       for (i = 0; i < ng; i++) {
           fh[0*ng+i] = (gam-1.0)*uhg[(1+nd)*ng+i] - p0_target;
           fh[1*ng+i] = uhg[(1+nd)*ng+i]/uhg[0*ng+i] - T0_target;
           for (j = 0; j < nd; j++) {
               fh[0*ng+i] += 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               fh[1*ng+i] -= 0.5*((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           }
           fh[0*ng+i] *= M2;
           fh[1*ng+i] *= M2;

           fh[2*ng+i] = uhg[1*ng+i]*sin(alpha_target) - uhg[2*ng+i]*cos(alpha_target);

            fh[3*ng+i] = sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i])*cos(theta_target) - uhg[3*ng+i];

           // Extrapolate J- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: J- = u_n - 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[(1+nd)*ng+i] = - 2.0*sqrt((gam/(gam-1.0))*aux1) + 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[(1+nd)*ng+i] += -nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate J- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: J- = u_n - sqrt(gam*e/(gam-1))
           }
       }

       for (i = 0; i < ng * nch * nc; i++)
           fh_u[i] = 0.0;

       for (i = 0; i < ng * nch * nch; i++)
           fh_uh[i] = 0.0;

       for (i = 0; i < ng; i++) {
           fh_uh[(1+nd) * sz2 + 0 * ng + i] = gam-1.0;
           fh_uh[(1+nd) * sz2 + 1 * ng + i] = 1.0/uhg[0*ng+i];
           fh_uh[0 * sz2 + 1 * ng + i] = - uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           for (j = 0; j < nd; j++) {
               fh_uh[0 * sz2 + 0 * ng + i] -= 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               fh_uh[(1+j) * sz2 + 0 * ng + i] += (2.0-gam)*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               fh_uh[0 * sz2 + 1 * ng + i] += ((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i]);
               fh_uh[(1+j) * sz2 + 1 * ng + i] -= ((gam-1.0)/gam)*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           }
           fh_uh[0*sz2 + 0*ng+i] *= M2;
           fh_uh[(1+nd)*sz2 + 0*ng+i] *= M2;
           fh_uh[0*sz2 + 1*ng+i] *= M2;
           fh_uh[(1+nd)*sz2 + 1*ng+i] *= M2;
           for (j = 0; j < nd; j++) {
               fh_uh[(1+j) * sz2 + 0 * ng + i] *= M2;
               fh_uh[(1+j) * sz2 + 1 * ng + i] *= M2;
           }

           fh_uh[1 * sz2 + 2*ng+i] = sin(alpha_target);
           fh_uh[2 * sz2 + 2*ng+i] = - cos(alpha_target);

           fh_uh[1 * sz2 + 3*ng + i] = (uhg[1*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
           fh_uh[2 * sz2 + 3*ng + i] = (uhg[2*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target);
           fh_uh[3 * sz2 + 3*ng + i] = (uhg[3*ng+i]/sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i]))*cos(theta_target) - 1.0;

           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }

           fh_uh[0*sz2 + (1+nd)*ng+i] = - (- uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
           fh_u[0*sz2 + (1+nd)*ng+i] = (- udg[(1+nd)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);

           for (j=0; j<nd; j++) {
               fh_uh[0*sz2 + (1+nd)*ng+i] += - (uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[0*sz2 + (1+nd)*ng+i] += (udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
               fh_uh[0*sz2 + (1+nd)*ng+i] += nl[j*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               fh_u[0*sz2 + (1+nd)*ng+i] += - nl[j*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
               fh_uh[(1+j) * sz2 + (1+nd)*ng+i] += - nl[j*ng+i]/uhg[0*ng+i] + (uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[(1+j) * sz2 + (1+nd)*ng+i] += nl[j*ng+i]/udg[0*ng+i] - (udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
           }

           fh_uh[(1+nd)*sz2 + (1+nd)*ng+i] = - (1.0/uhg[0*ng+i])*sqrt((gam/(gam-1.0))/aux1);
           fh_u[(1+nd)*sz2 + (1+nd)*ng+i] = (1.0/udg[0*ng+i])*sqrt((gam/(gam-1.0))/aux2);
       }
   }
   else if (ib==9) { // Subsonic turbomachinery outflow: Impose p. Extrapolate R+ = u + 2c/(gamma-1), s and u_t
       // TODO: Code SA BC. Code ALE.

       double Vn, p_target, kinEnergy, kinEnergyHat, aux1, aux2;

       p_target = ui[0];

       for (i = 0; i < ng; i++) {
           kinEnergyHat = 0.0;
           kinEnergy = 0.0;
           for (j=0; j<nd; j++) {
               kinEnergyHat += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               kinEnergy += 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/udg[0*ng+i];
           }
           fh[0*ng+i] = M2 * ((gam-1.0) * (uhg[(1+nd)*ng+i]-kinEnergyHat) - p_target);

           // Extrapolate J+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: J+ = u_n + 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[1*ng+i] = 2.0*sqrt((gam/(gam-1.0))*aux1) - 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[1*ng+i] += nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate J+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: J+ = u_n + sqrt(gam*e/(gam-1))
           }


           // Extrapolate s = p / rho^gam [Dimensionless entropy: s = p / rho^gam , since we define s_ref := u_ref^2 / rho_ref^(gam-1)]
           // Note: p = (gam-1)*(rE - 0.5*r*u^2 - 0.5*r*v^2 - 0.5*r*w^2) [both dimensional and non-dimensional]
           fh[2*ng+i] = M2 * (uhg[(1+nd)*ng+i]-kinEnergyHat) / pow(uhg[0*ng+i],gam) - M2 * (udg[(1+nd)*ng+i]-kinEnergy) / pow(udg[0*ng+i],gam);


           // Extrapolate v_t:
           double t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, t2_norm, nl_norm;
           nl_norm = sqrt(nl[0*ng+i]*nl[0*ng+i] + nl[1*ng+i]*nl[1*ng+i] + nl[2*ng+i]*nl[2*ng+i]);
           if (abs(nl[0*ng+i])/nl_norm < 0.9) {
               t1x = 0.0;
               t1y = -nl[2*ng+i];
               t1z = nl[1*ng+i];
           }
           else {
               t1x = -nl[2*ng+i];
               t1y = 0.0;
               t1z = nl[0*ng+i];
           }
           t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
           t1x = t1x / t1_norm;
           t1y = t1y / t1_norm;
           t1z = t1z / t1_norm;

           t2x = t1z*nl[1*ng+i] - t1y*nl[2*ng+i];
           t2y = t1x*nl[2*ng+i] - t1z*nl[0*ng+i];
           t2z = t1y*nl[0*ng+i] - t1x*nl[1*ng+i];
           t2_norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
           t2x = t2x / t2_norm;
           t2y = t2y / t2_norm;
           t2z = t2z / t2_norm;

           fh[3*ng+i] = t1x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                        t1y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                        t1z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);

           fh[4*ng+i] = t2x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                        t2y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                        t2z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);
       }

       for (i = 0; i < ng * nch * nc; i++)
           fh_u[i] = 0.0;

       for (i = 0; i < ng * nch * nch; i++)
           fh_uh[i] = 0.0;

       for (i = 0; i < ng; i++) {
           kinEnergyHat = 0.0;
           kinEnergy = 0.0;
           for (j=0; j<nd; j++) {
               kinEnergyHat += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               kinEnergy += 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/udg[0*ng+i];

               fh_uh[0*sz2     + 0*ng+i] += M2 * 0.5*(gam-1.0)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               fh_uh[(1+j)*sz2 + 0*ng+i] = - M2 * (gam-1.0)*uhg[(1+j)*ng+i]/uhg[0*ng+i];
           }
           fh_uh[(1+nd)*sz2 + 0*ng+i] = M2 * (gam-1.0);


           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }

           fh_uh[0*sz2 + 1*ng+i] = (- uhg[(1+nd)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
           fh_u[0*sz2 + 1*ng+i] = - (- udg[(1+nd)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);

           for (j=0; j<nd; j++) {
               fh_uh[0*sz2 + 1*ng+i] += (uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[0*sz2 + 1*ng+i] += - (udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
               fh_uh[0*sz2 + 1*ng+i] += - nl[j*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               fh_u[0*sz2 + 1*ng+i] += nl[j*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
               fh_uh[(1+j) * sz2 + 1*ng+i] += nl[j*ng+i]/uhg[0*ng+i] - (uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux1);
               fh_u[(1+j) * sz2 + 1*ng+i] += - nl[j*ng+i]/udg[0*ng+i] + (udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i])) * sqrt((gam/(gam-1.0))/aux2);
           }

           fh_uh[(1+nd)*sz2 + 1*ng+i] = (1.0/uhg[0*ng+i])*sqrt((gam/(gam-1.0))/aux1);
           fh_u[(1+nd)*sz2 + 1*ng+i] = - (1.0/udg[0*ng+i])*sqrt((gam/(gam-1.0))/aux2);


           fh_uh[0*sz2 + 2*ng+i] = - M2 * gam * uhg[(1+nd)*ng+i] / pow(uhg[0*ng+i],gam+1.0) +
                                     M2 * (gam+1.0)*kinEnergyHat / pow(uhg[0*ng+i],gam+1.0);
           fh_u[0*sz2 + 2*ng+i] = M2 * gam * udg[(1+nd)*ng+i] / pow(udg[0*ng+i],gam+1.0) -
                                  M2 * (gam+1.0)*kinEnergy / pow(udg[0*ng+i],gam+1.0);

           for (k=0; k<nd; k++) {
               fh_uh[(1+k)*sz2 + 2*ng+i] = - M2 * uhg[(1+k)*ng+i] / pow(uhg[0*ng+i],gam+1.0);
               fh_u[(1+k)*sz2 + 2*ng+i] = M2 * udg[(1+k)*ng+i] / pow(udg[0*ng+i],gam+1.0);
           }

           fh_uh[(1+nd)*sz2 + 2*ng+i] = M2 / pow(uhg[0*ng+i],gam);
           fh_u[(1+nd)*sz2 + 2*ng+i] = - M2 / pow(udg[0*ng+i],gam);

           double t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, t2_norm, nl_norm;
           nl_norm = sqrt(nl[0*ng+i]*nl[0*ng+i] + nl[1*ng+i]*nl[1*ng+i] + nl[2*ng+i]*nl[2*ng+i]);
           if (abs(nl[0*ng+i])/nl_norm < 0.9) {
               t1x = 0.0;
               t1y = -nl[2*ng+i];
               t1z = nl[1*ng+i];
           }
           else {
               t1x = -nl[2*ng+i];
               t1y = 0.0;
               t1z = nl[0*ng+i];
           }
           t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
           t1x = t1x / t1_norm;
           t1y = t1y / t1_norm;
           t1z = t1z / t1_norm;

           t2x = t1z*nl[1*ng+i] - t1y*nl[2*ng+i];
           t2y = t1x*nl[2*ng+i] - t1z*nl[0*ng+i];
           t2z = t1y*nl[0*ng+i] - t1x*nl[1*ng+i];
           t2_norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
           t2x = t2x / t2_norm;
           t2y = t2y / t2_norm;
           t2z = t2z / t2_norm;

           fh_uh[0*sz2 + 3*ng+i] = - t1x*uhg[1*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                     t1y*uhg[2*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                     t1z*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           fh_uh[1*sz2 + 3*ng+i] = t1x/uhg[0*ng+i];
           fh_uh[2*sz2 + 3*ng+i] = t1y/uhg[0*ng+i];
           fh_uh[3*sz2 + 3*ng+i] = t1z/uhg[0*ng+i];

           fh_u[0*sz2 + 3*ng+i] = t1x*udg[1*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                  t1y*udg[2*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                  t1z*udg[3*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           fh_u[1*sz2 + 3*ng+i] = - t1x/udg[0*ng+i];
           fh_u[2*sz2 + 3*ng+i] = - t1y/udg[0*ng+i];
           fh_u[3*sz2 + 3*ng+i] = - t1z/udg[0*ng+i];

           fh_uh[0*sz2 + 4*ng+i] = - t2x*uhg[1*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                     t2y*uhg[2*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]) -
                                     t2z*uhg[3*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           fh_uh[1*sz2 + 4*ng+i] = t2x/uhg[0*ng+i];
           fh_uh[2*sz2 + 4*ng+i] = t2y/uhg[0*ng+i];
           fh_uh[3*sz2 + 4*ng+i] = t2z/uhg[0*ng+i];

           fh_u[0*sz2 + 4*ng+i] = t2x*udg[1*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                  t2y*udg[2*ng+i]/(udg[0*ng+i]*udg[0*ng+i]) +
                                  t2z*udg[3*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           fh_u[1*sz2 + 4*ng+i] = - t2x/udg[0*ng+i];
           fh_u[2*sz2 + 4*ng+i] = - t2y/udg[0*ng+i];
           fh_u[3*sz2 + 4*ng+i] = - t2z/udg[0*ng+i];
       }
   }
    else {                                
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);           
        exit(-1);                                    
    }                
}


void fbouonly_euler3d(double *fh, double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{        
    
    double gam, epslm, gam1, Minf, M2, Re;    
    double nx, ny, nz, tm;        
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
    else if(ib==2) { /* inviscid wall */

        double ur, uru, urv, urw, urE, hr, hru, hrv, hrw, hrE, unu;
        double uinf1, uinf2, uinf3, uinf4, uinf5;
        
        for (i=0; i<ng; i++) {
            // normal vector
            nx  = nl[0*ng+i];              
            ny  = nl[1*ng+i]; 
            nz  = nl[2*ng+i]; 
            
            // UDG
            ur  = udg[0*ng+i];
            uru = udg[1*ng+i];
            urv = udg[2*ng+i];
            urw = udg[3*ng+i];
            urE = udg[4*ng+i];
            
            // UH
            hr  = uhg[0*ng+i];
            hru = uhg[1*ng+i];
            hrv = uhg[2*ng+i];
            hrw = uhg[3*ng+i];
            hrE = uhg[4*ng+i];
                              
            // make uinf
            unu   = uru*nx + urv*ny + urw*nz;
            uinf1 = ur;
            uinf2 = uru - unu*nx;
            uinf3 = urv - unu*ny;
            uinf4 = urw - unu*nz;
            uinf5 = urE;            
            
            // fh = uinf - uh;
            fh[0*ng+i] = uinf1 - hr;
            fh[1*ng+i] = uinf2 - hru;    
            fh[2*ng+i] = uinf3 - hrv;
            fh[3*ng+i] = uinf4 - hrw;
            fh[4*ng+i] = uinf5 - hrE;
        }                    
    }    
    else if (ib==4) { /* pressure */
        
        double *an  = new double [ng*nch*nch];
        double *An  = new double [ng*nch*nch];
        double *uinf  = new double [ng*nch];
        
        double pinf = (ui[4]-0.5)*(gam-1);
        
        getanonly(an,uhg,nl,param,0,ng,nch,nd);
        getanonly(An,uhg,nl,param,1,ng,nch,nd);                        
        
        for (i=0; i<ng*(nch-1); i++) 
            uinf[i] = udg[i];        
        for (i=0; i<ng; i++) {
            n = (nch-1)*ng+i;                            
            uinf[n] = pinf/(gam-1) + 0.5*udg[ng+i]*udg[ng+i]/udg[i] 
                    + 0.5*udg[2*ng+i]*udg[2*ng+i]/udg[i] + 0.5*udg[3*ng+i]*udg[3*ng+i]/udg[i];
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
   else if (ib==8) { // Subsonic turbomachinery inlet: Impose p0, T0, alpha, theta. Extrapolate R- = u - 2c/(gamma-1)
       // TODO: Code SA BC. Code ALE.

       double Vn, p0_target, T0_target, alpha_target, theta_target, aux1, aux2;

       p0_target = ui[0];
       T0_target = ui[1];
       alpha_target = ui[2];
       if (nd == 3) {
           theta_target = ui[3];
       }

       for (i = 0; i < ng; i++) {
           fh[0*ng+i] = (gam-1.0)*uhg[(1+nd)*ng+i] - p0_target;
           fh[1*ng+i] = uhg[(1+nd)*ng+i]/uhg[0*ng+i] - T0_target;
           for (j = 0; j < nd; j++) {
               fh[0*ng+i] += 0.5*(2.0-gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               fh[1*ng+i] -= 0.5*((gam-1.0)/gam)*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
           }
           fh[0*ng+i] *= M2;
           fh[1*ng+i] *= M2;

           fh[2*ng+i] = uhg[1*ng+i]*sin(alpha_target) - uhg[2*ng+i]*cos(alpha_target);

            fh[3*ng+i] = sqrt(uhg[1*ng+i]*uhg[1*ng+i]+uhg[2*ng+i]*uhg[2*ng+i]+uhg[3*ng+i]*uhg[3*ng+i])*cos(theta_target) - uhg[3*ng+i];

           // Extrapolate J- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: J- = u_n - 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[(1+nd)*ng+i] = - 2.0*sqrt((gam/(gam-1.0))*aux1) + 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[(1+nd)*ng+i] += -nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate J- = u_n - 2c/(gamma-1) [Dimensionless Riemann invariant: J- = u_n - sqrt(gam*e/(gam-1))
           }
       }
   }
   else if (ib==9) { // Subsonic turbomachinery outflow: Impose p. Extrapolate R+ = u + 2c/(gamma-1), s and u_t
       // TODO: Code SA BC. Code ALE.

       double Vn, p_target, kinEnergy, kinEnergyHat, aux1, aux2;

       p_target = ui[0];

       for (i = 0; i < ng; i++) {
           kinEnergyHat = 0.0;
           kinEnergy = 0.0;
           for (j=0; j<nd; j++) {
               kinEnergyHat += 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/uhg[0*ng+i];
               kinEnergy += 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/udg[0*ng+i];
           }
           fh[0*ng+i] = M2 * ((gam-1.0) * (uhg[(1+nd)*ng+i]-kinEnergyHat) - p_target);

           // Extrapolate J+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: J+ = u_n + 2*sqrt(gam*e/(gam-1))
           aux1 = uhg[(1+nd)*ng+i]/uhg[0*ng+i];
           aux2 = udg[(1+nd)*ng+i]/udg[0*ng+i];
           for (j=0; j<nd; j++) {
               aux1 -= 0.5*uhg[(1+j)*ng+i]*uhg[(1+j)*ng+i]/(uhg[0*ng+i]*uhg[0*ng+i]);
               aux2 -= 0.5*udg[(1+j)*ng+i]*udg[(1+j)*ng+i]/(udg[0*ng+i]*udg[0*ng+i]);
           }
           fh[1*ng+i] = 2.0*sqrt((gam/(gam-1.0))*aux1) - 2.0*sqrt((gam/(gam-1.0))*aux2);
           for (j=0; j<nd; j++) {
               fh[1*ng+i] += nl[j*ng+i]*(uhg[(1+j)*ng+i]/uhg[0*ng+i] - udg[(1+j)*ng+i]/udg[0*ng+i]);        // Extrapolate J+ = u_n + 2c/(gamma-1) [Dimensionless Riemann invariant: J+ = u_n + sqrt(gam*e/(gam-1))
           }


           // Extrapolate s = p / rho^gam [Dimensionless entropy: s = p / rho^gam , since we define s_ref := u_ref^2 / rho_ref^(gam-1)]
           // Note: p = (gam-1)*(rE - 0.5*r*u^2 - 0.5*r*v^2 - 0.5*r*w^2) [both dimensional and non-dimensional]
           fh[2*ng+i] = M2 * (uhg[(1+nd)*ng+i]-kinEnergyHat) / pow(uhg[0*ng+i],gam) - M2 * (udg[(1+nd)*ng+i]-kinEnergy) / pow(udg[0*ng+i],gam);


           // Extrapolate v_t:
           double t1x, t1y, t1z, t2x, t2y, t2z, t1_norm, t2_norm, nl_norm;
           nl_norm = sqrt(nl[0*ng+i]*nl[0*ng+i] + nl[1*ng+i]*nl[1*ng+i] + nl[2*ng+i]*nl[2*ng+i]);
           if (abs(nl[0*ng+i])/nl_norm < 0.9) {
               t1x = 0.0;
               t1y = -nl[2*ng+i];
               t1z = nl[1*ng+i];
           }
           else {
               t1x = -nl[2*ng+i];
               t1y = 0.0;
               t1z = nl[0*ng+i];
           }
           t1_norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
           t1x = t1x / t1_norm;
           t1y = t1y / t1_norm;
           t1z = t1z / t1_norm;

           t2x = t1z*nl[1*ng+i] - t1y*nl[2*ng+i];
           t2y = t1x*nl[2*ng+i] - t1z*nl[0*ng+i];
           t2z = t1y*nl[0*ng+i] - t1x*nl[1*ng+i];
           t2_norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
           t2x = t2x / t2_norm;
           t2y = t2y / t2_norm;
           t2z = t2z / t2_norm;

           fh[3*ng+i] = t1x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                        t1y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                        t1z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);

           fh[4*ng+i] = t2x*(uhg[1*ng+i]/uhg[0*ng+i] - udg[1*ng+i]/udg[0*ng+i]) +
                        t2y*(uhg[2*ng+i]/uhg[0*ng+i] - udg[2*ng+i]/udg[0*ng+i]) +
                        t2z*(uhg[3*ng+i]/uhg[0*ng+i] - udg[3*ng+i]/udg[0*ng+i]);
       }
   }
    else {                        
        printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        printf("Boundary condition %d is not implemented yet.", ib);            
        exit(-1);                                    
    }                
}


