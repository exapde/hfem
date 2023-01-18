
#include <random>

void syntheticTurb(double Ix, double Iy, Int nx, Int ny, Int px, Int py, Int dt)
{
    // INPUTS / OUTPUTS:
    // Ix: Integral length scale in the x-direction.
    // Iy: Integral length scale in the y-direction.
    // nx: Number of nodes in Ix.
    // ny: Number of nodes in Iy.
    // px: Number of nodes in the x-direction.
    // py: Number of nodes in the y-direction.
    // dt: Time-step size.
    
    // COMMENTS:
    // 1) Without loss of generality, we assume the inlet plane is (x,y).
    // 2) Approach to deal with non-constant integral length scales (Ix, Iy) in different locations of the inlet plane:
    //    The distance between points x_{i,j} and x_{i+1,j+1} is proportional to the local (Ix, Iy). 
    //    I'm hoping the b's to achieve the desired two-point-correlation are this way the same everywhere...
    
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    
    switch (STmethod)
    case 1:
        {
        // Digital filter approach by Touber and Sandham.
        // Ref: E.Touber, N.D. Sandham, Large-eddy simulation of low-frequency unsteadiness in a turbulent shock-induced separation bubble, Theor. Comput. Fluid Dyn. (2009) 23: 79?107
        
        double I = sqrt(Ix * Iy);

        Int Nx = 2*nx;      // We set x2 here, but other integers >=2 are ok too
        Int Ny = 2*ny;      // We set x2 here, but other integers >=2 are ok too
        
        // 1. Compute b's (in first time step only):
        // TODO: If we restart the simulation, b needs to be computed even though it is not the first time step
        if (firstStep == 1) {
            for (Int i=0; i<=2*Nx; i++)
                for (Int j=0; j<=2*Ny; j++) {
                    bxTilde[i][j] = exp(-PI * (double) ((i-Nx)*(i-Nx)) / (double) (nx*nx));
                    byTilde[i][j] = exp(-PI * (double) ((j-Ny)*(j-Ny)) / (double) (ny*ny));
                }
            
            for (Int i=0; i<=2*Nx; i++)
                for (Int j=0; j<=2*Ny; j++) {
                    double den2_x = 0, den2_y = 0;
                    for (Int k=-Nx; k<=Nx; k++)
                        den2_x += bxTilde[i+k][j] * bxTilde[i+k][j];
                    for (Int k=-Ny; k<=Ny; k++)
                        den2_y += byTilde[i][j+k] * byTilde[i][j+k];
                    bx_tmp[i][j] = bxTilde[i][j] / sqrt(den2_x);
                    by_tmp[i][j] = byTilde[i][j] / sqrt(den2_y);
                }
            for (Int i=0; i<=2*Nx; i++)
                for (Int j=0; j<=2*Ny; j++)
                    b[i][j] = bx[i][j] * by[i][j];
        }
        
        // 2. Generate random r's
        for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                rx[i][j] = distribution(generator);
                ry[i][j] = distribution(generator);
                rz[i][j] = distribution(generator);
            }
        
        // 3. Compute perturbation velocities
        for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                uPrime[i][j] = 0.0;
                vPrime[i][j] = 0.0;
                wPrime[i][j] = 0.0;
                for (Int ii=0; ii<=2*Nx; ii++) {
                    if (i-Nx+ii >= 0 && i-Nx+ii < px) {
                        for (Int jj=0; jj<=2*Ny; jj++) {
                            if (j-Ny+jj >= 0 && j-Ny+jj < py) {
                                uPrime[i][j] += bx[ii][jj]*rx[i-Nx+ii][j-Ny+jj];
                                vPrime[i][j] += bx[ii][jj]*ry[i-Nx+ii][j-Ny+jj];
                                wPrime[i][j] += bx[ii][jj]*rz[i-Nx+ii][j-Ny+jj];
                            }
                        }
                    }
                }
            }
        
        // 4. Compute velocities at next time step:
        if (firstStep == 0) {
            for (Int i=0; i<px; i++)
                for (Int j=0; j<py; j++) {
                    double tau = I / uMean[i][j];       // Lagrangian time scale
                    uPrime[i][j] = uPrimeOld[i][j] * exp(-PI * dt / (2.0 * tau)) + uPrime[i][j] * sqrt(1.0 - exp(-PI * dt / tau));
                    vPrime[i][j] = vPrimeOld[i][j] * exp(-PI * dt / (2.0 * tau)) + vPrime[i][j] * sqrt(1.0 - exp(-PI * dt / tau));
                    wPrime[i][j] = wPrimeOld[i][j] * exp(-PI * dt / (2.0 * tau)) + wPrime[i][j] * sqrt(1.0 - exp(-PI * dt / tau));
                }
        }
       for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                u[i][j] = uMean[i][j] + uPrime[i][j];
                v[i][j] = vMean[i][j] + vPrime[i][j];
                w[i][j] = wMean[i][j] + wPrime[i][j];
            }
    
        // 5. Compute temperature at next time step - Strong Reynolds analogy (SRA)
        for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                Tprime[i][j] = - gam1 * M2 * Tmean[i][j] * uPrime[i][j] / uMean[i][j];
                T[i][j] = Tmean[i][j] + Tprime[i][j];
            }
    
        // 6. Compute density at next time step
        // From the BLA, pressure is constant across the BL and hence pressure fluctuations are negligible compared to the velocity, density, and temperature fluctuations.
        for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                rhoPrime[i][j] = - Tprime[i][j] * (rhoMean[i][j] / Tmean[i][j]);
                rho[i][j] = rhoMean[i][j] + rhoPrime[i][j];
            }
        
        // 7. Compute UH:
        for (Int i=0; i<px; i++)
            for (Int j=0; j<py; j++) {
                UH[0][i][j] = rho[i][j];
                UH[1][i][j] = rho[i][j] * u[i][j];
                UH[2][i][j] = rho[i][j] * v[i][j];
                UH[3][i][j] = rho[i][j] * w[i][j];
                UH[4][i][j] = rho[i][j] * T[i][j] + 0.5 * rho[i][j] * (u[i][j]*u[i][j] + v[i][j]*v[i][j] + w[i][j]*w[i][j]);
            }
        
        // 8. Map UH to desired format in the true UH vector:
                
        break;
    }
    
    
}