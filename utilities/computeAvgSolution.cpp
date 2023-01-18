#ifndef __COMPUTEAVERAGESOLUTION
#define __COMPUTEAVERAGESOLUTION

// Written by: C. Nguyen & P. Fernandez

void computeAvgSolution(solstruct &sol, Int timeStep)
{
    Int inc = 1, i;
    Int lenUDG = sol.UDG.size();
    Int lenUH = sol.UH.size();
    
    if (timeStep >= 1) {
        for (i = 0; i < lenUDG; i++) {
            sol.UDG_avg[i] = sol.UDG_avg[i] * (timeStep-1) / timeStep;
            sol.UDG_avg[i] += sol.UDG[i] / timeStep;
        }
        for (i = 0; i < lenUH; i++) {
            sol.UH_avg[i] = sol.UH_avg[i] * (timeStep-1) / timeStep;
            sol.UH_avg[i] += sol.UH[i] / timeStep;
        } 
    }
    else
        error("Error No. G6HJ9C in computeAvgSolution.\n");
}

#endif
