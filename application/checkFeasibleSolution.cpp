#ifndef __CHECKFEASIBLESOLUTION
#define __CHECKFEASIBLESOLUTION

#include "FM/checkFeasibleSolution_FM.cpp"

// Written by: C. Nguyen & P. Fernandez

Int checkFeasibleSolution(double *UDG, double *UH, appstruct &app, Int* ndims)
{
    int feasibleSolution_local, feasibleSolution_global;
    switch (app.appname) {
        case 0:
            feasibleSolution_local = (int) checkFeasibleSolution_FM(UDG, UH, app, ndims);
            break;
        case 1:
            feasibleSolution_local = (int) checkFeasibleSolution_FM(UDG, UH, app, ndims);
            break;
        case 3:
            feasibleSolution_local = (int) checkFeasibleSolution_FM(UDG, UH, app, ndims);
            break;
        default: {
            printf("Application not implemented (appname = %d)\n",app.appname);
            exit(-1);
        }
    }
    
#ifdef  HAVE_MPI
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int root = 0, count = 1;
    int my_rank = (int) app.my_rank;
    int * feasibleSolution_locals;
    if (my_rank == root)
        feasibleSolution_locals = new int[world_size];
    
    MPI_Gather(&feasibleSolution_local, count, MPI_INT, &feasibleSolution_locals[0], count, MPI_INT, root, MPI_COMM_WORLD);
    
    if (my_rank == root) {
        feasibleSolution_global = 1;
        for (Int i = 0; i < world_size; i++) {
            if (feasibleSolution_locals[i] == 0) {
                feasibleSolution_global = 0;
                break;
            }
        }
        delete[] feasibleSolution_locals;
    }
    
    MPI_Bcast(&feasibleSolution_global, count, MPI_INT, root, MPI_COMM_WORLD);
#endif
    
    return (Int) feasibleSolution_global;
}

#endif
