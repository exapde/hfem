#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdlib.h>

// Written by: C. Nguyen & P. Fernandez

using namespace std;

#include "../preprocessing/digasopre.h"
#include "../preprocessing/preprocessing.cpp"
#include "../assembly/assembly.cpp"
#include "../solvers/dgsolvers.cpp"

vector<Int> getndimslegacy(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, sysstruct &sys);

int main(int argc, char** argv) 
{
   
    if( argc == 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }            
        
    appstruct app;
    vector< masterstruct > master;
    meshstruct mesh;
    solstruct sol;
    sysstruct sys;
    elemstruct elem;
    tempstruct temp;
        
    app.filein  = string(argv[1]); 
    app.fileout  = string(argv[2]);    
    app.my_rank = 0;        
    app.nproc = 1;        
    app.debugmode = 0;     
    sys.my_rank = 0; 
    sys.noThreads = 1;
    mesh.my_rank = 0;
    
    preprocessing(mesh, master, app, sol, elem, sys, temp, string(argv[1]), string(argv[2]));            
    
    vector<elemstruct> elems;
    elems.resize(sys.noThreads);
    vector<tempstruct> temps;
    temps.resize(sys.noThreads);
    elems[0] = elem;
    temps[0] = temp;        
    for (Int i = 0; i < sys.noThreads; i++)
        elems[i].my_rank = 0;     
    
    vector<Int> ndims = getndimslegacy(mesh, master, app, sol, sys);
    
    Int e = mesh.elementtype[0];
    solveProblem(sys, &elems[0], mesh, master[e], sol, app, &temps[0], &ndims[0]);
                
    return 0;             
}

vector<Int> getndimslegacy(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, sysstruct &sys)
{
    /* Get dimensions from ndims array */
    vector<Int> ndims(100,0);
    
    ndims[0] = app.nd;
    ndims[1] = app.ncd;
    ndims[2] = mesh.nfemax;
    ndims[3] = mesh.nvemax;
    ndims[4] = mesh.nvfmax;
    ndims[5] = mesh.ne;
    ndims[6] = mesh.nf;
    ndims[7] = mesh.nv;
    ndims[8] = mesh.ndofuh;
    ndims[9] = mesh.npemax;
    ndims[10] = mesh.npfmax;
    ndims[11] = mesh.nmemax;
    ndims[12] = mesh.nmfmax;
    ndims[13] = mesh.ngemax;
    ndims[14] = mesh.ngfmax;
    ndims[15] = app.porder[0];
    ndims[16] = app.porder[0];
    ndims[17] = app.dirkOrder;
    ndims[18] = app.dirkStage;
    ndims[19] = app.nc;
    ndims[20] = app.ncu;
    ndims[21] = app.ncq;
    ndims[22] = app.ncp;
    ndims[23] = app.nch;
    //Int ns  = ndims[24];
    //Int nb  = ndims[25];
    ndims[26] = app.dt.size();
    ndims[27] = app.physicsparam.size();
    ndims[28] = app.flag.size();
    ndims[29] = app.factor.size();
    ndims[30] = app.problem.size();
    ndims[31] = sys.numEntities;
    ndims[32] = sys.numBlocks;
    ndims[33] = sys.maxBlocksPerRow;
    ndims[34] = sys.blkSize;
    ndims[35] = sys.BJ_nrows;        
    ndims[39] = sys.nproc;
    ndims[40] = sys.entpartpts.size();        
    ndims[41] = sys.nentrecv;
    ndims[42] = sys.nentsend;
    ndims[43] = sys.elempartpts.size();        
    ndims[44] = sys.nelemrecv;
    ndims[45] = sys.nelemsend;
    ndims[46] = sys.nnbsd;
//        Int nglobalEnt2entStart = ndims[47];
//        Int nglobalEnt2ent = ndims[48];
    ndims[49] = sys.nmatrecv;
    ndims[50] = sys.nmatsend;
    ndims[75] = mesh.ncfmax;
           
    return ndims;
}


