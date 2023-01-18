#ifndef __PREPROCESSING
#define __PREPROCESSING

#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "readbinaryfiles.cpp"
#include "initializestructs.cpp"

#ifdef HAVE_MPI
void dmdparallel(sysstruct &sys, meshstruct &mesh, appstruct &app, solstruct &sol, dmdstruct &dmd);
#endif

#ifndef HAVE_MPI              
void dmdserial(sysstruct &sys, meshstruct &mesh, appstruct &app);
#endif

void preprocessing(meshstruct &mesh, vector< masterstruct > &master, appstruct &app, solstruct &sol, elemstruct &elems, sysstruct &sys, tempstruct &temps, string filein, string fileout) 
{   
    // read from binary files to construct app, master, mesh, sol
    readInput(app, master, mesh, sol, filein);                    
        
    if (app.nproc>1) {
#ifdef HAVE_MPI           
        dmdstruct dmd;
        dmdparallel(sys, mesh, app, sol, dmd);    
        if (app.debugmode == 1) 
            writeOutput(app, master, mesh, sol, dmd, fileout);

            //vector<Int> ndims = getndims(mesh, master, app, sol, sys);
            //printiarray(ndims);     

//         Int N = dmd.elempartpts[0]+dmd.elempartpts[1];
//         string filename1 = fileout + "elempart" + ".bin";
//         mpiwriteiarray2file(filename1, &dmd.elempart[0], N, app.my_rank, app.nproc);                
// 
//         N = dmd.entpartpts[0]+dmd.entpartpts[1];
//         string filename2 = fileout + "entpart" + ".bin";
//         mpiwriteiarray2file(filename2, &dmd.entpart[0], N, app.my_rank, app.nproc);

        Int N = dmd.elempartpts[0]+dmd.elempartpts[1];
        string filename1 = fileout + "elempart" +  "_np" + NumberToString(app.my_rank) + ".bin";            
        writeiarray2file(filename1, &dmd.elempart[0], N);
        
        N = dmd.entpartpts[0]+dmd.entpartpts[1];
        string filename2 = fileout + "entpart" +  "_np" + NumberToString(app.my_rank) + ".bin";            
        writeiarray2file(filename2, &dmd.entpart[0], N);
        
        // clear memory
        cleardmdstruct(dmd); 
#endif                
    }
    else {
#ifndef HAVE_MPI              
        dmdserial(sys, mesh, app);    
#endif    
    }
    
    // allocate memory     
    initializeStructs(mesh, master, app, sol, elems, sys, temps);        
}


#endif