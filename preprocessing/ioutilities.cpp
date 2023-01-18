#ifndef __IOUTILITIES
#define __IOUTILITIES

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

void printiarray(vector<Int> &a)
{    
    Int m = (Int) a.size();
    for (Int i=0; i<m; i++)
        cout << a[i] << "   ";
    cout << endl;
}

void print1iarray(Int* a, Int m)
{    
    for (Int i=0; i<m; i++)
        cout << a[i] << "   ";
    cout << endl;
}

void print2iarray(Int* a, Int m, Int n)
{
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3iarray(Int* a, Int m, Int n, Int p)
{    
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void printdarray(vector<double> &a, int prec)
{
    Int m = (Int) a.size();
    //cout.precision(prec);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void print1darray(double* a, Int m)
{
    //cout.precision(4);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void print2darray(double* a, Int m, Int n)
{
    //cout.precision(4);
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << scientific << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3darray(double* a, Int m, Int n, Int p)
{
    cout.precision(8);
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << scientific << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}


void readcarray(ifstream &in, vector<char> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(char)*N );
    }
}

void readdarray(ifstream &in, vector<double> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(double)*N );
    }
}

void readiarray(ifstream &in, vector<Int> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        in.read( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
    }
}

void readiarrayfromdouble(ifstream &in, vector<Int> &a, Int N)
{
    if (N>0) {
        a.resize(N);
        double read;
        for (unsigned i = 0; i < N; i++) {
            in.read( reinterpret_cast<char*>( &read ), sizeof read );
            a[i] = (Int) round(read);
        }
    }
}

void writeiarray(ofstream &out, vector<Int> &a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
}

void writedarray(ofstream &out, vector<double> &a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

void writeiarray(ofstream &out, vector<Int> &a)
{
    Int N = (Int) a.size();
    if (N > 0)
        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );
}

void writedarray(ofstream &out, vector<double> &a)
{
    Int N = (Int) a.size();
    if (N > 0)
        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

void writeiarraytodouble(ofstream &out, vector<Int> &a, Int N)
{
    double b;
    for (unsigned i = 0; i < N; i++) {
        b = (double) a[i];
        out.write( reinterpret_cast<char*>( &b ), sizeof(double) );
    }
}

void writedarray2out(ofstream &out, double *a, Int N)
{
    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );
}

void writedarray2file(string filename, double *a, Int N)
{
    // Open file to read
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }

    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(double) * N );

    out.close();
}

void writeiarray2file(string filename, Int *a, Int N)
{
    // Open file to read
    ofstream out(filename.c_str(), ios::out | ios::binary);
            
    if (!out) {
        error("Unable to open file " + filename);
    }

    out.write( reinterpret_cast<char*>( &a[0] ), sizeof(Int) * N );

    out.close();
}

#ifdef HAVE_MPI       
void mpiwriteiarray2file(string filename, Int *a, Int N, Int my_rank, Int nproc)
{
    // Determine Nall and Nstart
    Int *Nall = new int[nproc];
    MPI_Allgather(&N, 1, MPI_INT, Nall, 1, MPI_INT, MPI_COMM_WORLD);    
    Int Nstart = 0;
    for (int i = 0; i<my_rank; i++)
        Nstart = Nstart + Nall[i];
    
    // Write Nall into file
    writeiarray2file(filename, &Nall[0], nproc);
    
    // then write array a into the same file in parallel
    char *fileout;
    fileout = new char[filename.size() + 1];
    //memcpy(fileout, filename.c_str(), filename.size() + 1);   
    filename.copy(fileout,filename.size() + 1);    

    MPI_File file;
    MPI_Status status;          
    MPI_File_open(MPI_COMM_WORLD, fileout, MPI_MODE_CREATE|MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &file);                    

    MPI_Offset offset = sizeof(Int)*(Nstart + nproc);  // shift by nproc to account for Nall  
    MPI_File_seek(file, offset, MPI_SEEK_SET);        
    MPI_File_write(file, &a[0], N, MPI_INT, &status);
    MPI_File_close(&file);                                           
    
    delete [] fileout;
    delete [] Nall;
}

void mpiwritedarray2file(string filename, double *a, Int N, Int my_rank, Int nproc)
{
    // Determine Nall and Nstart
    Int *Nall = new int[nproc];
    MPI_Allgather(&N, 1, MPI_INT, Nall, 1, MPI_INT, MPI_COMM_WORLD);    
    Int Nstart = 0;
    for (int i = 0; i<my_rank; i++)
        Nstart = Nstart + Nall[i];
    
    // then write array a into the same file in parallel
    char *fileout;
    fileout = new char[filename.size() + 1];
    //memcpy(fileout, filename.c_str(), filename.size() + 1);   
    filename.copy(fileout,filename.size() + 1);    
    MPI_File file;
    MPI_Status status;          
    MPI_File_open(MPI_COMM_WORLD, fileout, MPI_MODE_CREATE|MPI_MODE_WRONLY,
                      MPI_INFO_NULL, &file);                    

    MPI_Offset offset = sizeof(double)*Nstart;  
    MPI_File_seek(file, offset, MPI_SEEK_SET);        
    MPI_File_write(file, &a[0], N, MPI_DOUBLE, &status);
    MPI_File_close(&file);                                           
    
    delete [] fileout;
    delete [] Nall;
}
#endif         

void solstruct::writeSol2File(string filename, Int ne, Int n, Int* ndims, Int writeQflag)
{
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];

    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        cout <<"Unable to open file" << filename << endl;
    }

    if (out) {
        // Write udg
        if (writeQflag == 0) {
            Int len = npv*ncu, inc = 1;
            double a = (double) npv*ncu*ne;
            double b = (double) n;
            out.write( reinterpret_cast<char*>( &a ), sizeof(double) );
            out.write( reinterpret_cast<char*>( &b ), sizeof(double) );
            for (int i=0; i<ne; i++) {
                DCOPY(&len, &UDG[i*npv*nc], &inc, &UDG2Write[0], &inc);            
                out.write( reinterpret_cast<char*>( &UDG2Write[0] ), sizeof(double) * npv*ncu );
            }
        }
        else if (writeQflag == 1) {
            double a = (double) npv*nc*ne;
            double b = (double) n;
            out.write( reinterpret_cast<char*>( &a ), sizeof(double) );
            out.write( reinterpret_cast<char*>( &b ), sizeof(double) );
            out.write( reinterpret_cast<char*>( &UDG[0] ), sizeof(double) * npv*nc*ne );
        }
        else {
            printf("writeQflag has invalid value in writeSol2File.\n");
            exit(-1);
        }

        // Write uh
        out.write( reinterpret_cast<char*>( &UH[0] ), sizeof(double) * n );
    }

    // Close file
    out.close();
}

void solstruct::writeAvgSol2File(string filename_avg, Int ne, Int n, Int* ndims, Int writeQflag)
{
    Int npv = ndims[9];
    Int nc  = ndims[19];
    Int ncu = ndims[20];

    // Open file to write
    ofstream out(filename_avg.c_str(), ios::out | ios::binary);

    if (!out) {
        cout <<"Unable to open file" << filename_avg << endl;
    }

    // If this is a valid file
    if (out) {
        // write udg_avg
        if (writeQflag == 0) {
            Int len = npv*ncu, inc = 1;
            for (int i=0; i<ne; i++) {
                DCOPY(&len, &UDG_avg[i*npv*nc], &inc, &UDG2Write[i*npv*ncu], &inc);
            }
            out.write( reinterpret_cast<char*>( &UDG2Write[0] ), sizeof(double) * npv*ncu*ne );
        }
        else if (writeQflag == 1) {
            out.write( reinterpret_cast<char*>( &UDG_avg[0] ), sizeof(double) * npv*nc*ne );
        }
        else {
            printf("writeQflag has invalid value in writeSol2File.\n");
            exit(-1);
        }

        // write uh_avg
        out.write( reinterpret_cast<char*>( &UH_avg[0] ), sizeof(double) * n );
    }

    // Close file
    out.close();
}

void writeTimeStepSize2File(string filename, Int timeStep, Int DIRKstage, double time, double dt)
{
    ofstream out(filename.c_str(), ios::out | ios::app);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else
        out << NumberToString(timeStep) + "\t" + NumberToString(DIRKstage) + "\t" + NumberToString(time) + "\t" + NumberToString(dt) + "\n";
    out.close();
}

void writeScalarField2File(string filename, double* field, Int ne, Int* ndims)
{
    Int npv = ndims[9];

    ofstream out(filename.c_str(), ios::out | ios::binary);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else if (out)
        out.write( reinterpret_cast<char*>( &field[0] ), sizeof(double) * npv*ne );
    out.close();
}

#endif