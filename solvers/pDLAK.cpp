#ifndef __PDLAK
#define __PDLAK

// Written by: C. Nguyen & P. Fernandez

#define THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE 2500       // Ratio "parallel overhead per thread" vs. "time to move a double to cache"

// Some useful info to define "THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE":
// - Overhead of "parallel for" and "parallel" in OpenMP: 1-10 um / thread (depending on source)
// - 2 GB/s-core of effective memory bandwidth corresponds to 4e-3 us/double (also, 2.5e5 doubles/s)
// - Hence: overhead "parallel for" / "time to load double in cache" = 250-2500 / thread

// Take-away: Threading outperforms if no_th \in (sz -/+ sqrt(sz*sz - 4*RATIO*sz)) / 2*RATIO, where no_th = "No. threads", sz = "No. doubles", and RATIO = THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE

/* DOUBLE PRECISION ROUTINES */

Int decideThreadingLevel1(Int sz1, Int noThreads, Int precision)
{
    // precision: Double: 1. Float: 2
    
    Int RATIO = THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE, threading;
    double sz = (double) sz1 / (double) precision;
    
    if ((noThreads < 2) || (sz < 4*RATIO) || (noThreads < (sz - sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (noThreads > (sz + sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (sz1 < 5*noThreads))
        threading = 0;
    else
        threading = 1;
    return threading;
}

Int decideThreadingLevel2(Int noThreads, Int sz1, Int sz2, Int precision)
{
    // precision: Double: 1. Float: 2
    
    Int RATIO = THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE, threading;
    double sz = (double) ((double) sz1*sz2) / (double) precision;
    
    if ((noThreads < 2) || (sz < 4*RATIO) || (noThreads < (sz - sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (noThreads > (sz + sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (max(sz1,sz2) < 5*noThreads))
        threading = 0;
    else
        threading = 1;
    return threading;
}

Int decideThreadingLevel3(Int noThreads, Int sz1, Int sz2, Int sz3, Int precision)
{
    // precision: Double: 1. Float: 2
    
    Int RATIO = THREAD_OVERHEAD_VS_TIME_CACHE_DOUBLE, threading;
    double sz = ((double) sz1*sz2*sz3) / (double) precision;
    
    if ((noThreads < 2) || (sz < 4*RATIO) || (noThreads < (sz - sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (noThreads > (sz + sqrt(sz*sz-4*sz*RATIO)) / (2*RATIO)) || (max(max(sz1,sz2),sz3) < 5*noThreads))
        threading = 0;
    else
        threading = 1;
    return threading;
}

void splitArray1D(Int arrayLen, Int this_thread, Int noThreads, Int* chunkLen, Int* startPosition)
{
    *startPosition = ((arrayLen * this_thread) / noThreads);
    Int endPosition = ((arrayLen * (this_thread+1)) / noThreads);
    
    if (this_thread == (noThreads-1)) {   // This is reduntant unless we go over the range in which double represent integers exactly
        if (endPosition != arrayLen) {
            printf("Warning 4TG67N in splitArray1D.\n");
            endPosition = arrayLen;
        }
    }
    
    *chunkLen = endPosition - *startPosition;
}

void DGEMV_OpenMP(char* ch, Int* m, Int* n, double* alpha, double* A, Int* lda, 
        double* x, Int* incx, double* beta, double* y, Int* incy, vector<vector<double> >* ys, Int noThreads) 
{
    /* y <- alpha*X^T * x + beta y */
    
    // TODO: Implement more sophisticated ways of dividing the matrix-vector product (e.g. dividing the matrix in both dimensions)
    Int oneInt = 1, i, j, postProc, dim2div, lenDim2div, lenDimNot2div;
    Int threading = decideThreadingLevel2(*m, *n, noThreads, 1);
    
    if (!threading)     // Do not divide work among threads
        DGEMV(ch, m, n, alpha, A, lda, x, incx, beta, y, incy);
    else {
        // Do divide work among threads
        
        // TODO: Improve the division criterion.
        if (*m > *n) {    // Divide work along first dimension
            dim2div = 1;
            lenDim2div = *m;
            lenDimNot2div = *n;
        }
        else {            // Divide work along second dimension
            dim2div = 2;
            lenDim2div = *n;
            lenDimNot2div = *m;
        }
        
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            Int *m_, *n_, *incy_;
            double *A_, *x_, *beta_, *y_;
            
            splitArray1D(lenDim2div, (Int) this_thread, noThreads, &myChunkSz, &startPt);
            
            if (dim2div == 1) {
                m_ = &myChunkSz;
                n_ = n;
                A_ = &A[startPt];
                if (*ch == 'N') {
                    x_ = x;
                    y_ = &y[startPt*(*incy)];
                    beta_ = beta;
                    incy_ = incy;
                    postProc = 0;
                }
                else if (*ch == 'T') {
                    x_ = &x[startPt*(*incx)];
                    y_ = &ys[0][this_thread][0];
                    incy_ = &oneInt;
                    *beta_ = 0.0;
                    postProc = 1;
                }
            }
            else if (dim2div == 2) {
                m_ = m;
                n_ = &myChunkSz;
                A_ = &A[startPt*(*lda)];
                if (*ch == 'N') {
                    x_ = &x[startPt*(*incx)];
                    y_ = &ys[0][this_thread][0];
                    incy_ = &oneInt;
                    *beta_ = 0.0;
                    postProc = 1;
                }
                else if (*ch == 'T') {
                    x_ = x;
                    y_ = &y[startPt*(*incy)];
                    incy_ = incy;
                    beta_ = beta;
                    postProc = 0;
                }
            }
            
            // Compute thread contribution to matrix-vector product:
            DGEMV(ch, m_, n_, alpha, A_, lda, x_, incx, beta_, y_, incy_);
// // //         }
        
        // Gather final result of matrix-vector product:
        if (postProc == 1) {
            for (i = 0; i < lenDimNot2div; i++)
                y[i*(*incy)] = (*beta) * y[i*(*incy)] + ys[0][0][i];
            for (j = 1; j < noThreads; j++)
                for (i = 0; i < lenDimNot2div; i++)
                    y[i*(*incy)] += ys[0][j][i];
        }
    }
}

// void DGEMVn_OpenMP(char* chn, Int* m, Int* n, double* alpha, double* A, Int* lda, 
//         double* x, Int* incx, double* beta, double* y, Int* incy, Int noThreads)
// {
//     /* y = alpha*X * x + beta y */
//     
//     // This OpenMP implementation assumes m >>> n
//     #pragma omp parallel num_threads(noThreads)
//     {
//         Int this_thread = (Int) omp_get_thread_num();
//         Int myChunkSz, startPt;
//         
//         splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);
// 
//         DGEMV(chn, &myChunkSz, n, alpha, &A[startPt], lda, x,
//               incx, beta, &y[startPt*(*incy)], incy);
//     }
// }
// 
// void DGEMVt_OpenMP(char* cht, Int* m, Int* n, double* alpha, double* A, Int* lda, 
//         double* x, Int* incx, double* beta, double* y, Int* incy, vector<vector<double> >* ys, Int noThreads) 
// {
//     // TODO: Not sure if this routines works properly, but that's irrelevant since it is not used anymore.
// 
//     /* y = alpha*X^T * x + beta y */
//     
//     // This OpenMP implementation assumes m >>> n
//     Int oneInt = 1, i, j;
//     double zero = 0.0;
//     
//     #pragma omp parallel num_threads(noThreads)
//     {
//         Int this_thread = (Int) omp_get_thread_num();
//         Int myChunkSz, startPt;
//         
//         splitArray1D((*m), (Int) this_thread, noThreads, &myChunkSz, &startPt);
// 
//         DGEMV(cht, &myChunkSz, n, alpha, &A[startPt], lda, 
//                 &x[startPt*(*incx)], incx, &zero, &ys[0][this_thread][0], &oneInt);
//     }
//     for (i = 0; i < (*n); i++)
//         y[i*(*incy)] = (*beta) * y[i*(*incy)] + ys[0][0][i];
//     for (j = 1; j < noThreads; j++)
//         for (i = 0; i < (*n); i++)
//             y[i*(*incy)] += ys[0][j][i];
// }
// 
// void DGEMV_OpenMP(char* ch, Int* m, Int* n, double* alpha, double* A, Int* lda, 
//         double* x, Int* incx, double* beta, double* y, Int* incy, vector<vector<double> >* ys, Int noThreads)
// {
//     if (*ch == 'N')
//         DGEMVn_OpenMP(ch, m, n, alpha, A, lda, 
//             x, incx, beta, y, incy, noThreads);
//     else if (*ch == 'T')
//         DGEMVt_OpenMP(ch, m, n, alpha, A, lda, 
//             x, incx, beta, y, incy, ys, noThreads);
// }

void DSCAL_OpenMP(Int* m, double* alpha, double* x, Int* inc, Int noThreads)
{
    /* x <- alpha*x */
    Int threading = decideThreadingLevel1(*m, noThreads, 1);
      
    if (!threading)
        DSCAL(m, alpha, x, inc);
    else {
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            DSCAL(&myChunkSz, alpha, &x[startPt*(*inc)], inc);
// // //         }
    }
}

void DAXPY_OpenMP(Int* m, double* alpha, double* x, Int* incx, double* y, Int* incy, Int noThreads)
{
    /* y <- alpha * x + y */
    
    Int threading = decideThreadingLevel1(*m, noThreads, 1);
    
    if (!threading)
        DAXPY(m, alpha, x, incx, y, incy);
    else {
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            DAXPY(&myChunkSz, alpha, &x[startPt*(*incx)], incx, &y[startPt*(*incy)], incy);
// // //         }
    }
}

double DNRM2_OpenMP(Int* m, double* x, Int* incx, Int noThreads)
{
    // TODO: This function has not been validated.
    Int threading = decideThreadingLevel1(*m, noThreads, 1);
    double x2 = 0.0;
    
    if (!threading) {
        x2 = DNRM2(m, x, incx);
        x2 = x2*x2;
    }
    else {
// // //         #pragma omp parallel num_threads(noThreads) reduction(+:x2)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            if (myChunkSz > 0)
                x2 = DNRM2(&myChunkSz, &x[startPt*(*incx)], incx);
            else
                x2 = 0.0;
            x2 = x2*x2;
// // //         }
    }
    
    return sqrt(x2);
}

double DDOT_OpenMP(Int* m, double* x, Int* incx, double* y, Int* incy, Int noThreads)
{
    Int threading = decideThreadingLevel1(*m, noThreads, 1);
    double dot = 0.0;
    
    if (!threading)
        dot = DDOT(m, x, incx, y, incy);
    else {
// // //         #pragma omp parallel num_threads(noThreads) reduction(+:dot)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D((*m), (Int) this_thread, noThreads, &myChunkSz, &startPt);
                
            if (myChunkSz > 0)
                dot = DDOT(&myChunkSz, &x[startPt*(*incx)], incx, &y[startPt*(*incy)], incy);
            else
                dot = 0.0;
// // //         }
    }
    
    return dot;
}

#ifdef  HAVE_MPI

double DDOTMPI(Int m, double* x, Int incx, double* y, Int incy, Int noThreads) 
{
//     double local_dot = DDOT(&m, x, &incx, y, &incy);
    double local_dot = DDOT_OpenMP(&m, x, &incx, y, &incy, noThreads);
    
    double global_dot;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_dot;
}

void DGEMVMPI(char ch, Int m, Int n, double alpha, double* A, Int lda, double* x, Int incx, 
        double beta, double* y, Int incy, double* ylocal, vector<vector<double> >* ylocals, Int noThreads) 
{
    // TODO: I think this function may not work if beta != 0. Same with other MPI functions here.

    /* y = alpha*X^T * x + beta y */
//     DGEMV(&cht, &m, &n, &alpha, A, &lda, x, &incx, &beta, ylocal, &incy);
//     DGEMVt_OpenMP(&cht, &m, &n, &alpha, A, &lda, x, &incx, &beta, ylocal, &incy, ylocals, noThreads);
    DGEMV_OpenMP(&ch, &m, &n, &alpha, A, &lda, x, &incx, &beta, ylocal, &incy, ylocals, noThreads);

    MPI_Allreduce(ylocal, y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void DGEMMMPI(char cha, char chb, Int m, Int n, Int k, double alpha, double* A, Int lda,
              double* B, Int ldb, double beta, double* C, Int ldc, double* Clocal, Int noThreads)
{
    // TODO: Implement a threaded version of this function.
    
    /* C = alpha* opt(A) * opt(B) + beta C */
    DGEMM(&cha, &chb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, Clocal, &ldc);
    
    Int p = m*n;
    MPI_Allreduce(Clocal, C, p, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

#endif


/* SINGLE PRECISION ROUTINES */

void SGEMV_OpenMP(char* ch, Int* m, Int* n, float* alpha, float* A, Int* lda, 
        float* x, Int* incx, float* beta, float* y, Int* incy, vector<vector<float> >* ys, Int noThreads) 
{
    /* y <- alpha*X^T * x + beta y */
    // TODO: Implement more sophisticated ways of dividing the matrix-vector product (e.g. dividing the matrix in both dimensions)
    
    Int oneInt = 1, i, j, postProc, dim2div, lenDim2div, lenDimNot2div;
    Int threading = decideThreadingLevel2(*m, *n, noThreads, 2);

    if (!threading)     // Do not divide work among threads
        SGEMV(ch, m, n, alpha, A, lda, x, incx, beta, y, incy);
    else {      // Do divide work among threads
        // TODO: Improve the division criterion.
        if (*m > *n) {    // Divide work along first dimension
            dim2div = 1;
            lenDim2div = *m;
            lenDimNot2div = *n;
        }
        else {            // Divide work along second dimension
            dim2div = 2;
            lenDim2div = *n;
            lenDimNot2div = *m;
        }
        
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            Int *m_, *n_, *incy_;
            float *A_, *x_, *beta_, *y_;

            splitArray1D(lenDim2div, (Int) this_thread, noThreads, &myChunkSz, &startPt);
            
            if (dim2div == 1) {
                m_ = &myChunkSz;
                n_ = n;
                A_ = &A[startPt];
                if (*ch == 'N') {
                    x_ = x;
                    y_ = &y[startPt*(*incy)];
                    incy_ = incy;
                    beta_ = beta;
                    postProc = 0;
                }
                else if (*ch == 'T') {
                    x_ = &x[startPt*(*incx)];
                    y_ = &ys[0][this_thread][0];
                    incy_ = &oneInt;
                    *beta_ = 0.0;
                    postProc = 1;
                }
            }
            else if (dim2div == 2) {
                m_ = m;
                n_ = &myChunkSz;
                A_ = &A[this_thread*(*lda)];
                if (*ch == 'N') {
                    x_ = &x[startPt*(*incx)];
                    y_ = &ys[0][this_thread][0];
                    incy_ = &oneInt;
                    *beta_ = 0.0;
                    postProc = 1;
                }
                else if (*ch == 'T') {
                    x_ = x;
                    y_ = &y[startPt*(*incy)];
                    incy_ = incy;
                    beta_ = beta;
                    postProc = 0;
                }
            }
            
            // Compute thread contribution to matrix-vector product:
            SGEMV(ch, m_, n_, alpha, A_, lda, x_, incx, beta_, y_, incy_);
// // //         }
        
        // Gather final result of matrix-vector product:
        if (postProc == 1) {
            for (i = 0; i < lenDimNot2div; i++)
                y[i*(*incy)] = (*beta) * y[i*(*incy)] + ys[0][0][i];
            for (j = 1; j < noThreads; j++)
                for (i = 0; i < lenDimNot2div; i++)
                    y[i*(*incy)] += ys[0][j][i];
        }
    }
}

// void SGEMVn_OpenMP(char* chn, Int* m, Int* n, float* alpha, float* A, Int* lda, 
//         float* x, Int* incx, float* beta, float* y, Int* incy, Int noThreads)
// {
//     /* y = alpha*X * x + beta y */
//     
//     // This OpenMP implementation assumes m >>> n
//     #pragma omp parallel num_threads(noThreads)
//     {
//         Int this_thread = (Int) omp_get_thread_num();
//         Int myChunkSz, startPt;
//         
//         splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);
// 
//         SGEMV(chn, &myChunkSz, n, alpha, &A[startPt], lda, x,
//               incx, beta, &y[startPt*(*incy)], incy);
//     }
// }
// 
// void SGEMVt_OpenMP(char* cht, Int* m, Int* n, float* alpha, float* A, Int* lda, 
//         float* x, Int* incx, float* beta, float* y, Int* incy, vector<vector<float> >* ys, Int noThreads) 
// {
//     // TODO: Not sure if this routines works properly, but that's irrelevant since it is not used anymore.

//     /* y = alpha*X^T * x + beta y */
//     
//     // This OpenMP implementation assumes m >>> n
//     Int oneInt = 1, i, j;
//     float zero = 0.0;
//     #pragma omp parallel num_threads(noThreads)
//     {
//         Int this_thread = (Int) omp_get_thread_num();
//         Int myChunkSz, startPt;
//         
//         splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);
// 
//         SGEMV(cht, &myChunkSz, n, alpha, &A[startPt], lda, 
//                 &x[startPt*(*incx)], incx, &zero, &ys[0][this_thread][0], &oneInt);
//     }
//     for (i = 0; i < *n; i++)
//         y[i*(*incy)] = (*beta) * y[i*(*incy)] + ys[0][0][i];
//     for (j = 1; j < noThreads; j++)
//         for (i = 0; i < *n; i++)
//             y[i*(*incy)] += ys[0][j][i];
// }
// 
// void SGEMV_OpenMP(char* ch, Int* m, Int* n, float* alpha, float* A, Int* lda, 
//         float* x, Int* incx, float* beta, float* y, Int* incy, vector<vector<float> >* ys, Int noThreads)
// {
//     if (*ch == 'N')
//         SGEMVn_OpenMP(ch, m, n, alpha, A, lda, 
//             x, incx, beta, y, incy, noThreads);
//     else if (*ch == 'T')
//         SGEMVt_OpenMP(ch, m, n, alpha, A, lda, 
//             x, incx, beta, y, incy, ys, noThreads);
// }

void SSCAL_OpenMP(Int* m, float* alpha, float* x, Int* inc, Int noThreads)
{
    /* x <- alpha*x */
    
    Int threading = decideThreadingLevel1(*m, noThreads, 2);
    
    if (!threading)
        SSCAL(m, alpha, x, inc);
    else {
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            SSCAL(&myChunkSz, alpha, &x[startPt*(*inc)], inc);
// // //         }
    }
}

void SAXPY_OpenMP(Int* m, float* alpha, float* x, Int* incx, float* y, Int* incy, Int noThreads)
{
    /* y <- alpha * x + y */
    
    Int threading = decideThreadingLevel1(*m, noThreads, 2);
    
    if (!threading)
        SAXPY(m, alpha, x, incx, y, incy);
    else {
// // //         #pragma omp parallel num_threads(noThreads)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            SAXPY(&myChunkSz, alpha, &x[startPt*(*incx)], incx, &y[startPt*(*incy)], incy);
// // //         }
    }
}

float SNRM2_OpenMP(Int *m, float* x, Int* incx, Int noThreads)
{
    float x2 = 0.0;
    Int threading = decideThreadingLevel1(*m, noThreads, 2);
    
    if (!threading) {
        x2 = SNRM2(m, x, incx);
        x2 = x2*x2;
    }
    else {
// // //         #pragma omp parallel num_threads(noThreads) reduction(+:x2)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            if (myChunkSz > 0)
                x2 = SNRM2(&myChunkSz, &x[startPt*(*incx)], incx);
            else
                x2 = 0.0;
            x2 = x2*x2;
// // //         }
    }
    
    return sqrt(x2);
}

float SDOT_OpenMP(Int *m, float* x, Int* incx, float* y, Int* incy, Int noThreads)
{
    float dot = 0.0;
    Int threading = decideThreadingLevel1(*m, noThreads, 2);
    
    if (!threading)
        dot = SDOT(m, x, incx, y, incy);
    else {
// // //         #pragma omp parallel num_threads(noThreads) reduction(+:dot)
// // //         {
// // //             Int this_thread = (Int) omp_get_thread_num();
            Int this_thread = 0;
            Int myChunkSz, startPt;
            
            splitArray1D(*m, (Int) this_thread, noThreads, &myChunkSz, &startPt);

            if (myChunkSz > 0)
                dot = SDOT(&myChunkSz, &x[startPt*(*incx)], incx, &y[startPt*(*incy)], incy);
            else
                dot = 0.0;
// // //         }
    }
    
    return dot;
}

#ifdef  HAVE_MPI

double SDOTMPI(Int m, float* x, Int incx, float* y, Int incy, Int noThreads) 
{
//     float local_dot = SDOT(&m, x, &incx, y, &incy);
    float local_dot = SDOT_OpenMP(&m, x, &incx, y, &incy, noThreads);
    
    float global_dot;
    MPI_Allreduce(&local_dot, &global_dot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    return global_dot;
}

void SGEMVMPI(char ch, Int m, Int n, float alpha, float* A, Int lda, float* x, Int incx, 
        float beta, float* y, Int incy, float* ylocal, vector<vector<float> >* ylocals, Int noThreads) 
{
    // TODO: I think this function may not work if beta != 0. Same with other MPI functions here.

    /* y = alpha*X^T * x + beta y */
//     SGEMV(&cht, &m, &n, &alpha, A, &lda, x, &incx, &beta, ylocal, &incy);
    SGEMV_OpenMP(&ch, &m, &n, &alpha, A, &lda, x, &incx, &beta, ylocal, &incy, ylocals, noThreads);
    
    MPI_Allreduce(ylocal, y, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}

void SGEMMMPI(char cha, char chb, Int m, Int n, Int k, float alpha, float* A, Int lda,
              float* B, Int ldb, float beta, float* C, Int ldc, float* Clocal, Int noThreads)
{
    // TODO: Implement a threaded version of this function.
    
    /* C = alpha* opt(A) * opt(B) + beta C */
    SGEMM(&cha, &chb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, Clocal, &ldc);

    Int p = m*n;
    MPI_Allreduce(Clocal, C, p, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}

#endif

#endif
