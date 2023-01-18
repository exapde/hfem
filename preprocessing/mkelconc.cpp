#include "mex.h"
#include <math.h>
//#include <stdint.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// Written by: C. Nguyen & P. Fernandez

using namespace std;
typedef int Int;

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
    cout.precision(prec);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void print1darray(double* a, Int m)
{
    cout.precision(4);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void print2darray(double* a, Int m, Int n)
{
    cout.precision(4);
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << scientific << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3darray(double* a, Int m, Int n, Int p)
{
    cout.precision(4);
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

vector<Int> xiny(vector< vector<double> > &x, vector< vector<double> > &y, Int m, Int n)
{
// Determine if each row of x is a member of y
// If row i of x is a member of y and x(i,:) = y(j,:) then in(i) = j
// Else in(i) = -1    
    Int i, j, k;
    Int dim = x.size();
    vector<Int> in(m,-1);    
    double d;    
    for (i=0; i<m; i++) 
        for (j=0; j<n; j++) {
            d = (x[0][i]-y[0][j])*(x[0][i]-y[0][j]);
            for (k=1; k<dim; k++)
                d += (x[k][i]-y[k][j])*(x[k][i]-y[k][j]);
            if (d<1e-12) {
                in[i] = j;
                break;
            }
        }
            
    return in;
}                

vector<Int> find(const vector<Int> &x, Int a, Int compare)
{
    Int i, j = 0;
    Int n = x.size();        
    vector<Int> idx(n,0);        
    
    switch (compare) {
        case 0: // equal to
            for (i=0; i<n; i++) 
                if (x[i] == a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case 1: // greater than
            for (i=0; i<n; i++) 
                if (x[i] > a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;
        case -1: // less than
            for (i=0; i<n; i++) 
                if (x[i] < a) {
                    idx[j] = i;
                    j += 1;
                }    
            break;                 
        default:
            mexErrMsgTxt("Comparison operator not implemented.\n");     
    }
    idx.erase(idx.begin()+j,idx.end());    
    
    return idx;
}


void elementnodes(vector< vector<double> > &elemnodes, vector< vector<double> > &pv, vector<double> &philocal, Int dim, Int npe, Int nve)
{
    Int j, k, m; 
    //Int nve = round(((double) philocal.size())/npe);            
    for (j=0; j<dim; j++)
        for (k=0; k<npe; k++) {
            elemnodes[j][k] = 0.0;    
            for (m=0; m<nve; m++)
                elemnodes[j][k] += philocal[m*npe+k]*pv[j][m];
        }
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{        
    double* p = mxGetPr(prhs[0]);
    double* t = mxGetPr(prhs[1]);
    double* f = mxGetPr(prhs[2]);
    double* t2f = mxGetPr(prhs[3]);
    double* elementtype = mxGetPr(prhs[4]);   
    double* elem = mxGetPr(prhs[5]);   
    double* isEDGface = mxGetPr(prhs[6]);    
    double* phivl = mxGetPr(prhs[7]);        
    double* phifc = mxGetPr(prhs[8]);
    double* permface = mxGetPr(prhs[9]);
    double* npes = mxGetPr(prhs[10]);
    double* nves = mxGetPr(prhs[11]);
    double* nfes = mxGetPr(prhs[12]);
    double* npfs = mxGetPr(prhs[13]);
    double* nvfs = mxGetPr(prhs[14]);
    double* ndims = mxGetPr(prhs[15]);
        
    Int i, j, k, l, m, n, q, r, ei, npe, nve, nfe, npf, nvf, ne, nf, np;
    Int dim, npemax, nvemax, nfemax, npfmax, nvfmax, nelem;    
    Int e1, e2, if1, if2, ii, neg, ndof, nedg;
    
    dim = (Int) ndims[0];     
    ne = (Int) ndims[1];     
    nf = (Int) ndims[2];     
    np = (Int) ndims[3];     
    npemax  = (Int) ndims[4];
    nvemax  = (Int) ndims[5];
    nfemax = (Int) ndims[6]; 
    npfmax  = (Int) ndims[7];
    nvfmax = (Int) ndims[8]; 
    nelem  = (Int) ndims[9]; 
        
    vector< vector<double> > philocvl(nelem, vector<double>());   
    vector< vector< vector<double> > > philocfc(nelem, vector< vector<double> >());      
    vector< vector< vector<Int> > > perm(nelem, vector< vector<Int> >());      
    vector< vector<double> > pf(dim, vector<double>(nvfmax));   
    vector< vector<double> > pv(dim, vector<double>(nvemax));   
    vector< vector<double> > facenodes(dim, vector<double>(npfmax)); 
    vector< vector<double> > dgface(dim, vector<double>(npfmax)); 
    vector< vector<double> > dgnodes(dim, vector<double>(npemax));   
    vector< vector<double> > edg(dim, vector<double>(nf*npfmax));   
    vector<Int> ieg(nf*npfmax, 0);    
    vector<Int> fi(nvfmax, 0);      
    vector<Int> tof(npfmax, 0);      
    vector<Int> iedg(npfmax, 0);      
    vector<Int> in, jn, kn;        
    
//     print1darray(ndims, 10);
//     print1darray(npes, 1);
//     print1darray(nves, 1);
//     print1darray(nfes, 1);
//     print1darray(nvfs, 4);
//     print1darray(npfs, 4);
//     mexPrintf("%d\n",nelem);
    //mexErrMsgTxt("Comparison operator not implemented.\n");
    
    m = 0; r = 0; q = 0;
    for (n=0; n<nelem; n++) {
        //i = (Int) elem[n];        
        npe = (Int) npes[n];
        nve = (Int) nves[n];        
        philocvl[n].resize(npe*nve);
        for (j=0; j<nve; j++)
            for (k=0; k<npe; k++)
                philocvl[n][j*npe+k] = phivl[m+j*npe+k];
        m += npe*nve;                                 
        
        nfe = (Int) nfes[n];
        philocfc[n].resize(nfe);
        perm[n].resize(nfe);
        for (l=0; l<nfe; l++) {
            npf = (Int) npfs[l*nelem+n];
            nvf = (Int) nvfs[l*nelem+n];
            philocfc[n][l].resize(npf*nvf);            
            for (j=0; j<nvf; j++)
                for (k=0; k<npf; k++) 
                    philocfc[n][l][j*npf+k] = phifc[r+j*npf+k];                            
            r += npf*nvf;                                    
            
            perm[n][l].resize(npf);
            for (k=0; k<npf; k++) 
                perm[n][l][k] = (Int) permface[q+k];                            
            q += npf;             
        }                        
    }         
    
    
    //mexErrMsgTxt("Comparison operator not implemented.\n");
                       
    mwSize szQ[3]; 
    szQ[0]=npfmax; szQ[1]=nfemax; szQ[2] = ne; 
    plhs[0]=mxCreateNumericArray(3,szQ,mxDOUBLE_CLASS,mxREAL);                  
    double *elcon;
    elcon = mxGetPr(plhs[0]);        
    for (i=0; i<ne*nfemax*npfmax; i++)
        elcon[i] = -1;
    
    mwSize szR[2]; 
    szR[0]=npfmax; szR[1]=nf;
    plhs[1]=mxCreateNumericArray(2,szR,mxDOUBLE_CLASS,mxREAL);              
    double *facecon;
    facecon = mxGetPr(plhs[1]);        
    for (i=0; i<nf*npfmax; i++)
        facecon[i] = -1;
        
    Int nmax = (1+nvfmax+2);     
    //print1darray(f, nmax);
    //mexErrMsgTxt("Comparison operator not implemented.\n");
    //mexPrintf("%d\n",nmax);    
    ndof = 0;
    nedg = 0;
    for (i=0; i<nf; i++) {      
        nvf = (Int) f[i*nmax]; //nvf = elements[ei].face[if1].size();
        for (j=0; j<nvf; j++)
            fi[j] = (Int) f[i*nmax+j+1];
        e1 = (Int) f[i*nmax+nvfmax+1]; 
        e2 = (Int) f[i*nmax+nvfmax+2];        
        ei = (Int) elementtype[e1];            
        
//         print1darray(f, nmax);
//         printiarray(fi);        
//         mexPrintf("%d  %d   %d   %d\n",nvf,e1,e2,ei);    
//         mexErrMsgTxt("Comparison operator not implemented.\n");
        
        // index of ei in elem array
        for (j=0; j<nelem; j++)
            if ((elem[j]-ei)==0)
                break;
        ii = j;
        
        npe = (Int) npes[ii];             
        nfe = (Int) nfes[ii]; 
        nve = (Int) nves[ii];     
        
        // index of face i on element e1
        for (j=0; j<nfe; j++)
            if ((t2f[e1*nfemax+j]-i)==0)
                break;
        if1 = j;          
        npf = (Int) npfs[if1*nelem+ii];
        //nvf = (Int) nvfs[if1*nelem+ii];
                            
        // coordinates of the face i         
        for (m=0; m<nvf; m++)             
            for (n=0; n<dim; n++)
                pf[n][m] = p[n*np+fi[m]];        
        
        // nodal points on face i
        for (j=0; j<dim; j++)
            for (k=0; k<npf; k++) {
                facenodes[j][k] = 0;    
                for (m=0; m<nvf; m++)
                    facenodes[j][k] += philocfc[ii][if1][m*npf+k]*pf[j][m];
            }
        
//         cout.precision(4);
//         for (j=0; j<dim; j++) {
//             for (k=0; k<nvf; k++) 
//                 cout<< scientific <<pf[j][k]<< "   ";    
//             cout << endl;
//         }
//         for (m=0; m<nvf; m++) {
//             for (k=0; k<npf; k++) 
//                 cout<< scientific <<philocfc[ii][if1][m*npf+k]<< "   ";    
//             cout << endl;
//         }
//         for (j=0; j<dim; j++) {
//             for (k=0; k<npf; k++) 
//                 cout<< scientific <<facenodes[j][k]<< "   ";    
//             cout << endl;
//         }
//          mexErrMsgTxt("Comparison operator not implemented.\n");
                                
        if (isEDGface[i]==0) { // hdg face
            // dof numbering on face i
            for (k=0; k<npf; k++)
                tof[k] = ndof + k;
            ndof += npf;  // number of degrees of freedom 
        }        
        else { // edg face
            if (nedg==0) {
                for (j=0; j<dim; j++)
                    for (k=0; k<npf; k++)
                        edg[j][k] = facenodes[j][k]; // edg nodes on faces
                for (k=0; k<npf; k++) {
                    tof[k] = ndof + k;  // dof numbering on face i    
                    ieg[k] = tof[k];    // dof numbering of edg nodes
                }
                ndof += npf; 
                nedg += npf;
            }
            else {
                in = xiny(facenodes, edg, npf, nedg);  // find which rows of pf are in edg 
                jn = find(in,-1,0);  // new edg nodes have in==-1       
                kn = find(in,-1,1);  // old edg nodes have in>-1                  
                neg = jn.size();         // number of new edg nodes                         
                
                for (j=0; j<dim; j++)
                    for (k=0; k<neg; k++)
                        edg[j][nedg+k] = facenodes[j][jn[k]]; // update edg with new edg nodes            
                
                for (k=0; k<neg; k++)
                    tof[jn[k]] = ndof+k; // tof for new edge nodes            
                for (k=0; k<kn.size(); k++)
                    tof[kn[k]] = ieg[in[kn[k]]]; // tof for old edge nodes                          
//                 for (k=0; k<npf; k++) 
//                     if (in[k]>-1)                        
//                         tof[k] = ieg[in[k]];                
                                
                for (k=0; k<neg; k++) 
                    ieg[nedg+k] = ndof+k; // update ieg with dof numbering of new edge nodes                 
                
                ndof += neg; // update ndof                     
                nedg += neg; // update nedg                                    
            }                        
        }        
        
        // dgnodes on element e1 
        for (m=0; m<nve; m++) {      
            k = (Int) t[e1*nvemax+m];
            for (n=0; n<dim; n++)                 
                pv[n][m] = p[n*np+k];                    
        }
                
        elementnodes(dgnodes, pv, philocvl[ii], dim, npe, nve);
        
//         cout.precision(4);
//         for (j=0; j<dim; j++) {
//             for (k=0; k<nve; k++) 
//                 cout<< scientific <<pv[j][k]<< "   ";    
//             cout << endl;
//         }
//         for (m=0; m<nve; m++) {
//             for (k=0; k<npe; k++) 
//                 cout<< scientific <<philocvl[ii][m*npe+k]<< "   ";    
//             cout << endl;
//         }        
//         for (j=0; j<dim; j++) {
//             for (k=0; k<npe; k++) 
//                 cout<< scientific <<dgnodes[j][k]<< "   ";    
//             cout << endl;
//         }
//          mexErrMsgTxt("Comparison operator not implemented.\n");
                                        
        // dg nodes on face i from element e1
        for (j=0; j<dim; j++)
            for (k=0; k<npf; k++)                
                dgface[j][k] = dgnodes[j][perm[ii][if1][k]];

        in = xiny(dgface, facenodes, npf, npf); // match facenodes to dg1   
        
//         for (j=0; j<dim; j++) {
//             for (k=0; k<npf; k++) 
//                 cout<< scientific <<dgface[j][k]<< "   ";    
//             cout << endl;
//         }        
//         for (j=0; j<dim; j++) {
//             for (k=0; k<npe; k++) 
//                 cout<< scientific <<dgnodes[j][k]<< "   ";    
//             cout << endl;
//         }
//         for (j=0; j<dim; j++) {
//             for (k=0; k<npf; k++) 
//                 cout<< scientific <<facenodes[j][k]<< "   ";    
//             cout << endl;
//         }               
//         printiarray(in);
//         mexErrMsgTxt("Comparison operator not implemented.\n");
        
        // assign dof numbering of face i to elcon from element e1                        
        for (k=0; k<npf; k++) {                
            elcon[npfmax*nfemax*e1+npfmax*if1+k] = tof[in[k]];
            facecon[npfmax*i+k] = tof[k];
        }
        
        if (e2>-1) {               
            ei = elementtype[e2];
            // index of ei in elem array
            for (j=0; j<nelem; j++)
                if ((elem[j]-ei)==0)
                    break;
            ii = j;            
            npe = (Int) npes[ii];             
            nfe = (Int) nfes[ii]; 
            nve = (Int) nves[ii];     
                        
            for (j=0; j<nfe; j++)
                if ((t2f[e2*nfemax+j]-i)==0)
                    break;
            if2 = j;  // location of face i on element e2
            npf = (Int) npfs[if2*nelem+ii];
            
            // dgnodes on element e2 
            for (m=0; m<nve; m++) {            
                k = (Int) t[e2*nvemax+m];
                for (n=0; n<dim; n++)
                    pv[n][m] = p[n*np+k];        
            }
            elementnodes(dgnodes, pv, philocvl[ii], dim, npe, nve);
            
            // dg nodes on face i from element e2
            for (j=0; j<dim; j++)
                for (k=0; k<npf; k++)                
                    dgface[j][k] = dgnodes[j][perm[ii][if2][k]];            
            in = xiny(dgface, facenodes, npf, npf); // match facenodes to dg1    
            // assign dof numbering of face i to elcon from element e1                        
            for (k=0; k<npf; k++)                 
                elcon[npfmax*nfemax*e2+npfmax*if2+k] = tof[in[k]];            
        }                                                                    
    }                    
    
    if (nlhs>2) {
        mwSize szW[2]; 
        szW[0]=nedg; szW[1]=dim;
        plhs[2]=mxCreateNumericArray(2,szW,mxDOUBLE_CLASS,mxREAL);              
        double *edgnodes = mxGetPr(plhs[2]);        
        for (j=0; j<dim; j++)
            for (i=0; i<nedg; i++) 
                edgnodes[j*nedg+i] = edg[j][i];
    }
    
}
