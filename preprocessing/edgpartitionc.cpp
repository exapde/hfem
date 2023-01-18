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

vector<Int> uniqueiarray(vector<Int> & b)
{	
    // input: integer array b
    // output: sorted integer array a without duplicate elements
                     
    // make a new copy of b
    vector<Int> a = b;
            
	// First Sort the given range to bring duplicate
	// elements at consecutive positions
	sort(a.begin(), a.end());
  
	vector<Int>::iterator newEnd;
 
	// Override duplicate elements
	newEnd = unique(a.begin(), a.end());
 
    // remove duplicate elements
	a.erase(newEnd, a.end());
    
    return a;
}

void uniqueiarray(vector<Int> & a, vector<Int> & b)
{	
    // input: integer array b
    // output: sorted integer array a without duplicate elements
    
    // make a new copy of b
    a = b;
            
	// First Sort the given range to bring duplicate
	// elements at consecutive positions
	sort(a.begin(), a.end());
  
	vector<Int>::iterator newEnd;
 
	// Override duplicate elements
	newEnd = unique(a.begin(), a.end());
 
    // remove duplicate elements
	a.erase(newEnd, a.end());
}

vector<Int> setdifference(vector<Int> v1, vector<Int> v2)
{

    vector<Int> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
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

vector<Int> find(double* x, Int a,  Int n)
{
    Int i, j = 0;    
    double b = a;
    vector<Int> idx(n,0);        
    for (i=0; i<n; i++) 
        if (x[i] == b) {
            idx[j] = i;
            j += 1;
        }    
    idx.erase(idx.begin()+j,idx.end());    
    
    return idx;
}


vector<Int> findcol(const vector<Int> &x, Int a, Int m, Int n)
{
    Int i, k, j = 0;      
    
    vector<Int> idx(n,-1);        
    for (i=0; i<n; i++) 
        for (k=0; k<m; k++)
            if (x[i*m+k] == a) {
                idx[j] = i;
                j += 1;
                break;
            }    
    idx.erase(idx.begin()+j,idx.end());    
    
    return idx;
}

    
vector<Int> iarrayatindex(vector<Int> &a, vector<Int> &ind)
{
    Int i, n = ind.size();
    vector<Int> b(n,0); 
    for (i = 0; i<n; i++)
        b[i] = a[ind[i]];    
    return b;
}

vector<Int> darrayatindex(double *a, vector<Int> &ind)
{
    Int i, n = ind.size();
    vector<Int> b(n,0); 
    for (i = 0; i<n; i++)
        b[i] = (Int) a[ind[i]];    
    return b;
}

vector<Int> iarray2datindex(vector<Int> &a, vector<Int> &ind, Int m)
{
    Int i, j, n = ind.size();
    vector<Int> b(m*n,0); 
    for (i = 0; i<n; i++)
        for (j = 0; j<m; j++)
            b[i*m+j] = a[ind[i]*m+j];    
    return b;
}

vector<Int> darray2datindex(double* a, vector<Int> &ind, Int m)
{
    Int i, j, n = ind.size();
    vector<Int> b(m*n,0); 
    for (i = 0; i<n; i++)
        for (j = 0; j<m; j++)
            b[i*m+j] = (Int) a[ind[i]*m+j];    
    return b;
}


vector<double> darrayatindex(vector<double> &a, vector<Int> &ind)
{
    Int n = ind.size();
    vector<double> b(n,0); 
    for (Int i = 0; i<n; i++)
        b[i] = a[ind[i]];    
    return b;
}

vector<Int> mkintent(vector<Int> &facecon, vector<Int> &face2cpu, Int nin, Int npfmax, Int nf, Int my_rank)
{            
    Int i, j, m, icpu, imax;
    vector<Int> ci, cpus, a, b, ent;
    ci.resize(nin,0);
    for (i=0; i<nin; i++)
        ci[i] = i;
    a = iarray2datindex(facecon,ci,npfmax);
    ent = uniqueiarray(a);    
    if (ent[0]<0)
        ent.erase(ent.begin());                
    Int nent = ent.size();
    
    //cout<<nent<<endl;    
    vector<Int> intent(nent,-1);
    for (i=0; i<nent; i++) {
        ci = findcol(facecon,ent[i],npfmax,nf);
        cpus = iarrayatindex(face2cpu,ci);
        uniqueiarray(a, cpus);    
        b = a;
        m = b.size();
        for (j=0; j<m; j++)
            b[j] = find(cpus,a[j],0).size();
        imax = distance(b.begin(), max_element(b.begin(), b.end()));
        icpu = a[imax];
        if (icpu==my_rank)
            intent[i] = 1;
    }
    ci = find(intent,1,0);
    intent = iarrayatindex(ent,ci);  
    
    return intent;
}

mxArray* to_mxArray(const vector<Int> &a) 
{
   Int n = a.size();   
   mwSize sz1[1]; sz1[0]=n; 
   mxArray *b = mxCreateNumericArray(1,sz1,mxDOUBLE_CLASS,mxREAL);                         
   double *data = mxGetPr(b);
   for (Int i = 0; i < n; i++)       
         data[i] = a[i];         
   return b;
}

vector<Int> from_mxArray(const mxArray *a) 
{      
   mwSize nelem = mxGetNumberOfElements(a); 
   double* c = mxGetPr(a);
   Int n = (Int) nelem;   
   vector<Int> b(n,-1);      
   for (Int i = 0; i < n; i++)       
         b[i] = (Int) c[i];         
   return b;
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{            
    const mxArray *extintelem;
    extintelem = prhs[0];    
    double* elcon = mxGetPr(prhs[1]);
    double* facecon = mxGetPr(prhs[2]);
    double* t2f = mxGetPr(prhs[3]);
    double* face2cpu = mxGetPr(prhs[4]);
    double* ndims = mxGetPr(prhs[5]);
        
    Int i, j, ne, nf, nent, nfemax, npfmax, nproc;         
    ne = (Int) ndims[0];       
    nf = (Int) ndims[1];       
    nent = (Int) ndims[2];         
    nfemax = (Int) ndims[3];
    npfmax = (Int) ndims[4];
    nproc = (Int) ndims[5];
                        
    //cout<<ne<<" "<<nf<<" "<<nent<<" "<<nfemax<<" "<<npfmax<<" "<<nproc<<endl;
    //mexErrMsgTxt("here.\n");           
    
    mwSize sz0[1]; sz0[0]=nent; 
    plhs[0]=mxCreateNumericArray(1,sz0,mxDOUBLE_CLASS,mxREAL);                  
    double *ent2cpu = mxGetPr(plhs[0]);            

    mxArray *out1 = mxCreateCellMatrix((mwSize)nproc,1);
        
    mwSize sz2[2]; sz2[0]=nproc; sz2[1]=2; 
    plhs[2]=mxCreateNumericArray(2,sz2,mxDOUBLE_CLASS,mxREAL);                  
    double *out2 = mxGetPr(plhs[2]);        
    
    vector<Int> elem, extintent, extintface, temp1, temp2;   
    Int nintface, nextface, nintent, nextent;
    for (i=0; i<nproc; i++) {
        elem = from_mxArray(mxGetCell(extintelem,(mwSize) i));
        
        temp1 = darray2datindex(t2f,elem,nfemax);
        temp2 = uniqueiarray(temp1); // all faces on subdomain i    
        if (temp2[0]<0)
            temp2.erase(temp2.begin());                   
        //interior faces       
        extintface = find(face2cpu,i,nf); 
        nintface = extintface.size(); 
        temp1 = setdifference(temp2, extintface); // exterior faces on subdomain i      
        nextface = temp1.size();
        extintface.insert(extintface.end(), temp1.begin(), temp1.end());                                     

        temp1 = darray2datindex(elcon,elem,nfemax*npfmax);        
        temp2 = uniqueiarray(temp1); // all entities on subdomain i     
        if (temp2[0]<0)
            temp2.erase(temp2.begin());   
        
        // interior entities on subdomain i     
        elem = darray2datindex(facecon,extintface,npfmax);
        temp1 = darrayatindex(face2cpu,extintface);
        extintent = mkintent(elem, temp1, nintface, npfmax, nintface+nextface, i);    
        nintent = extintent.size();         
        temp1 = setdifference(temp2, extintent);  // exterior entities on subdomain i     
        nextent = temp1.size();
        extintent.insert(extintent.end(), temp1.begin(), temp1.end());                  
        mxSetCell(out1,i,to_mxArray(extintent));        
        out2[i] = nintent;
        out2[nproc+i] = nextent;                   
        for (j=0; j<nintent; j++)
            ent2cpu[extintent[j]] = i;                    
    }                      
    plhs[1] = out1;    
}
