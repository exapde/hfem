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

vector<Int> darray2datindex(double* a, vector<Int> &ind, Int m)
{
    Int i, j, n = ind.size();
    vector<Int> b(m*n,0); 
    for (i = 0; i<n; i++)
        for (j = 0; j<m; j++)
            b[i*m+j] = (Int) a[ind[i]*m+j];    
    return b;
}

void elcon2entcon(vector<Int> &rowent2elem, vector<Int> &colent2elem, double* elem2ent, Int np, Int ne)
{
    Int i, j, ind;
    //Int ndof =  *max_element(elem2ent.begin(), elem2ent.end());
    double emax = 0;
    for (i=0; i<np*ne; i++)
        if (elem2ent[i]>emax)
            emax = elem2ent[i];
    Int ndof = (Int) (emax+1);
    
    // store number of neighboring elements for each entity    
    rowent2elem.resize(ndof,0);
//     for (i=0; i<ndof+1; i++)
//         rowent2elem[i] = 0;
    
    vector<Int> elc, k;
    elc.resize(np);
    for (i=0; i<ne; i++) { // for each element i
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = (Int) elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = k[j];  // get entity index
            if (ind>=0)
                rowent2elem[ind] += 1; // increase the number of elements by 1
        }            
    }     
    // cummulative sum
    partial_sum(rowent2elem.begin(), rowent2elem.end(), rowent2elem.begin());
    rowent2elem.insert(rowent2elem.begin(),0);
    
    //printiarray(rowent2elem);    
    //mexErrMsgTxt("Here.\n");    
    
    // store neighboring-element indices for each entity
    colent2elem.resize(rowent2elem.back(),0);
//     for (i=0; i<rowent2elem.back(); i++)
//         colent2elem[i] = 0;
    
    vector<Int> inc(ndof,0);
    for (i=0; i<ne; i++) {
        for (j=0;j<np;j++) // entities on element i  
            elc[j] = (Int) elem2ent[i*np+j];
        uniqueiarray(k,elc); // remove duplicate entities on element i      
        for (j=0; j<k.size(); j++) {
            ind = k[j];  // get entity index
            if (ind>=0) {
                colent2elem[rowent2elem[ind]+inc[ind]] = i;
                inc[ind] += 1; // pointer to the next element
            }
        }                    
    }       
}

void ent2elem(vector<Int> &elems, vector<Int> &rowent2elem, vector<Int> &colent2elem, vector<Int> &ents)
{        
    // sort and remove duplications
    vector<Int>  ent;
    uniqueiarray(ent, ents);    
    if (ent[0]<0)
        ent.erase(ent.begin());            
    
    Int nent = ent.size();
    Int i, j, k, ei, ni, ri, nt = 0;
    for (i=0; i<nent; i++)
        nt += rowent2elem[ent[i]+1]-rowent2elem[ent[i]];
    
    vector<Int> elem(nt,-1);
    j = 0;
    for (i=0; i<nent; i++) {
        ei = ent[i];
        ni = rowent2elem[ei+1]-rowent2elem[ei];
        for (k=0; k<ni; k++) {
            ri = rowent2elem[ei]+k;
            elem[j] = colent2elem[ri];
            j += 1;
        }
    }    
    uniqueiarray(elems,elem);
}

void mkextintelem(vector<Int> &extintelem, vector<Int> &extintelempts, vector<Int> &re, vector<Int> &ce, 
        double* elem2cpu, double* elem2ent, Int np, Int ne, Int overlappinglevel, Int my_rank)
{            
    Int j, m, k, nelem;
    //extintelem = find(elem2cpu,my_rank,0);    
    extintelem = find(elem2cpu,my_rank,ne);    
    extintelempts[0] = extintelem.size();
    vector<Int> elems = extintelem;
    
    vector<Int> ents;
    for (j=0; j<overlappinglevel; j++) {
        nelem = elems.size();
        ents.resize(np*nelem);
        for (m=0; m<nelem; m++)
            for (k=0; k<np; k++)
                ents[m*np+k] = (Int) elem2ent[elems[m]*np+k];
        ent2elem(elems, re, ce, ents);
    }
    
    vector<Int> extelem = setdifference(elems, extintelem);   
    extintelempts[1] = extelem.size();    
    extintelem.insert(extintelem.end(), extelem.begin(), extelem.end());                  
}
    
mxArray* to_mxArray(const vector<Int> &a) {
   Int n = a.size();   
   mwSize sz1[1]; sz1[0]=n; 
   mxArray *b = mxCreateNumericArray(1,sz1,mxDOUBLE_CLASS,mxREAL);                         
   double *data = mxGetPr(b);
   for (Int i = 0; i < n; i++)       
         data[i] = a[i];         
   return b;
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{            
    double* t = mxGetPr(prhs[0]);
    double* t2f = mxGetPr(prhs[1]);
    double* elem2cpu = mxGetPr(prhs[2]);        
    double* face2cpu = mxGetPr(prhs[3]);
    double* ndims = mxGetPr(prhs[4]);
    
    Int i, j, k, m, n, ei, ne, np, nlevel, nproc, nf, nfemax;         
    np = (Int) ndims[0];     
    ne = (Int) ndims[1];   
    nf = (Int) ndims[2];
    nfemax = (Int) ndims[3];
    nlevel = (Int) ndims[4];
    nproc = (Int) ndims[5];    
    
    mxArray *out0 = mxCreateCellMatrix((mwSize)nproc,1);
    mxArray *out2 = mxCreateCellMatrix((mwSize)nproc,1);
    mxArray *out4 = mxCreateCellMatrix((mwSize)nproc,1);
    
    mwSize sz2[2]; sz2[0]=nproc; sz2[1]=2; 
    plhs[1]=mxCreateNumericArray(2,sz2,mxDOUBLE_CLASS,mxREAL);                  
    double *out1 = mxGetPr(plhs[1]);        
    plhs[3]=mxCreateNumericArray(2,sz2,mxDOUBLE_CLASS,mxREAL);                  
    double *out3 = mxGetPr(plhs[3]);        
    
    vector<Int> re, ce;    
    elcon2entcon(re, ce, t, np, ne);        
    //printiarray(re);
    //mexErrMsgTxt("Here.\n");        
    vector<Int> extintelem, extintface, temp1, temp2;                    
    vector<Int> extintelempts(2,0);    
    Int nintface, nextface;
    for (i=0; i<nproc; i++) {
        mkextintelem(temp1, extintelempts, re, ce, elem2cpu, t, np, ne, nlevel+1, i);                             
        mkextintelem(extintelem, extintelempts, re, ce, elem2cpu, t, np, ne, nlevel, i);      
        temp2 = setdifference(temp1,extintelem);
        mxSetCell(out4,i,to_mxArray(temp2));        
        mxSetCell(out0,i,to_mxArray(extintelem));
        out1[i] = extintelempts[0];
        out1[nproc+i] = extintelempts[1];  
        
        temp1 = darray2datindex(t2f,extintelem,nfemax);
        temp2 = uniqueiarray(temp1); // all faces on subdomain i    
        if (temp2[0]<0)
            temp2.erase(temp2.begin());                   
        //interior faces       
        extintface = find(face2cpu,i,nf); 
        nintface = extintface.size(); 
        temp1 = setdifference(temp2, extintface); // exterior faces on subdomain i      
        nextface = temp1.size();
        extintface.insert(extintface.end(), temp1.begin(), temp1.end());                  
        mxSetCell(out2,i,to_mxArray(extintface));
        out3[i] = nintface;
        out3[nproc+i] = nextface;                   
    }          
    plhs[0] = out0;
    plhs[2] = out2;    
    plhs[4] = out4;    
}
