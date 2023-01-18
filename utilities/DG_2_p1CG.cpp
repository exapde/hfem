#ifndef __DG2P1CG
#define __DG2P1CG
 
// Written by C. Nguyen and P. Fernandez
 
void DG_2_p1CG_2d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
{
    // avgMethod:   0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
 
    Int i, j, currentElem, currentVertex, currentNode;
    Int nd = ndims[0];
    Int nve = ndims[3];
    Int ne = ndims[5];
//    Int nv = ndims[7];
    Int npv = ndims[9];
    Int porder = app.porder;
    Int elemtype = getElemtype(nd, nve);
     
        Int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element
 
    // Get number of total vertices in the mesh (vertices on periodic boundaries are counted twice, i.e., separately)
    Int numTotalVertices = 0;
    for (i = 0; i < ne*nve; i++) {
        if (mesh.t[i] > numTotalVertices)
            numTotalVertices = mesh.t[i];
    }
    numTotalVertices += 1;
//    if (nv != numTotalVertices)
//        error("There are ghost vertices in the mesh.\n");
     
    // Creation of a matrix that contains the vertex-to-element connectivity:
    Int *numElementsAdj2Vertex = new Int[numTotalVertices];
    Int *vertex2elem = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    Int *vertex2node = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    double *fieldVertices = new double[numTotalVertices];
    Int *activeElement = new Int[ne];
    for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++)
        numElementsAdj2Vertex[currentVertex] = 0;
    for (currentElem = 0; currentElem < ne; currentElem++) {
        activeElement[currentElem] = 1;
        for (j = 0; j < nve; j++) {
            currentVertex = mesh.t[j*ne+currentElem];
            getNodalIndexOfVertex(&currentNode, j, nd, elemtype, porder);
            numElementsAdj2Vertex[currentVertex] += 1;
            if (numElementsAdj2Vertex[currentVertex] > 50)
                error("Over 50 elements are connected to a vertex.\n");
            vertex2elem[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentElem;
            vertex2node[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentNode;
        }
        for (j = 0; j < npv; j++) {
            if (DGfield[currentElem*npv + j] <= 0.0) {
                activeElement[currentElem] = 0;
                break;
            }
        }
    }
     
    double *avgFieldElem = new double[ne];
    if (avgMethod == 0) {      // 0: Vertex value given by David's reconstruction
        // Compute average field in each element:
        for (currentElem = 0; currentElem < ne; currentElem++) {
            avgFieldElem[currentElem] = 0.0;
            for (j = 0; j < npv; j++)
                avgFieldElem[currentElem] += DGfield[currentElem*npv + j];     // TODO: This is only an approximation to the true average. Some weighting is required for the exact average.
            avgFieldElem[currentElem] /= ((double) npv);
        }
 
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
            fieldVertices[currentVertex] = 0.0;
            for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
                currentElem = vertex2elem[currentVertex + j*numTotalVertices];
                fieldVertices[currentVertex] += avgFieldElem[currentElem];
            }
            fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
        }
    }
    else if (avgMethod == 1) {     // 1: Vertex value is average of nodes on vertex at different elements
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
            fieldVertices[currentVertex] = 0.0;
            for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
                currentElem = vertex2elem[currentVertex + j*numTotalVertices];
                currentNode = vertex2node[currentVertex + j*numTotalVertices];
                fieldVertices[currentVertex] += DGfield[currentElem*npv + currentNode];
            }
            fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
        }
    }
    else
        error("Averaging method not recognized in DG_2_p1CG_2d.\n");
     
    // Compute p1CG field at each DG node:
    double xi1, xi2, a, b, c;
    double f1, f2, f3, f4, f12, f43;
    if (elemtype == 0) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[0*ne+currentElem]];
            f2 = fieldVertices[mesh.t[1*ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
 
            a = -f1+f2;
            b = -f1+f3;
            c = f1;
 
            for (j = 0; j < npv; j++) {
                xi1 = master.plocvl[0*npv+j];
                xi2 = master.plocvl[1*npv+j];
 
                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
    else if (elemtype == 1) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[0*ne+currentElem]];
            f2 = fieldVertices[mesh.t[1*ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];
 
            for (j = 0; j < npv; j++) {
                xi1 = master.plocvl[0*npv+j];
                xi2 = master.plocvl[1*npv+j];
 
                f12 = (1.0-xi1)*f1 + xi1*f2;
                f43 = (1.0-xi1)*f4 + xi1*f3;
 
                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = (1.0-xi2)*f12 + xi2*f43;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
 
    delete[] numElementsAdj2Vertex; delete[] vertex2elem; delete[] vertex2node; delete[] avgFieldElem; delete[] fieldVertices; delete[] activeElement;
}
 
void DG_2_p1CG_3d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
{
    // avgMethod:   0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
 
    Int i, j, currentElem, currentVertex, currentNode;
    Int nd = ndims[0];
    Int nve = ndims[3];
    Int ne = ndims[5];
//    Int nv = ndims[7];
    Int npv = ndims[9];
    Int porder = app.porder;
    Int elemtype = getElemtype(nd, nve);
 
    Int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element
 
    // Get number of total vertices in the mesh (vertices on periodic boundaries are counted twice, i.e., separately)
    Int numTotalVertices = 0;
    for (i = 0; i < ne*nve; i++) {
        if (mesh.t[i] > numTotalVertices) {
            numTotalVertices = mesh.t[i];
        }
    }
    numTotalVertices += 1;
//    if (nv != numTotalVertices) {
//        error("There are ghost vertices in the mesh.\n");
//    }
 
    // Creation of a matrix that contains the vertex-to-element connectivity:
    Int *numElementsAdj2Vertex = new Int[numTotalVertices];
    Int *vertex2elem = new Int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
    Int *vertex2node = new Int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
    double *fieldVertices = new double[numTotalVertices];
    Int *activeElement = new Int[ne];
    for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
        numElementsAdj2Vertex[currentVertex] = 0;
    }
    for (currentElem = 0; currentElem < ne; currentElem++) {
        activeElement[currentElem] = 1;
        for (j = 0; j < nve; j++) {
            currentVertex = mesh.t[j*ne+currentElem];
            getNodalIndexOfVertex(&currentNode, j, nd, elemtype, porder);
            numElementsAdj2Vertex[currentVertex] += 1;
            if (numElementsAdj2Vertex[currentVertex] > 100) {
                error("Over 100 elements are connected to a vertex.\n");
            }
            vertex2elem[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentElem;
            vertex2node[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentNode;
        }
        for (j = 0; j < npv; j++) {
            if (DGfield[currentElem*npv + j] <= 0.0) {
                activeElement[currentElem] = 0;
                break;
            }
        }
    }
 
    double *avgFieldElem = new double[ne];
    if (avgMethod == 0) {      // 0: Vertex value given by David's reconstruction
        // Compute average field in each element:
        for (currentElem = 0; currentElem < ne; currentElem++) {
            avgFieldElem[currentElem] = 0.0;
            for (j = 0; j < npv; j++) {
                avgFieldElem[currentElem] += DGfield[currentElem*npv + j];     // TODO: This is only an approximation to the true average. Some weighting is required for the exact average.
            }
            avgFieldElem[currentElem] = avgFieldElem[currentElem] / (double) npv;
        }
 
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
            fieldVertices[currentVertex] = 0.0;
            for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
                currentElem = vertex2elem[currentVertex + j*numTotalVertices];
                fieldVertices[currentVertex] += avgFieldElem[currentElem];
            }
            fieldVertices[currentVertex] = fieldVertices[currentVertex] / max((double) numElementsAdj2Vertex[currentVertex],1.0);
        }
    }
    else if (avgMethod == 1) {     // 1: Vertex value is average of nodes on vertex at different elements
        // Compute field value at each vertex (average of the average field in the adjacent elements):
        for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
            fieldVertices[currentVertex] = 0.0;
            for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
                currentElem = vertex2elem[currentVertex + j*numTotalVertices];
                currentNode = vertex2node[currentVertex + j*numTotalVertices];
                fieldVertices[currentVertex] += DGfield[currentElem*npv + currentNode];
            }
            fieldVertices[currentVertex] = fieldVertices[currentVertex] / max((double) numElementsAdj2Vertex[currentVertex],1.0);
        }
    }
    else
        error("Averaging method not recognized in DG_2_p1CG_3d.\n");
 
    // Compute p1CG field   at each DG node:
    double xi1, xi2, xi3, a, b, c, d;
    double f1, f2, f3, f4, f5, f6, f7, f8, f12, f43, f56, f87, f1234, f5678;
    if (elemtype == 0) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[0*ne+currentElem]];
            f2 = fieldVertices[mesh.t[1*ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];
 
            a = -f1+f2;
            b = -f1+f3;
            c = -f1+f4;
            d = f1;
 
            for (j = 0; j < npv; j++) {
                xi1 = master.plocvl[0*npv+j];
                xi2 = master.plocvl[1*npv+j];
                xi3 = master.plocvl[2*npv+j];
 
                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c*xi3 + d;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
    else if (elemtype == 1) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[0*ne+currentElem]];
            f2 = fieldVertices[mesh.t[1*ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];
            f5 = fieldVertices[mesh.t[4*ne+currentElem]];
            f6 = fieldVertices[mesh.t[5*ne+currentElem]];
            f7 = fieldVertices[mesh.t[6*ne+currentElem]];
            f8 = fieldVertices[mesh.t[7*ne+currentElem]];
 
            for (j = 0; j < npv; j++) {
                xi1 = master.plocvl[0*npv+j];
                xi2 = master.plocvl[1*npv+j];
                xi3 = master.plocvl[2*npv+j];
 
                f12 = (1.0-xi1)*f1 + xi1*f2;
                f43 = (1.0-xi1)*f4 + xi1*f3;
                f1234 = (1.0-xi2)*f12 + xi2*f43;
 
                f56 = (1.0-xi1)*f5 + xi1*f6;
                f87 = (1.0-xi1)*f8 + xi1*f7;
                f5678 = (1.0-xi2)*f56 + xi2*f87;
 
                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = (1.0-xi3)*f1234 + xi3*f5678;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
     
    delete[] numElementsAdj2Vertex; delete[] vertex2elem; delete[] vertex2node; delete[] avgFieldElem; delete[] fieldVertices; delete[] activeElement;
}
 
void DG_2_p1CG(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
{
    Int nd = ndims[0];
 
    if (nd == 2)
        DG_2_p1CG_2d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
    else if (nd == 3)
        DG_2_p1CG_3d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
    else
        error("Invalid number of dimensions.\n");
}
 
#endif


// #ifndef __DG2P1CG
// #define __DG2P1CG
// 
// // Written by C. Nguyen and P. Fernandez
// 
// void DG_2_p1CG_2d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
// {
//     // avgMethod:   0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
// 
//     Int i, j, currentElem, currentVertex, currentNode;
//     Int nd = ndims[0];
//     Int nve = ndims[3];
//     Int ne = ndims[5];
//     Int npv = ndims[9];
//     Int nqvR = master.nqvR;
//     Int porder = app.porder;
//     Int elemtype = getElemtype(nd, nve);
//     double neighMearures;
//     
//     Int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element
//     
//     // Get number of total vertices in the mesh (vertices on periodic boundaries are counted twice, i.e., separately)
//     Int numTotalVertices = 0;
//     for (i = 0; i < ne*nve; i++) {
//         if (mesh.t[i] > numTotalVertices)
//             numTotalVertices = mesh.t[i];
//     }
//     numTotalVertices += 1;
// //    if (nv != numTotalVertices)
// //        error("There are ghost vertices in the mesh.\n");
//     
//     // Creation of a matrix that contains the vertex-to-element connectivity:
//     Int *numElementsAdj2Vertex = new Int[numTotalVertices];
//     Int *vertex2elem = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
//     Int *vertex2node = new Int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
//     double *fieldVertices = new double[numTotalVertices];
//     Int *activeElement = new Int[ne];
//     for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++)
//         numElementsAdj2Vertex[currentVertex] = 0;
//     for (currentElem = 0; currentElem < ne; currentElem++) {
//         activeElement[currentElem] = 1;
//         for (j = 0; j < nve; j++) {
//             currentVertex = mesh.t[j*ne+currentElem];
//             getNodalIndexOfVertex(&currentNode, j, nd, elemtype, porder);
//             numElementsAdj2Vertex[currentVertex] += 1;
//             if (numElementsAdj2Vertex[currentVertex] > 50)
//                 error("Over 50 elements are connected to a vertex.\n");
//             vertex2elem[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentElem;
//             vertex2node[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentNode;
//         }
//         for (j = 0; j < npv; j++) {
//             if (DGfield[currentElem*npv + j] <= 0.0) {
//                 activeElement[currentElem] = 0;
//                 break;
//             }
//         }
//     }
//     
//     double *avgFieldElem = new double[ne];
//     if (avgMethod == 0) {      // 0: Vertex value given by David's reconstruction
//         // Compute average field in each element:
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             avgFieldElem[currentElem] = 0.0;
//             for (j = 0; j < npv; j++) {
// //                 for (k = 0; k < nqvR; k++)
// //                     avgFieldElem[currentElem] += master.shapvgR[k*npv+j] * DGfield[currentElem*npv + j];
//                 avgFieldElem[currentElem] += DGfield[currentElem*npv + j];     // TODO: This is only an approximation to the true average. Some weighting is required for the exact average.
//             }
//             avgFieldElem[currentElem] /= ((double) npv);
//         }
// 
//         // Compute field value at each vertex (average of the average field in the adjacent elements):
//         for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
//             fieldVertices[currentVertex] = 0.0;
// //             neighMearures = 0.0;
//             for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
//                 currentElem = vertex2elem[currentVertex + j*numTotalVertices];
//                 fieldVertices[currentVertex] += avgFieldElem[currentElem];
// //                 fieldVertices[currentVertex] += elemMeasure[currentElem] * avgFieldElem[currentElem];
// //                 neighMearures += elemMeasure[currentElem];
//             }
//             fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
// //             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
//         }
//     }
//     else if (avgMethod == 1) {     // 1: Vertex value is average of nodes on vertex at different elements
//         // Compute field value at each vertex (average of the average field in the adjacent elements):
//         for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
//             fieldVertices[currentVertex] = 0.0;
// //             neighMearures = 0.0;
//             for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
//                 currentElem = vertex2elem[currentVertex + j*numTotalVertices];
//                 currentNode = vertex2node[currentVertex + j*numTotalVertices];
//                 fieldVertices[currentVertex] += DGfield[currentElem*npv + currentNode];
// //                 fieldVertices[currentVertex] += elemMeasure[currentElem] * DGfield[currentElem*npv + currentNode];
// //                 neighMearures += elemMeasure[currentElem];
//             }
//             fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
// //             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
//         }
//     }
//     else
//         error("Averaging method not recognized in DG_2_p1CG_2d.\n");
//     
//     // Compute p1CG field at each DG node:
//     double xi1, xi2, a, b, c;
//     double f1, f2, f3, f4, f12, f43;
//     if (elemtype == 0) {
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             f1 = fieldVertices[mesh.t[0*ne+currentElem]];
//             f2 = fieldVertices[mesh.t[1*ne+currentElem]];
//             f3 = fieldVertices[mesh.t[2*ne+currentElem]];
// 
//             a = -f1+f2;
//             b = -f1+f3;
//             c = f1;
// 
//             for (j = 0; j < npv; j++) {
//                 xi1 = master.plocvl[0*npv+j];
//                 xi2 = master.plocvl[1*npv+j];
// 
//                 if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
//                     p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c;
//                 else
//                     p1CGfield[currentElem*npv+j] = 0.0;
//             }
//         }
//     }
//     else if (elemtype == 1) {
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             f1 = fieldVertices[mesh.t[0*ne+currentElem]];
//             f2 = fieldVertices[mesh.t[1*ne+currentElem]];
//             f3 = fieldVertices[mesh.t[2*ne+currentElem]];
//             f4 = fieldVertices[mesh.t[3*ne+currentElem]];
// 
//             for (j = 0; j < npv; j++) {
//                 xi1 = master.plocvl[0*npv+j];
//                 xi2 = master.plocvl[1*npv+j];
// 
//                 f12 = (1.0-xi1)*f1 + xi1*f2;
//                 f43 = (1.0-xi1)*f4 + xi1*f3;
// 
//                 if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
//                     p1CGfield[currentElem*npv+j] = (1.0-xi2)*f12 + xi2*f43;
//                 else
//                     p1CGfield[currentElem*npv+j] = 0.0;
//             }
//         }
//     }
// 
//     delete[] numElementsAdj2Vertex; delete[] vertex2elem; delete[] vertex2node; delete[] avgFieldElem; delete[] fieldVertices; delete[] activeElement;
// }
// 
// void DG_2_p1CG_3d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
// {
//     // avgMethod:   0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
// 
//     Int i, j, currentElem, currentVertex, currentNode;
//     Int nd = ndims[0];
//     Int nve = ndims[3];
//     Int ne = ndims[5];
// //    Int nv = ndims[7];
//     Int npv = ndims[9];
//     Int nqvR = master.nqvR;
//     Int porder = app.porder;
//     Int elemtype = getElemtype(nd, nve);
//     double neighMearures;
// 
//     Int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element
// 
//     // Get number of total vertices in the mesh (vertices on periodic boundaries are counted twice, i.e., separately)
//     Int numTotalVertices = 0;
//     for (i = 0; i < ne*nve; i++) {
//         if (mesh.t[i] > numTotalVertices) {
//             numTotalVertices = mesh.t[i];
//         }
//     }
//     numTotalVertices += 1;
// //    if (nv != numTotalVertices) {
// //        error("There are ghost vertices in the mesh.\n");
// //    }
// 
//     // Creation of a matrix that contains the vertex-to-element connectivity:
//     Int *numElementsAdj2Vertex = new Int[numTotalVertices];
//     Int *vertex2elem = new Int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
//     Int *vertex2node = new Int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
//     double *fieldVertices = new double[numTotalVertices];
//     Int *activeElement = new Int[ne];
//     for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
//         numElementsAdj2Vertex[currentVertex] = 0;
//     }
//     for (currentElem = 0; currentElem < ne; currentElem++) {
//         activeElement[currentElem] = 1;
//         for (j = 0; j < nve; j++) {
//             currentVertex = mesh.t[j*ne+currentElem];
//             getNodalIndexOfVertex(&currentNode, j, nd, elemtype, porder);
//             numElementsAdj2Vertex[currentVertex] += 1;
//             if (numElementsAdj2Vertex[currentVertex] > 100) {
//                 error("Over 100 elements are connected to a vertex.\n");
//             }
//             vertex2elem[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentElem;
//             vertex2node[currentVertex + numTotalVertices*(numElementsAdj2Vertex[currentVertex]-1)] = currentNode;
//         }
//         for (j = 0; j < npv; j++) {
//             if (DGfield[currentElem*npv + j] <= 0.0) {
//                 activeElement[currentElem] = 0;
//                 break;
//             }
//         }
//     }
// 
//     double *avgFieldElem = new double[ne];
//     if (avgMethod == 0) {      // 0: Vertex value given by David's reconstruction
//         // Compute average field in each element:
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             avgFieldElem[currentElem] = 0.0;
//             for (j = 0; j < npv; j++) {
// //                 for (k = 0; k < nqvR; k++)
// //                     avgFieldElem[currentElem] += master.shapvgR[k*npv+j] * DGfield[currentElem*npv + j];
//                 avgFieldElem[currentElem] += DGfield[currentElem*npv + j];     // TODO: This is only an approximation to the true average. Some weighting is required for the exact average.
//             }
//             avgFieldElem[currentElem] = avgFieldElem[currentElem] / (double) npv;
//         }
// 
//         // Compute field value at each vertex (average of the average field in the adjacent elements):
//         for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
//             fieldVertices[currentVertex] = 0.0;
// //             neighMearures = 0.0;
//             for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
//                 currentElem = vertex2elem[currentVertex + j*numTotalVertices];
//                 fieldVertices[currentVertex] += avgFieldElem[currentElem];
// //                 fieldVertices[currentVertex] += elemMeasure[currentElem] * avgFieldElem[currentElem];
// //                 neighMearures += elemMeasure[currentElem];
//             }
//             fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
// //             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
//         }
//     }
//     else if (avgMethod == 1) {     // 1: Vertex value is average of nodes on vertex at different elements
//         // Compute field value at each vertex (average of the average field in the adjacent elements):
//         for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
//             fieldVertices[currentVertex] = 0.0;
// //             neighMearures = 0.0;
//             for (j = 0; j < numElementsAdj2Vertex[currentVertex]; j++) {
//                 currentElem = vertex2elem[currentVertex + j*numTotalVertices];
//                 currentNode = vertex2node[currentVertex + j*numTotalVertices];
//                 fieldVertices[currentVertex] += DGfield[currentElem*npv + currentNode];
// //                 fieldVertices[currentVertex] += elemMeasure[currentElem] * DGfield[currentElem*npv + currentNode];
// //                 neighMearures += elemMeasure[currentElem];
//             }
//             fieldVertices[currentVertex] /= max((double) numElementsAdj2Vertex[currentVertex],1.0);
// //             fieldVertices[currentVertex] /= max(neighMearures,1.0e-8);
//         }
//     }
//     else
//         error("Averaging method not recognized in DG_2_p1CG_3d.\n");
// 
//     // Compute p1CG field   at each DG node:
//     double xi1, xi2, xi3, a, b, c, d;
//     double f1, f2, f3, f4, f5, f6, f7, f8, f12, f43, f56, f87, f1234, f5678;
//     if (elemtype == 0) {
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             f1 = fieldVertices[mesh.t[0*ne+currentElem]];
//             f2 = fieldVertices[mesh.t[1*ne+currentElem]];
//             f3 = fieldVertices[mesh.t[2*ne+currentElem]];
//             f4 = fieldVertices[mesh.t[3*ne+currentElem]];
// 
//             a = -f1+f2;
//             b = -f1+f3;
//             c = -f1+f4;
//             d = f1;
// 
//             for (j = 0; j < npv; j++) {
//                 xi1 = master.plocvl[0*npv+j];
//                 xi2 = master.plocvl[1*npv+j];
//                 xi3 = master.plocvl[2*npv+j];
// 
//                 if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
//                     p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c*xi3 + d;
//                 else
//                     p1CGfield[currentElem*npv+j] = 0.0;
//             }
//         }
//     }
//     else if (elemtype == 1) {
//         for (currentElem = 0; currentElem < ne; currentElem++) {
//             f1 = fieldVertices[mesh.t[0*ne+currentElem]];
//             f2 = fieldVertices[mesh.t[1*ne+currentElem]];
//             f3 = fieldVertices[mesh.t[2*ne+currentElem]];
//             f4 = fieldVertices[mesh.t[3*ne+currentElem]];
//             f5 = fieldVertices[mesh.t[4*ne+currentElem]];
//             f6 = fieldVertices[mesh.t[5*ne+currentElem]];
//             f7 = fieldVertices[mesh.t[6*ne+currentElem]];
//             f8 = fieldVertices[mesh.t[7*ne+currentElem]];
// 
//             for (j = 0; j < npv; j++) {
//                 xi1 = master.plocvl[0*npv+j];
//                 xi2 = master.plocvl[1*npv+j];
//                 xi3 = master.plocvl[2*npv+j];
// 
//                 f12 = (1.0-xi1)*f1 + xi1*f2;
//                 f43 = (1.0-xi1)*f4 + xi1*f3;
//                 f1234 = (1.0-xi2)*f12 + xi2*f43;
// 
//                 f56 = (1.0-xi1)*f5 + xi1*f6;
//                 f87 = (1.0-xi1)*f8 + xi1*f7;
//                 f5678 = (1.0-xi2)*f56 + xi2*f87;
// 
//                 if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
//                     p1CGfield[currentElem*npv+j] = (1.0-xi3)*f1234 + xi3*f5678;
//                 else
//                     p1CGfield[currentElem*npv+j] = 0.0;
//             }
//         }
//     }
//     
//     delete[] numElementsAdj2Vertex; delete[] vertex2elem; delete[] vertex2node; delete[] avgFieldElem; delete[] fieldVertices; delete[] activeElement;
// }
// 
// void DG_2_p1CG(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims, Int avgMethod)
// {
//     Int nd = ndims[0];
// 
//     if (nd == 2)
//         DG_2_p1CG_2d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
//     else if (nd == 3)
//         DG_2_p1CG_3d(p1CGfield, DGfield, mesh, master, app, ndims, avgMethod);
//     else
//         error("Invalid number of dimensions.\n");
// }
// 
// #endif
