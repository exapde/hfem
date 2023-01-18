#ifndef __DG2P1CG
#define __DG2P1CG

// Written by C. Nguyen and P. Fernandez

void DG_2_p1CG_2d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims)
{
    int i, j, currentElem, currentVertex, currentNode, elemtype;
    int nd = (int) ndims[0];
    int nve = (int) ndims[3];
    int ne = (int) ndims[5];
//    int nv = (int) ndims[7];
    int npv = (int) ndims[9];
    int porder = (int) app.porder;

    int avgMethod = 1;        // 0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
    int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element

    if (nve == 3) {     // Tri
        elemtype = 0;
    }
    else if (nve == 4) {        // Quad
        elemtype = 1;
    }
    else {
        error("Invalid number of vertices per element. Neither tri nor quad.\n");
    }

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
    int *numElementsAdj2Vertex = new int[numTotalVertices];
    int *vertex2elem = new int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    int *vertex2node = new int[numTotalVertices*50];      // We assume at most 50 elements are connected to a vertex.
    double *fieldVertices = new double[numTotalVertices];
    int *activeElement = new int[ne];
    for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
        numElementsAdj2Vertex[currentVertex] = 0;
    }
    for (currentElem = 0; currentElem < ne; currentElem++) {
        activeElement[currentElem] = 1;
        for (j = 0; j < nve; j++) {
            currentVertex = mesh.t[j*ne+currentElem];
            getVertexIndex(&currentNode, j, nd, elemtype, porder);
            numElementsAdj2Vertex[currentVertex] += 1;
            if (numElementsAdj2Vertex[currentVertex] > 50) {
                error("Over 50 elements are connected to a vertex.\n");
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

    // Compute p1CG field at each DG node:
    double xi1, xi2, a, b, c;
    double f1, f2, f3, f4, f12, f43;
    if (elemtype == 0) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[currentElem]];
            f2 = fieldVertices[mesh.t[ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];

            a = -f1+f2;
            b = -f1+f3;
            c = f1;

            for (j = 0; j < npv; j++) {
                xi1 = master.plocnvl[j];
                xi2 = master.plocnvl[npv+j];

                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
    else if (elemtype == 1) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[currentElem]];
            f2 = fieldVertices[mesh.t[ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];

            for (j = 0; j < npv; j++) {
                xi1 = master.plocnvl[j];
                xi2 = master.plocnvl[npv+j];

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

void DG_2_p1CG_3d(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims)
{
    int i, j, currentElem, currentVertex, currentNode, elemtype;
    int nd = (int) ndims[0];
    int nve = (int) ndims[3];
    int ne = (int) ndims[5];
//    int nv = (int) ndims[7];
    int npv = (int) ndims[9];
    int porder = (int) app.porder;

    int avgMethod = 1;        // 0: Vertex value given by David's reconstruction. 1: Vertex value is average of nodes on vertex at different elements
    int zeroPassiveElemets = 0;     // 1: Zero out elements for which the field is nonpositive on a nonempty portion of the element

    if (nve == 4) {     // Tet
        elemtype = 0;
    }
    else if (nve == 8) {        // Hexa
        elemtype = 1;
    }
    else {
        error("Invalid number of vertices per element. Neither tet nor hexa.\n");
    }

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
    int *numElementsAdj2Vertex = new int[numTotalVertices];
    int *vertex2elem = new int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
    int *vertex2node = new int[numTotalVertices*100];      // We assume at most 50 elements are connected to a vertex.
    double *fieldVertices = new double[numTotalVertices];
    int *activeElement = new int[ne];
    for (currentVertex = 0; currentVertex < numTotalVertices; currentVertex++) {
        numElementsAdj2Vertex[currentVertex] = 0;
    }
    for (currentElem = 0; currentElem < ne; currentElem++) {
        activeElement[currentElem] = 1;
        for (j = 0; j < nve; j++) {
            currentVertex = mesh.t[j*ne+currentElem];
            getVertexIndex(&currentNode, j, nd, elemtype, porder);
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

    // Compute p1CG field at each DG node:
    double xi1, xi2, xi3, a, b, c, d;
    double f1, f2, f3, f4, f5, f6, f7, f8, f12, f43, f56, f87, f1234, f5678;
    if (elemtype == 0) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[currentElem]];
            f2 = fieldVertices[mesh.t[ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];

            a = -f1+f2;
            b = -f1+f3;
            c = -f1+f4;
            d = f1;

            for (j = 0; j < npv; j++) {
                xi1 = master.plocnvl[j];
                xi2 = master.plocnvl[npv+j];
                xi3 = master.plocnvl[2*npv+j];

                if (activeElement[currentElem] == 1 || zeroPassiveElemets == 0)
                    p1CGfield[currentElem*npv+j] = a*xi1 + b*xi2 + c*xi3 + d;
                else
                    p1CGfield[currentElem*npv+j] = 0.0;
            }
        }
    }
    else if (elemtype == 1) {
        for (currentElem = 0; currentElem < ne; currentElem++) {
            f1 = fieldVertices[mesh.t[currentElem]];
            f2 = fieldVertices[mesh.t[ne+currentElem]];
            f3 = fieldVertices[mesh.t[2*ne+currentElem]];
            f4 = fieldVertices[mesh.t[3*ne+currentElem]];
            f5 = fieldVertices[mesh.t[4*ne+currentElem]];
            f6 = fieldVertices[mesh.t[5*ne+currentElem]];
            f7 = fieldVertices[mesh.t[6*ne+currentElem]];
            f8 = fieldVertices[mesh.t[7*ne+currentElem]];

            for (j = 0; j < npv; j++) {
                xi1 = master.plocnvl[j];
                xi2 = master.plocnvl[npv+j];
                xi3 = master.plocnvl[2*npv+j];

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

void DG_2_p1CG(double *p1CGfield, double *DGfield, meshstruct &mesh, masterstruct &master, appstruct &app, Int* ndims)
{
    int nd = (int) ndims[0];

    if (nd == 2) {
        DG_2_p1CG_2d(p1CGfield, DGfield, mesh, master, app, ndims);
    }
    else if (nd == 3) {
        DG_2_p1CG_3d(p1CGfield, DGfield, mesh, master, app, ndims);
    }
    else {
        error("Invalid number of dimensions.\n");
    }
    
}

#endif
