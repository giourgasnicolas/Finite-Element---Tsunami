# include "tsunami.h"

double* E;
double* U;
double* V;
double* FE;
double* FU;
double* FV;
double* bath;
int* MeshElem;
double* X;
double* Y;
int nElem, nNode, nEdge;
int** node;
int** elem; 

#define interpolate_edge(phi,U, map) ((phi)[0]*(U)[(map)[0]] + (phi)[1]*(U)[(map)[1]])
#define interpolate_triangle(phi,U,map) ((phi)[0]*(U)[(map)[0]] + (phi)[1]*(U)[(map)[1]] + (phi)[2]*(U)[(map)[2]])

void femShallowEdgeMap(int index, int map[2][2])
{
    int i, j, k, FuncNode, FuncElem;
    for (j=0; j < 2; ++j) {
        FuncNode = node[index][j];
        for (k=0; k < 2; k++) {
            FuncElem = elem[index][k];
            map[k][j] = (nElem)*3;
            if (FuncElem >= 0) {
                for (i=0; i < 3; i++) {
                    if (MeshElem[FuncElem*3 + i] == FuncNode) {
                        map[k][j] = FuncElem*3 + i;  
                    }
                }
            }
        }
    }
}
void femShallowMultiplyInverseMatrix()
{         
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}}; 
    double BEloc[3], BUloc[3], BVloc[3], xLoc[3], yLoc[3], jac; 
    int    FuncElem,i,j,mapElem[3];

    for (FuncElem=0; FuncElem < nElem; FuncElem++) {
        
        for (j = 0; j < 3; ++j)
            mapElem[j] = FuncElem * 3 + j;

        int *mapCoord = &(MeshElem[FuncElem*3]);
        
        for (i=0; i < 3; ++i) {
        	  xLoc[i] = X[mapCoord[i]];
        	  yLoc[i] = Y[mapCoord[i]]; 
              BEloc[i] = FE[mapElem[i]];
              BUloc[i] = FU[mapElem[i]];
              BVloc[i] = FV[mapElem[i]];
              FE[mapElem[i]] = 0.0;
              FU[mapElem[i]] = 0.0;
              FV[mapElem[i]] = 0.0;
        }
        
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                FE[mapElem[i]] += invA[i][j] * BEloc[j] / jac; 
                FU[mapElem[i]] += invA[i][j] * BUloc[j] / jac; 
                FV[mapElem[i]] += invA[i][j] * BVloc[j] / jac; 
            }
        }
    }
}
void femShallowAddIntegralsElements()
{    
    double  xLoc[3], yLoc[3], dphidx[3], dphidy[3], phi[3], FuncBath[3], xsi, eta, weight, jac, eta_h, u, v, y, x, h, z3d, f, fact;
    int     i,j,k,FuncElem,mapElem[3];
    for (FuncElem=0; FuncElem < nElem; FuncElem++){
        
        for (j = 0; j < 3; ++j)
            mapElem[j] = FuncElem * 3 + j;
        
        int *mapCoord = &(MeshElem[FuncElem*3]);
        for (j=0; j < 3; ++j){
            xLoc[j] = X[mapCoord[j]];
            yLoc[j] = Y[mapCoord[j]];
            FuncBath[j] = bath[mapCoord[j]];
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac; 
        
        for (k=0; k < 3; k++){
            xsi = gaussTriangleXsi[k];
            eta = gaussTriangleEta[k];
            weight = gaussTriangleWeight[k];     
            
            phi[0] = 1 - xsi - eta;  
            phi[1] = xsi;
            phi[2] = eta;
            
            u = interpolate_triangle(phi,U,mapElem); 
            v = interpolate_triangle(phi,V,mapElem); 
            eta_h = interpolate_triangle(phi,E,mapElem);
            
            y = phi[0]*yLoc[0] + phi[1]*yLoc[1] + phi[2]*yLoc[2]; 
            x = phi[0]*xLoc[0] + phi[1]*xLoc[1] + phi[2]*xLoc[2];
            h = phi[0]*FuncBath[0] + phi[1]*FuncBath[1] + phi[2]*FuncBath[2];
 
            z3d   = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
            f     = 2*Omega*(z3d/R);
            fact  = (4*R*R + x*x + y*y)/(4*R*R);

            for (i=0; i < 3; i++){
                FE[mapElem[i]] += (((dphidx[i]*h*u + dphidy[i]*h*v)*(fact)) + phi[i]*((h*(x*u+y*v))/(R*R)))*jac*weight;
                FU[mapElem[i]] += (phi[i]*(f*v - Gamma*u) + dphidx[i]*g*eta_h*(fact) + phi[i]*((g*x*eta_h)/(2*R*R)))*jac*weight;
                FV[mapElem[i]] += (phi[i]*(-f*u - Gamma*v) + dphidy[i]*g*eta_h*(fact)+ phi[i]*((g*y*eta_h)/(2*R*R)))*jac*weight;
            }
        }
    }
}
void femShallowAddIntegralsEdges()
{
    double  xEdge[2], yEdge[2], phiEdge[2], FuncBath[2], xsi, weight, jac, eta_rix, u, v, eta_L, eta_R, u_L, u_R, dxdxsi, dydxsi, norm, y, x, h, fact, fact2;
    int     i,j,k,FuncEdge,mapEdge[2][2];
    for (FuncEdge=0; FuncEdge < nEdge; FuncEdge++){
        femShallowEdgeMap(FuncEdge,mapEdge);
        for (j=0; j < 2; ++j){
            int FuncNode = node[FuncEdge][j];
            xEdge[j] = X[FuncNode];
            yEdge[j] = Y[FuncNode];
            FuncBath[j] = bath[FuncNode];
        }

        dxdxsi = xEdge[1] - xEdge[0]; 
        dydxsi = yEdge[1] - yEdge[0]; 
        norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double normal[2] = {dydxsi/norm, -dxdxsi/norm};   //ICI
        jac = norm / 2.0;
        
        for (k=0; k < 2; k++){
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k];     
            phiEdge[0] = (1.0 - xsi)/2.0; phiEdge[1] = (1.0 + xsi)/2.0; 
            
            y = phiEdge[0]*yEdge[0] + phiEdge[1]*yEdge[1]; 
            x = phiEdge[0]*xEdge[0] + phiEdge[1]*xEdge[1];
            h = phiEdge[0]*FuncBath[0] + phiEdge[1]*FuncBath[1];

            if (elem[FuncEdge][1] == -1){ 
                eta_L = interpolate_edge(phiEdge, E, mapEdge[0]);
                eta_R = eta_L;
                u_L   = interpolate_edge(phiEdge, U, mapEdge[0])*normal[0] + interpolate_edge(phiEdge, V, mapEdge[0])*normal[1];
                u_R   = - u_L;
            }

            else{
                eta_L = interpolate_edge(phiEdge, E, mapEdge[0]);
                eta_R = interpolate_edge(phiEdge, E, mapEdge[1]);
                u_L   = interpolate_edge(phiEdge, U, mapEdge[0])*normal[0] + interpolate_edge(phiEdge, V, mapEdge[0])*normal[1];
                u_R   = interpolate_edge(phiEdge, U, mapEdge[1])*normal[0] + interpolate_edge(phiEdge, V, mapEdge[1])*normal[1];
            }

            eta_rix = (eta_L + eta_R)/2.0 + sqrt(h/g)*(u_L - u_R)/2.0;   //solveur de Rieman
            u       = (u_L + u_R)/2.0     + sqrt(g/h)*(eta_L - eta_R)/2.0;
            
            fact  = ((4*R*R + x*x + y*y)/(4*R*R))*jac*weight;
            fact2 = g*eta_rix*(fact);
            
            for (i=0; i < 2; i++){
                FE[mapEdge[0][i]] -= phiEdge[i]*h*u*(fact);
                FE[mapEdge[1][i]] += phiEdge[i]*h*u*(fact);
                FU[mapEdge[0][i]] -= phiEdge[i]*normal[0]*fact2;
                FU[mapEdge[1][i]] += phiEdge[i]*normal[0]*fact2;
                FV[mapEdge[0][i]] -= phiEdge[i]*normal[1]*fact2;
                FV[mapEdge[1][i]] += phiEdge[i]*normal[1]*fact2;
            }
        }
    }
}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    int i, j, n, trash;
    double dtrash;
    
    FILE* file = fopen(meshFileName,"r");

    fscanf(file, "Number of nodes %d \n", &nNode); 
    bath = malloc(sizeof(double) * nNode);
    X    = malloc(sizeof(double) * nNode); 
    Y    = malloc(sizeof(double) * nNode); 
    for (i = 0; i < nNode; i++)
        fscanf(file, "%d : %le %le %le \n", &trash, &X[i], &Y[i], &bath[i]);

    fscanf(file, "Number of triangles %d \n", &nElem);
    MeshElem = malloc(sizeof(int) * 3 * nElem);
    for (i = 0; i < nElem; i++)
        fscanf(file, "%d : %d %d %d \n", &trash, &MeshElem[i * 3], &MeshElem[i * 3 + 1], &MeshElem[i * 3 + 2]);

    fscanf(file, "Number of edges %d \n", &nEdge);
    elem = malloc(sizeof(int*) * nEdge);
    node = malloc(sizeof(int*) * nEdge);
    for (i = 0; i < nEdge; i++) {
        elem[i] = malloc(sizeof(int) * 2);
        node[i] = malloc(sizeof(int) * 2);
        fscanf(file, "%d : %d %d : %d %d \n", &trash, &node[i][0], &node[i][1], &elem[i][0], &elem[i][1]);
    }

    E  = calloc(nElem * 3, sizeof(double));
    V  = calloc(nElem * 3, sizeof(double));
    U  = calloc(nElem * 3, sizeof(double));

    double* PE  = calloc(nElem * 3, sizeof(double));
    double* K1E  = calloc(1+nElem * 3, sizeof(double));
    double* K2E  = calloc(1+nElem * 3, sizeof(double));

    double* PU  = calloc(nElem * 3, sizeof(double));
    double* K1U  = calloc(1 + nElem * 3, sizeof(double));
    double* K2U  = calloc(1 + nElem * 3, sizeof(double));

    double* PV  = calloc(nElem * 3, sizeof(double));
    double* K1V  = calloc(1 + nElem * 3, sizeof(double));
    double* K2V  = calloc(1 + nElem * 3, sizeof(double));

    double* Eold = E;
    double* Uold = U;
    double* Vold = V;
    
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i * 3 + j] = tsunamiInitialConditionOkada(X[MeshElem[i * 3 + j]], Y[MeshElem[i * 3 + j]]);

    for(n = 1; n <= nmax; n++){
        FE = K1E;
        FU = K1U;
        FV = K1V;

        femShallowAddIntegralsElements();
        femShallowAddIntegralsEdges();
        femShallowMultiplyInverseMatrix(); 
        for(i=0; i < 3*nElem; i++){
            PE[i] = E[i] + dt*K1E[i];
            PU[i] = U[i] + dt*K1U[i];
            PV[i] = V[i] + dt*K1V[i];
        }
        E = PE;
        U = PU;
        V = PV;
        FE = K2E;
        FU = K2U;
        FV = K2V;
        femShallowAddIntegralsElements();
        femShallowAddIntegralsEdges();
        femShallowMultiplyInverseMatrix(); 
        E = Eold;
        U = Uold;
        V = Vold;
        for (i=0; i < 3*nElem; i++) {
            E[i] += (dt/2.0)*(K1E[i]+K2E[i]);
            U[i] += (dt/2.0)*(K1U[i]+K2U[i]);
            V[i] += (dt/2.0)*(K1V[i]+K2V[i]);
            K1E[i] = 0; K1U[i] = 0;K1V[i] = 0;
            K2E[i] = 0; K2U[i] = 0;K2V[i] = 0;
        }
        FU = K1U;
        FV = K1V;
        FE = K1E;
        if ((n)%sub == 0)
            tsunamiWriteFile(baseResultName,n,U,V,E,nElem,3); 
    }
    free(E);
    free(U);
    free(V);
    free(PV);
    free(PE);
    free(PU);
    free(K1U);
    free(K1E);
    free(K1V);
    free(K2E);
    free(K2U);
    free(K2V);
    free(bath);
    free(MeshElem);
    free(X);
    free(Y);
    for (i = 0; i < nEdge; i++) {
        free(elem[i]);
        free(node[i]);
    }
    free(elem);
    free(node); 
}