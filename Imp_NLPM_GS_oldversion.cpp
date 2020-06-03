// mex file, using gauss-seidel algorithm for NLPM implicit scheme
//  input: g, dt, f
// output: u 
/*
(1-dt*C)*U^(n+1) = U^n
A*U^(n+1) = U^n
Aij   = 1+dt[(gij+gi-1j)/2 + (gij+gi+1j)/2 + (gij+gij-1)/2 + (gij+gij+1)/2]
Ai-1j = -dt(gij+gi-1j)/2
Ai+1j = -dt(gij+gi+1j)/2
Aij-1 = -dt(gij+gij-1)/2
Aij+1 = -dt(gij+gij+1)/2
*/

/*
 *Input: 
             *g: 
             *V_rt: ...................... The right part 
             dt:       
 *Output:   
             *V: ........................ V^(n+1)
             *A   //no
             *Lap__ //no
             *R__   //no
*/
/* Revision: 1.0 */
/*Copyright(C) Ying Wen 05/22/2020*/
/*Can not deal with large image, because A is N*N*/

// mex -I'/usr/local/opt/openblas/include' -L'/usr/local/opt/openblas/lib' -lblas Imp_NLPM_GS.cpp

#include "mex.h"
#include "matrix.h"
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
// #include "Nystrom.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <cblas.h>
// #include <omp.h>
using namespace std;



void mexFunction(int nlhs, mxArray *plhs[], //ouput
                int nrhs, const mxArray *prhs[]) //input
{
    if (nrhs != 3)
        mexErrMsgTxt("Five inputs required.");
    if (nlhs != 2)
        mexErrMsgTxt("One output required.");

    double *g, *V_rt;
    double dt;
    mwSize NX, NY;
    g = mxGetPr(prhs[0]);
    V_rt = mxGetPr(prhs[1]);//R
    dt = mxGetScalar(prhs[2]);
    NX = mxGetM(prhs[1]);
    NY = mxGetN(prhs[1]);

    int Nx, Ny;
    Nx = NX;
    Ny = NY;
    double *V, *A;
    plhs[0] = mxCreateDoubleMatrix(Nx, Ny, mxREAL);
    V = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(Nx*Ny, Nx*Ny, mxREAL);
    A = mxGetPr(plhs[1]);
    int N = Nx*Ny;
    mexPrintf("N:%d, N*N:%d\n", N, N*N);
    // double *A = new double[N*N]; // Matrix A
    
    mexPrintf("*************\n");
    clock_t startTime, endTime;

    // Prepare Matrix A, A is symmetirc
    /*
    Aij   = 1+dt[(gij+gi-1j)/2 + (gij+gi+1j)/2 + (gij+gij-1)/2 + (gij+gij+1)/2]
    Ai-1j = -dt(gij+gi-1j)/2
    Ai+1j = -dt(gij+gi+1j)/2
    Aij-1 = -dt(gij+gij-1)/2
    Aij+1 = -dt(gij+gij+1)/2
    */
    startTime = clock();
    int Nrow, Ncl1, Ncl2, Ncl3, Ncl4, Ncl5;
    int i_,j_,i__, j__;
    int A_c;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            mexPrintf("i:%d, j:%d,  N:%d*************\n", i,j,N);
            // j*Nx+i row number of A
            Nrow = j*Nx+i;
            Ncl1 = j*Nx+i;
            i_ = (i==0)?i:i-1;
            i__ = (i==Nx-1)?i:i+1;
            j_ = (j==0)?j:j-1;
            j__ = (j==Ny-1)?j:j+1;

            Ncl2 = j*Nx+i__;
            Ncl3 = j*Nx+i_;
            Ncl4 = j__*Nx+i;
            Ncl5 = j_*Nx+i;
            // mexPrintf("%d,%d,%d,%d,%d,%d,%d,%d,%d\n", i_,i__,j_,j__, Ncl1, Ncl2, Ncl3, Ncl4, Ncl5);

            A_c = Nrow*N;
            mexPrintf("%d,%d,%d,%d,   %d,%d,%d,%d,%d  %d\n", i_,i__,j_,j__, Ncl1, Ncl2, Ncl3, Ncl4, Ncl5, A_c);
            double gij = g[Ncl1];
            double g1  = (gij + g[Ncl3])*dt/2;
            double g2  = (gij + g[Ncl2])*dt/2;
            double g3  = (gij + g[Ncl5])*dt/2;
            double g4  = (gij + g[Ncl4])*dt/2;


            A[A_c+Ncl1] += 1 + g1 + g2 + g3 + g4;
            A[A_c+Ncl2] += -1*g2;
            A[A_c+Ncl3] += -1*g1;
            A[A_c+Ncl4] += -1*g4;
            A[A_c+Ncl5] += -1*g3;
        }    
        // mexPrintf("i=%d \n",i);    
    }
    endTime = clock();
	mexPrintf("Prepare A time : %f s \n", (double)(endTime - startTime) / CLOCKS_PER_SEC );
    
    // solve V from AV=V_rt, using G-S
    startTime = clock();
    double* V_old = new double[N];
    double* V_dif = new double[N];
    double eps = 0.001;
    double infnm_, stp, V_infnm;//infnm__, 
    stp = eps;
    //initial V = V_rt/A(diag)
    for (int i = 0; i < N; i++)
    {
        V[i] = V_rt[i]/A[i*N+i];
    }
    double Rem; //Remainder
    for (int i = 0; i<200 & stp>=eps; i++)
    {
        cblas_dcopy(N, V, 1, V_old, 1);
        // mexPrintf("-1");
        infnm_ = 0;
        for (int i = 0; i < N; i++)
        {
            Rem = 0;//A[i*N+i]*V[i];
            if (i-1>=0)
            {
                Rem += A[i*N+i-1]*V[i-1];
            }
            if (i-Nx>=0)
            {
                Rem += A[i*N+i-Nx]*V[i-Nx];
            }
            if (i+1<N)
            {
                Rem += A[i*N+i+1]*V[i+1];
            }
            if (i+Nx<N)
            {
                Rem += A[i*N+i+Nx]*V[i+Nx];
            }

            V[i] = (V_rt[i] - Rem)/A[i*N+i];
            
            V_dif[i] = abs(V[i] - V_old[i]);
            infnm_ = (V_dif[i]>infnm_)?V_dif[i]:infnm_;
        }
        // compare the infinite norm   
        stp = infnm_;    ///infnm__
        mexPrintf("i: %d, stop infinite norm: %f, infnm_: %f \n", i, stp, infnm_); 
        // mexPrintf("-----------------------------------------\n");
    }
    endTime = clock();
    mexPrintf("G-S, V, AV=V_rt, time: %f s\n",  (double)(endTime - startTime) / CLOCKS_PER_SEC );

    free(V_old);
    free(V_dif);
    // free(A);
}



