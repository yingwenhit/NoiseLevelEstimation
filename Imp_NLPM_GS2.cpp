// mex file, using gauss-seidel algorithm for NLPM implicit scheme
//  input: g, dt, f
// output: u


/*
 *Input: 
             *g: 
             *V_rt: ...................... The right part 
             dt:       
 *Output:   
             *V: ........................ V^(n+1)
*/
/* Revision: 2.0 */
/*Copyright(C) Ying Wen 05/26/2020*/
/*Change A to a [5*N] matrix, to reduce the storage consuming*/

// mex -I'/usr/local/opt/openblas/include' -L'/usr/local/opt/openblas/lib' -lblas Imp_NLPM_GS2.cpp

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

void mexFunction(int nlhs, mxArray *plhs[],       //ouput
                 int nrhs, const mxArray *prhs[]) //input
{
    if (nrhs != 3)
        mexErrMsgTxt("Five inputs required.");
    if (nlhs != 1)
        mexErrMsgTxt("One output required.");

    double *g, *V_rt;
    double dt;
    mwSize NX, NY;
    g = mxGetPr(prhs[0]);
    V_rt = mxGetPr(prhs[1]); //R
    dt = mxGetScalar(prhs[2]);
    NX = mxGetM(prhs[1]);
    NY = mxGetN(prhs[1]);

    int Nx, Ny;
    Nx = NX;
    Ny = NY;
    int N = Nx * Ny;
    double *V, *A;
    plhs[0] = mxCreateDoubleMatrix(Nx, Ny, mxREAL);
    V = mxGetPr(plhs[0]);

    // plhs[1] = mxCreateDoubleMatrix(N, 5, mxREAL);
    // A = mxGetPr(plhs[1]);

    A = new double[5 * N]; // Matrix A

    // mexPrintf("*************\n");
    clock_t startTime, endTime;

    // Prepare Matrix A, A is symmetirc
    // startTime = clock();
    int Nrow, Ncl[5]; //Ncl1, Ncl2, Ncl3, Ncl4, Ncl5;
    int i_, j_, i__, j__;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            // j*Nx+i row number of A
            Nrow = j * Nx + i;
            Ncl[0] = j * Nx + i;
            i_ = (i == 0) ? i : i - 1;
            i__ = (i == Nx - 1) ? i : i + 1;
            j_ = (j == 0) ? j : j - 1;
            j__ = (j == Ny - 1) ? j : j + 1;

            Ncl[1] = j * Nx + i__; // +N
            Ncl[2] = j * Nx + i_; // +3N
            Ncl[3] = j__ * Nx + i; //0
            Ncl[4] = j_ * Nx + i; // +4N

            double g_[5];
            g_[0] = g[Ncl[0]] * dt;
            g_[1] = (g_[0] + g[Ncl[1]] * dt) / 2;
            g_[2] = (g_[0] + g[Ncl[2]] * dt) / 2;
            g_[3] = (g_[0] + g[Ncl[3]] * dt) / 2;
            g_[4] = (g_[0] + g[Ncl[4]] * dt) / 2;

            A[2 * N + Nrow] += 1 + g_[1] + g_[2] + g_[3] + g_[4];
            for (int k = 1; k < 5; k++)
            {
                int N_s=0;
                if (Ncl[k] - Nrow == 0)
                {
                    N_s = 2 * N;
                }
                if (Ncl[k] - Nrow == 1)
                    N_s = N;

                if (Ncl[k] - Nrow == -1)
                    N_s = 3*N;

                if (Ncl[k] - Nrow == Nx)
                    N_s = 0;

                if (Ncl[k] - Nrow == -1 * Nx)
                    N_s = 4*N;


                // if (N_s + Nrow>=5*N)
                //     mexPrintf("ERROR!\n");

                A[N_s + Nrow] += -1 * g_[k];
            }
        }
        // mexPrintf("i=%d \n",i);
    }
    // endTime = clock();
    // mexPrintf("Prepare A time : %f s \n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

    // solve V from AV=V_rt, using G-S
    // startTime = clock();
    double *V_old = new double[N];
    double *V_dif = new double[N];
    double eps = 0.001;
    double infnm_, stp, V_infnm; //infnm__,
    stp = eps;
    //initial V = V_rt/A(diag)
    for (int i = 0; i < N; i++)
    {
        V[i] = V_rt[i] / A[2 * N + i];
    }

    double Rem; //Remainder
    for (int i = 0; i < 1000 & stp >= eps; i++)
    {
        cblas_dcopy(N, V, 1, V_old, 1);
        infnm_ = 0;
        for (int i = 0; i < N; i++)
        {
            Rem = 0; 
            if (i-1>=0)
            {
                Rem += A[3*N+i]*V[i-1];// i-1: 3N
                
            }
            if (i-Nx>=0)
            {
                Rem += A[4*N+i]*V[i-Nx];// I-Nx 4N
            }
            if (i+1<N)
            {
                Rem += A[N+i]*V[i+1];//i+1 N
            }
            if (i+Nx<N)
            {
                Rem += A[i]*V[i+Nx];//i+Nx 0
            }

            V[i] = (V_rt[i] - Rem)/A[2*N+i];

            V_dif[i] = abs(V[i] - V_old[i]);
            infnm_ = (V_dif[i] > infnm_) ? V_dif[i] : infnm_;
        }
        // compare the infinite norm
        stp = infnm_; 
        // mexPrintf("i: %d, stop infinite norm: %f, infnm_: %f \n", i, stp, infnm_);
    }
    // endTime = clock();
    // mexPrintf("G-S, V, AV=V_rt, time: %f s\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

    free(V_old);
    free(V_dif);
    free(A);
}
