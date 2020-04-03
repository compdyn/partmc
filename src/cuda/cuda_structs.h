//
// Created by Christian on 01/04/2020.
//

#ifndef CAMPGPU_CUDA_STRUCTS_H
#define CAMPGPU_CUDA_STRUCTS_H

typedef struct
{
    //Init variables ("public")
    int threads,blocks;
    int maxIt;
    int mattype;
    int nrows;
    int nnz;
    int n_cells;
    double tolmax;
    double* ddiag;

    // Intermediate variables ("private")
    double * dr0;
    double * dr0h;
    double * dn0;
    double * dp0;
    double * dt;
    double * ds;
    double * dAx2;
    double * dy;
    double * dz;

    // Matrix data ("input")
    double* A;
    int*    jA;
    int*    iA;

    //GPU pointers ("input")
    double* dA;
    int*    djA;
    int*    diA;
    double* dx;
    double* aux;
    double* daux;

    // ODE solver variables
    double* dewt;
    double* dacor;
    double* dacor_init;
    double* dtempv;
    double* dftemp;
    double* dzn;
    double* dcv_y;

} itsolver;


#endif //CAMPGPU_CUDA_STRUCTS_H
