// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

#ifndef min
#define min(a,b) ((a < b) ? a : b)
#endif
#ifndef max
#define max(a,b) ((a > b) ? a : b)
#endif

typedef struct _matrixSize      // Optional Command-line multiplier for matrix sizes
{
    unsigned int uiWA, uiHA, uiWB, uiHB, uiWC, uiHC;
} sMatrixSize;
void
matrixMulCPU(float* C, const float* A, const float* B, unsigned int hA, unsigned int wA, unsigned int wB);

void randomInit(float* data, int size);

void printDiff(float* data1, float* data2, int width, int height, int iListLength, float fListTol);

void initializeCUDA_random(int argc, char** argv, int& devID, int& iSizeMultiple, sMatrixSize& matrix_size);
void initializeCUDA(int argc, char** argv, int& devID, int& numberPoints, sMatrixSize& matrix_size);
int matrixMultiply(int argc, char** argv, float* h_A, float* h_B, float* h_C, int devID, sMatrixSize& matrix_size);
int matrixMultiply_random(int argc, char** argv, int devID, sMatrixSize& matrix_size);
