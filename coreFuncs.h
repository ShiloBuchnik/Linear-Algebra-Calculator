#include <stdbool.h>

int negativeZerosAndFindRank(int numRow, int numColumn, double* matrix);
double deterCalc(int n, double* originMatrix);
void findAndPrintAdjoint(int n, double* matrix);
void gaussJordanAndFindInverse(int numRow, int numColumn, double* inputMatrix, double* inverseMatrix,
                               double* permutation, bool* isPermutationIdentity, double* L, double* U);
double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, double* matrix1, double* matrix2);
void LUDecomposition(int numRow, int numColumn, double* inputMatrix, double* Permutation, double* L, double* U);
void findDVFromUAndPrint(int numRow, double* U);
double* leastSquares(int numRow, int numColumn, double* A, double** x, double* b);
void QRDecomposition(int numRow, int numColumn, double* inputMatrix, double* Q, double* R, double* economy_Q, double* economy_R);
double frobeniusNorm(int numRow, int numColumn, double* matrix);
double* linearConvolution(int size1, double* x, int size2, double* h);
double* circularConvolution(int size1, double* x, int size2, double* h);