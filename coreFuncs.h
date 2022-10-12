int negativeZerosAndFindRank(int numRow, int numColumn, long double* matrix);
long double deterCalc(int n, long double* originMatrix);
void findAndPrintAdjoint(int n, long double* matrix);
void GaussJordanAndFindInverse(int numRow, int numColumn, long double* inputMatrix, long double* inverseMatrix,
                               long double* permutation, bool* isPermutationIdentity, long double* L, long double* U);
long double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, long double* matrix1, long double* matrix2);
void LUDecomposition(int numRow, int numColumn, long double* inputMatrix, long double* Permutation, long double* L, long double* U);
void findDVFromUAndPrint(int numRow, long double* U);
long double* leastSquares(int numRow, int numColumn, long double* A, long double** x, long double* b);
void QRDecomposition(int numRow, int numColumn, long double* inputMatrix, long double* Q, long double* R, long double* economy_Q, long double* economy_R);