int negativeZerosAndFindRank(int numRow, int numColumn, double* matrix);
double deterCalc(int n, double* originMatrix);
void findAndPrintAdjoint(int n, double* matrix);
int CramersRule(int size, double* scalarMatrix, double* solution, double* bColumn);
void GaussJordanAndFindInverse(int numRow, int numColumn, double* inputMatrix, double* inverseMatrix,
                               double* permutation, int* isPermutationIdentity, double* L, double* U);
double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, double* matrix1, double* matrix2);
void LUDecomposition(int numRow, int numColumn, double* inputMatrix, double* Permutation, double* L, double* U);
void findDVFromUAndPrint(int numRow, double* U);
double* pseudoInverseCalc(int numRow, int numColumn, double* A);