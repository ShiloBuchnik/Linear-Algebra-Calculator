int negativeZerosAndFindRank(int numRow, int numColumn, double *matrix);
double deterCalc(int n, double *originMatrix);
void findAndPrintAdjoint(int n, double *matrix);
int CramersRule(int size, double *solution, double *scalarMatrix, double *bColumn);
void GaussJordanAndFindInverse(int numRow, int numColumn, double *inputMatrix, double *inverseMatrix,
                               double *permutation, double *L, double *U);
double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, double *matrix1, double *matrix2);