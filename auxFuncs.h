#include <stdbool.h>

#define min(a, b) (a <= b) ? a : b
#define max(a, b) (a > b) ? a : b

int enableColor(void);
void swap(long double* num1, long double* num2);
int areEqual(long double a, long double b);
void positiveIntVerify(int* input, int* arr, int arrLen);
void doubleVerify(int numRow, int numColumn, long double* matrix, bool print);
void getRowAndColumnNum(int* numRow, int* numColumn);
void getVectorSize(int* size);
void displayMatrix(int numRow, int numColumn, long double* matrix);
void displaySolution(int size, long double* solution);
void setIdentityMatrix(int n, long double* matrix);
void matrixCopy(int numRow, int numColumn, long double* srcMatrix, long double* destMatrix);
void switchRows(int row1, int row2, int numColumn, long double* matrix);
long double* insertColumn(int numRow, int numColumn, int column, long double* matrix, long double* columnVector, char extend_or_overwrite);
void truncateMatrix(int numColumn, long double* srcMatrix, int newNumRow, int newNumColumn, long double* destMatrix);
long double* getStandardVector(int size, int k);
long double* separateRowOrColumn(int numRow, int numColumn, int row_or_column, long double* matrix, char rc);
int printOpeningAndGetNum(void);
void printText(int num);