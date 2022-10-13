#include <stdbool.h>

#define min(a, b) (a <= b) ? a : b
#define max(a, b) (a > b) ? a : b

int enableColor(void);
void swap(double* num1, double* num2);
int areEqual(double a, double b);
void positiveIntVerify(int* input, int* arr, int arrLen);
void doubleVerify(int numRow, int numColumn, double* matrix, bool print);
void getRowAndColumnNum(int* numRow, int* numColumn);
void getVectorSize(int* size);
void displayMatrix(int numRow, int numColumn, double* matrix);
void displaySolution(int size, double* solution);
void setIdentityMatrix(int n, double* matrix);
void matrixCopy(int numRow, int numColumn, double* srcMatrix, double* destMatrix);
void switchRows(int row1, int row2, int numColumn, double* matrix);
double* insertColumn(int numRow, int numColumn, int column, double* matrix, double* columnVector, char extend_or_overwrite);
void truncateMatrix(int numColumn, double* srcMatrix, int newNumRow, int newNumColumn, double* destMatrix);
double* getStandardVector(int size, int k);
double* returnRowOrColumn(int numRow, int numColumn, int row_or_column, double* matrix, char rc);
int modulo(int n, int N);
int printOpeningAndGetNum(void);
void printText(int num);