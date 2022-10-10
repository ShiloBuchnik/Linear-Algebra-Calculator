#include <stdbool.h>

int enableColor(void);
void swap(double* num1, double* num2);
int areEqual(double a, double b);
void positiveIntVerify(int* input, int* arr, int arrLen);
void doubleVerify(int numRow, int numColumn, double* matrix, bool print);
void getRowAndColumnNum(int* numRow, int* numColumn);
void displayMatrix(int numRow, int numColumn, double* matrix);
void displaySolution(int size, double* solution);
void setIdentityMatrix(int n, double* matrix);
void matrixCopy(int numRow, int numColumn, double* srcMatrix, double* destMatrix);
void switchRows(int row1, int row2, int numColumn, double* matrix);
double* insertColumn(int numRow, int numColumn, int column, double* matrix, double* columnVector);
int printOpeningAndGetNum(void);
void printText(int num);