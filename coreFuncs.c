#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "auxFuncsDec.h"

// Getting rid of those annoying negative zeros, and returning row rank of an echelon form matrix
int negativeZerosAndFindRank(int numRow, int numColumn, double* matrix)
{
    if (!matrix) return 0;

    bool flag = 0;
    int nonZeroRows = 0;
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < numColumn; j++)
        {
            if (*(matrix + i*numColumn + j) == 0.0) *(matrix + i*numColumn + j) = 0.0;
            else flag = 1;
        }
        if (flag) nonZeroRows++;
        flag = 0;
    }

    return nonZeroRows;
}

/* This function takes the original matrix, its size, and a point the defines the wanted minor
It copies said minor into 'dupMatrix' */
static void getMinor(int n, int rowValue, int columnValue, double* originMatrix, double* dupMatrix)
{
    int i, j; // For running on originArray
    int x = 0; // For running on dupArray

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == rowValue || j == columnValue) continue;

            *(dupMatrix + x) = *(originMatrix + i*n + j);
            x++;
        }
    }
}

/* Calculating determinant recursively (with minors), input array is a 1D matrix, 'n' is size
Yes, I know we can utilize the GaussJordanElimination for that, or even compute it with the LU decomposition,
but isn't recursion fun? */
double deterCalc(int n, double* originMatrix)
{
    if (n == 1) return *(originMatrix); // Base case, Det(a) = a

    double* dupMatrix = (double*) malloc( (n - 1) * (n - 1) * sizeof(double)); // This array represents a minor
    double sum = 0, minor;

    for (int j = 0; j < n; j++)
    {
        getMinor(n, 0, j, originMatrix, dupMatrix); // Running on 1st row to calculate determinant
        minor = deterCalc(n - 1, dupMatrix); // Calculating the minor

        sum = sum + pow(-1, j) * *(originMatrix + j) * minor; // Determinants got a +-+- sign pattern, hence (-1)^j
    }

    free(dupMatrix);

    return sum;
}

void findAndPrintAdjoint(int n, double* matrix)
{
    if (n == 1)
    {
      printf("%lf\n\n", 1.0);
      return;
    }

    double* dupArray = (double*) malloc( (n - 1) * (n - 1) * sizeof(double)); // This array stores the minor
    double* adjoint = (double*) malloc(n * n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            getMinor(n, j, i, matrix, dupArray); // We want the 'ji' minor

            *(adjoint + i*n + j) = deterCalc(n - 1, dupArray) * pow(-1, j + i);
        }
    }

    negativeZerosAndFindRank(n, n, adjoint);
    displayMatrix(n, n, adjoint);
    free(dupArray);
    free(adjoint);
}

/* This function takes the 'A' and 'b' in Ax=b, and stores the 'x', if there is a solution, in 'solution'
The calculation is done with Cramer's rule */
bool CramersRule(int size, double* scalarMatrix, double* solution, double* bColumn)
{
    double matrixDeter = deterCalc(size, scalarMatrix);

    if (areEqual(matrixDeter, 0.0)) return 0; // No single solution

    for (int i = 0; i < size; i++)
    {
        switchColumns(size, size, i, scalarMatrix, bColumn); // Switching to calculate determinant ('size' represents the b vector in Ax=b)

        *(solution + i) = deterCalc(size, scalarMatrix) / matrixDeter; // The actual formula

        switchColumns(size, size, i, scalarMatrix, bColumn); // Switching back
    }

    return 1; // A single solution found
}

/* This function multiplies a row or column, in a 1D matrix, by a scalar
We pass a 'start' and 'end' pointers, and a char saying if it's a row ('r') or column ('c').
Based on that we know what arithmetic needs to be done */
static void rowOrColumnMultiply(double* start, double* end, double scalar, int numColumn, char row_or_column)
{
    int interval = (row_or_column == 'r' ? 1 : numColumn);

    while (start != end)
    {
        *start *= scalar;
        start += interval;
    }
}


/* 'input' is passed to perform elimination on, and an identity matrix is passed as 'inverse'
Through elimination, as 'input' gets echelon form, 'inverse' becomes the inverse matrix (if 'input' is invertible)
If we don't want to find the inverse, we pass NULL

Since LU decomposition relies on Gauss Elimination (echelon form), it makes sense to do it in this function as well
BUT, this function does GaussJordan Elimination (*reduced* echelon form), so we had to add an addition to the algorithm:
When we divide the row by pivot (which doesn't happen on Gauss elim'), we multiply the corresponding pivot in 'L'
Then, in the end, we run on 'L's diagonal, and:
For every element 'x' which is not equal to 1, we divide its column by 'x',
and then go to 'U' and multiply the corresponding *row* by same 'x'. Think for yourself why we can do this.

*Remember that 'matrix[i][j] = matrix + i*numColumn + j' */
void GaussJordanAndFindInverse(int numRow, int numColumn, double* inputMatrix, double* inverseMatrix,
                               double* permutation, int* isPermutationIdentity, double* L, double* U)
{
    /* 'j' iterates over columns, 'i' iterates over rows, 'z' stores the current row in the process (pivot's row)
    'k' (row) and 'p' (column) are used to subtract the entire pivot's row from all the rows below it */
    int i = 0, j = 0, k = 0, p = 0, z = 0;
    double pivot, scalar;

    for (j = 0; j < numColumn; j++) // For the bottom triangle of the matrix. After this, the matrix is in echelon form
    {
        for (i = z; i < numRow; i++)
        {
            if (*(inputMatrix + i*numColumn + j) == 0) continue; // Searching for a non-zero pivot in the given column

            if (i != z) // In this case, our pivot has zeroed out (can't happen when finding LU), and we need to switch rows for a new pivot.
            {
                switchRows(i, z, numColumn, inputMatrix);
                if (inverseMatrix) switchRows(i, z, numColumn, inverseMatrix);

                if (permutation) switchRows(i, z, numColumn, permutation); // We need to "write" the switch in the permutation
                if (isPermutationIdentity) *isPermutationIdentity = 0;
            }

            pivot = *(inputMatrix + z * numColumn + j);
            if (L) *(L + z * numRow + z) *= pivot;
            rowOrColumnMultiply(inputMatrix + z*numColumn + 0, inputMatrix + z*numColumn + numColumn, 1 / pivot, numColumn, 'r');
            if (inverseMatrix) rowOrColumnMultiply(inverseMatrix + z*numColumn + 0, inverseMatrix + z*numColumn + numColumn, 1 / pivot, numColumn, 'r');
            // Dividing the row by leading entry (pivot)

            for (k = z + 1; k < numRow; k++) // Subtracting pivot's row ('z's row) from the rows below it
            {
                scalar = *(inputMatrix + k * numColumn + j);
                if (L) *(L + k * numRow + z) = scalar;
                for (p = 0; p < numColumn; p++) // Note to self: can't start from 'z' instead of 0, since we have to calculate inverse as well
                {
                   *(inputMatrix + k*numColumn + p) -= scalar * *(inputMatrix + z * numColumn + p);
                   if (inverseMatrix) *(inverseMatrix + k*numColumn + p) -= scalar * *(inverseMatrix + z * numColumn + p);

                   if (areEqual(*(inputMatrix + k*numColumn + p), 0.0)) *(inputMatrix + k*numColumn + p) = 0; // Handling rounding errors
                   if (inverseMatrix && areEqual(*(inverseMatrix + k*numColumn + p), 0.0)) *(inverseMatrix + k*numColumn + p) = 0;
                }
            }
            z++;
            break; // We've zeroed all elements below pivot, and we're done with this column
        }
    }

    if (U) arrayCopy(numRow, numColumn, inputMatrix, U);

    for (i = numRow - 1; i >= 0; i--) // For the upper triangle of the matrix. After this, the matrix is in *reduced* echelon form
    {
        for (j = 0; j < numColumn; j++)
        {
            if (*(inputMatrix + i*numColumn + j) == 0) continue; // Searching for the leading element in the row

            for (k = i - 1; k >= 0; k--) // Subtracting pivot's row ('i's row) from rows above
                {
                    scalar = *(inputMatrix + k * numColumn + j);
                    for (p = 0; p < numColumn; p++)
                    {
                        *(inputMatrix + k*numColumn + p) -= scalar * *(inputMatrix + i * numColumn + p);
                        if (inverseMatrix) *(inverseMatrix + k*numColumn + p) -= scalar * *(inverseMatrix + i * numColumn + p);
                    }
                }
                break; // We've zeroed all elements above pivot, and we're done with this column
        }
    }
}

static double dotProduct(int column1, int column2, double* vector1, double* vector2)
{
    double sum = 0;
    for (int i = 0; i < column1; i++) sum += *(vector1 + i) * *(vector2 + i*column2);

    return sum;
}

/* THIS FUNCTION ALLOCATES MEMORY
This function multiplies 'matrix1' and 'matrix2' (in this order)
Note that 'numColumn1 == numRow2', so no need for another variable
Note that order of matrices passed matters */
double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, double* matrix1, double* matrix2)
{
    double* product = (double*) malloc(numRow1 * numColumn2 * sizeof(double));

    for (int i = 0; i < numRow1; i++)
    {
        for (int j = 0; j < numColumn2; j++)
        {
            *(product + i*numColumn2 + j) = dotProduct(numColumn1, numColumn2, matrix1 + numColumn1*i, matrix2 + j);
        }
    }

    return product;
}

/* This function takes the permutation we acquired from 'GaussJordanAndFindInverse', and finds the LU decomposition of a given matrix
We call to 'GaussJordan' once to get the permutation, and then pass it to this function, and call 'GaussJordan' again to get 'L' and 'U' */
void LUDecomposition(int numRow, int numColumn, double* inputMatrix, double* Permutation, double* L, double* U)
{
    inputMatrix = matrixMultiplier(numRow, numRow, numColumn, Permutation, inputMatrix);
    GaussJordanAndFindInverse(numRow, numColumn, inputMatrix, NULL, NULL, NULL, L, U);
    free(inputMatrix);

    // Running on the 'L' matrix
    for (int j = 0; j < numRow; j++) // L is square, ergo we run up to 'numRow', not 'numColumn'
    {
        if (areEqual(*(L + j*numRow + j), 1.0)) continue;// Optimization: if the diagonal element is 1, we're moving on to the next one.
        double scalar = *(L + j*numRow + j);

        rowOrColumnMultiply(L + j * numRow + j, L + numRow * numRow + j, 1 / scalar, numRow, 'c');
        rowOrColumnMultiply(U + j * numColumn + j, U + (j + 1) * numColumn, scalar, numColumn, 'r');
    }
}

/* This function gets the 'D' and 'V' matrices from the 'U' matrix (from 'LU' decomposition). Thus, we get an 'LVD' decomposition
This decomposition only exists when the matrix is square, so 'U' is square, so no need for 'numColumn' */
void findDVFromUAndPrint(int numRow, double* U)
{
    double* D = (double*) malloc(numRow * numRow * sizeof(double));
    setIdentityMatrix(numRow, D);
    double* V = (double*) malloc(numRow * numRow * sizeof(double));

    for (int i = 0; i < numRow; i++)
    {
        double pivot = *(U + i*numRow + i);
        *(D + i*numRow + i) = pivot;

        for (int j = 0; j < numRow; j++) *(V + i*numRow + j) = *(U + i*numRow + j) / pivot;
    }

    printf("D:\n");
    displayMatrix(numRow, numRow, D);
    printf("V:\n");
    displayMatrix(numRow, numRow, V);
}

/* THIS FUNCTION ALLOCATES MEMORY
This function returns the transposition of 'inputMatrix'.
Remember that 'numRow' and 'numColumn' are switched for the transposition! */
static double* transpose(int numRow, int numColumn, double* inputMatrix)
{
    double* transpose = (double*) malloc(numRow * numColumn * sizeof(double));

    // Running on transpose, so 'numRow' and 'numColumn' roles are switched
    for (int i = 0; i < numColumn; i++)
    {
        for (int j = 0; j < numRow; j++)
        {
            *(transpose + i*numRow + j) = *(inputMatrix + j*numColumn + i);
        }
    }

    return transpose;
}

static bool areColumnsIndependent(int numRow, int numColumn, double* inputMatrix)
{
    bool columnsIndependent = 0;
    double* inputMatrixCopy = (double*) malloc(numRow * numColumn * sizeof(double));
    arrayCopy(numRow, numColumn, inputMatrix, inputMatrixCopy);

    if (numColumn <= numRow) // If there's more columns than rows, they're necessarily dependent, and there's no need to find the rank
    {
        GaussJordanAndFindInverse(numRow, numColumn, inputMatrixCopy, NULL, NULL, NULL, NULL, NULL);
        int rank = negativeZerosAndFindRank(numRow, numColumn, inputMatrixCopy);
        if (numColumn == rank) columnsIndependent = 1;
    }
    free(inputMatrixCopy);

    return columnsIndependent;
}


/* THIS FUNCTION ALLOCATES MEMORY
This function returns the pseudo inverse of the 'inputMatrix', *if and only if the columns are independent*, and NULL otherwise
Although a pseudo inverse exists for every matrix, when the columns are independent, there's a simple formula to calculate it;
and in Least Squares it's all we care about

Remember that 'numRow' and 'numColumn' are switched for the pseudo inverse! */
static double* pseudoInverseCalc(int numRow, int numColumn, double* A)
{
    if (!areColumnsIndependent(numRow, numColumn, A)) return NULL;

    double* A_T = transpose(numRow, numColumn, A); // A^T
    double* A_T_A = matrixMultiplier(numColumn, numRow, numColumn, A_T, A); // A^T * A
    double* A_T_A_inverse = (double*) malloc(numColumn * numColumn * sizeof(double)); // (A^T * A)^-1
    setIdentityMatrix(numColumn, A_T_A_inverse);
    GaussJordanAndFindInverse(numColumn, numColumn, A_T_A, A_T_A_inverse, NULL, NULL, NULL, NULL);

    double* pseudo_inverse = matrixMultiplier(numColumn, numColumn, numRow, A_T_A_inverse, A_T); // (A^T * A)^-1 * A^T

    free(A_T);
    free(A_T_A);
    free(A_T_A_inverse);
    return pseudo_inverse;
}

/* THIS FUNCTION ALLOCATES MEMORY
This function subtract two vectors of same size, and returns the result */
static double* vectorSubtraction(int size, double* vector1, double* vector2)
{
    double* result = (double*) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) *(result + i) = *(vector1 + i) - *(vector2 + i);

    return result;
}

/* THIS FUNCTION ALLOCATES MEMORY
This function solves a system of equations through the least squares algorithm
If there's a *single* minimum point, it stores it in 'x' (note it's a double pointer), and returns the remainder vector
(When the minimum point is a true solution, the remainder is 0, obviously)
If there's infinite solutions, it returns NULL */
double* leastSquares(int numRow, int numColumn, double* A, double** x, double* b)
{
    double* pseudo_inverse = pseudoInverseCalc(numRow, numColumn, A);
    if (!pseudo_inverse) return NULL;

    *x = matrixMultiplier(numColumn, numRow, 1, pseudo_inverse, b); // This is the minimal point

    double* temp = matrixMultiplier(numRow, numColumn, 1, A, *x);
    double* P_X = vectorSubtraction(numRow, temp, b); // This is the function's value ('p(x)') at minimal point

    free(pseudo_inverse);
    free(temp);

    return P_X;
}