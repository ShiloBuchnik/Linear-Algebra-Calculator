#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "auxFuncs.h"

// Getting rid of those annoying negative zeros, and returning row rank *of an echelon form matrix*
int negativeZerosAndFindRank(int numRow, int numColumn, long double* matrix)
{
    if (!matrix) return 0;

    bool flag = 0;
    int nonZeroRows = 0;
    for (int i = 0; i < numRow; i++)
    {
        for (int j = 0; j < numColumn; j++)
        {
            if (areEqual(*(matrix + i*numColumn + j), 0.0)) *(matrix + i*numColumn + j) = 0.0;
            else flag = 1;
        }
        if (flag) nonZeroRows++;
        flag = 0;
    }

    return nonZeroRows;
}

/* This function takes the original matrix, its size, and a point the defines the wanted minor
It copies said minor into 'dupMatrix' */
static void getMinor(int n, int rowValue, int columnValue, long double* originMatrix, long double* dupMatrix)
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
long double deterCalc(int n, long double* originMatrix)
{
    if (n == 1) return *(originMatrix); // Base case, Det(a) = a

    long double* dupMatrix = (long double*) malloc( (n - 1) * (n - 1) * sizeof(long double)); // This array represents a minor
    long double sum = 0, minor;

    for (int j = 0; j < n; j++)
    {
        getMinor(n, 0, j, originMatrix, dupMatrix); // Running on 1st row to calculate determinant
        minor = deterCalc(n - 1, dupMatrix); // Calculating the minor

        sum = sum + (long double) pow(-1, j) * *(originMatrix + j) * minor; // Determinants got a +-+- sign pattern, hence (-1)^j
    }

    free(dupMatrix);

    return sum;
}

void findAndPrintAdjoint(int n, long double* matrix)
{
    if (n == 1)
    {
      printf("%lf\n\n", 1.0);
      return;
    }

    long double* dupArray = (long double*) malloc( (n - 1) * (n - 1) * sizeof(long double)); // This array stores the minor
    long double* adjoint = (long double*) malloc(n * n * sizeof(long double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            getMinor(n, j, i, matrix, dupArray); // We want the 'ji' minor

            *(adjoint + i*n + j) = deterCalc(n - 1, dupArray) * (long double) pow(-1, j + i);
        }
    }

    negativeZerosAndFindRank(n, n, adjoint);
    displayMatrix(n, n, adjoint);
    free(dupArray);
    free(adjoint);
}

/* This function multiplies a row or column, in a 1D matrix, by a scalar
We pass a 'start' and 'end' pointers, and a char saying if it's a row ('r') or column ('c').
Based on that we know what arithmetic needs to be done */
static void rowOrColumnMultiply(long double* start, long double* end, long double scalar, int numColumn, char row_or_column)
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

The 'isPermutationIdentity' flag is set to 1 if there's no need for a permutation

*Remember that 'matrix[i][j] = matrix + i*numColumn + j' */
void gaussJordanAndFindInverse(int numRow, int numColumn, long double* inputMatrix, long double* inverseMatrix,
                               long double* permutation, bool* isPermutationIdentity, long double* L, long double* U)
{
    /* 'j' iterates over columns, 'i' iterates over rows, 'z' stores the current row in the process (pivot's row)
    'k' (row) and 'p' (column) are used to subtract the entire pivot's row from all the rows below it */
    int i = 0, j = 0, k = 0, p = 0, z = 0;
    long double pivot, scalar;

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

    if (U) matrixCopy(numRow, numColumn, inputMatrix, U);

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

// This function multiplies vectors element-wise (dot product)
static long double dotProduct(int size, long double* vector1, long double* vector2)
{
    long double sum = 0;
    for (int i = 0; i < size; i++) sum += *(vector1 + i) * *(vector2 + i);

    return sum;
}

/* THIS FUNCTION ALLOCATES MEMORY
This function multiplies 'matrix1' and 'matrix2' (in this order)
Note that 'numColumn1 == numRow2', so no need for another variable
Note that order of matrices passed matters */
long double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, long double* matrix1, long double* matrix2)
{
    long double* product = (long double*) malloc(numRow1 * numColumn2 * sizeof(long double));

    for (int i = 0; i < numRow1; i++)
    {
        for (int j = 0; j < numColumn2; j++)
        {
            long double* row_vector = separateRowOrColumn(numRow1, numColumn1, i, matrix1, 'r');
            long double* column_vector = separateRowOrColumn(numColumn1, numColumn2, j, matrix2, 'c');
            *(product + i*numColumn2 + j) = dotProduct(numColumn1, row_vector, column_vector);

            free(row_vector);
            free(column_vector);
        }
    }

    return product;
}

/* This function takes the permutation needed, and finds the LU decomposition of a given matrix
Use it as follows:
Call to 'GaussJordan' once to get the permutation and the flag telling us if no permutation is needed; with 'L' and 'U' as NULL
Now multiply the permutation with the matrix, and pass it to 'GaussJordan' to get 'L' and 'U'; with now permutation and the flag as NULL

Which makes sense, since in the algorithm we need to first find the permutation in elimination, and then eliminate again */
void LUDecomposition(int numRow, int numColumn, long double* inputMatrix, long double* Permutation, long double* L, long double* U)
{
    inputMatrix = matrixMultiplier(numRow, numRow, numColumn, Permutation, inputMatrix);
    gaussJordanAndFindInverse(numRow, numColumn, inputMatrix, NULL, NULL, NULL, L, U);
    free(inputMatrix);

    // Running on the 'L' matrix
    for (int j = 0; j < numRow; j++) // L is square, ergo we run up to 'numRow', not 'numColumn'
    {
        if (areEqual(*(L + j*numRow + j), 1.0)) continue;// Optimization: if the diagonal element is 1, we're moving on to the next one.
        long double scalar = *(L + j*numRow + j);

        rowOrColumnMultiply(L + j * numRow + j, L + numRow * numRow + j, 1 / scalar, numRow, 'c');
        rowOrColumnMultiply(U + j * numColumn + j, U + (j + 1) * numColumn, scalar, numColumn, 'r');
    }
}

/* This function gets the 'D' and 'V' matrices from the 'U' matrix (from 'LU' decomposition). Thus, we get an 'LVD' decomposition
This decomposition only exists when the matrix is square, so 'U' is square, so no need for 'numColumn' */
void findDVFromUAndPrint(int numRow, long double* U)
{
    long double* D = (long double*) malloc(numRow * numRow * sizeof(long double));
    setIdentityMatrix(numRow, D);
    long double* V = (long double*) malloc(numRow * numRow * sizeof(long double));

    for (int i = 0; i < numRow; i++)
    {
        long double pivot = *(U + i*numRow + i);
        *(D + i*numRow + i) = pivot;

        for (int j = 0; j < numRow; j++) *(V + i*numRow + j) = *(U + i*numRow + j) / pivot;
    }

    printf("D:\n");
    displayMatrix(numRow, numRow, D);
    printf("V:\n");
    displayMatrix(numRow, numRow, V);
}

/* THIS FUNCTION ALLOCATES MEMORY FOR 'transpose'
This function returns the transposition of 'inputMatrix'.
Remember that 'numRow' and 'numColumn' are switched for the transposition! */
static long double* transpose(int numRow, int numColumn, long double* inputMatrix)
{
    long double* transpose = (long double*) malloc(numRow * numColumn * sizeof(long double));

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

static bool areColumnsIndependent(int numRow, int numColumn, long double* inputMatrix)
{
    bool columnsIndependent = 0;
    long double* inputMatrixCopy = (long double*) malloc(numRow * numColumn * sizeof(long double));
    matrixCopy(numRow, numColumn, inputMatrix, inputMatrixCopy);

    if (numColumn <= numRow) // If there's more columns than rows, they're necessarily dependent, and there's no need to find the rank
    {
        gaussJordanAndFindInverse(numRow, numColumn, inputMatrixCopy, NULL, NULL, NULL, NULL, NULL);
        int rank = negativeZerosAndFindRank(numRow, numColumn, inputMatrixCopy);
        if (numColumn == rank) columnsIndependent = 1;
    }
    free(inputMatrixCopy);

    return columnsIndependent;
}

/* THIS FUNCTION MIGHT ALLOCATE MEMORY FOR 'result'
This function takes a flag and subtracts two vectors of same size
If alloc == 'y', then it allocates the result and returns it
If alloc == 'n', then it puts the result in 'vector1' */
static long double* vectorSubtraction(int size, long double* vector1, long double* vector2, char alloc)
{
    if (alloc != 'y' && alloc != 'n') return NULL;

    if (alloc == 'y')
    {
        long double* result = (long double*) malloc(size * sizeof(long double));
        for (int i = 0; i < size; i++) *(result + i) = *(vector1 + i) - *(vector2 + i);

        return result;
    }
    else
    {
        for (int i = 0; i < size; i++) *(vector1 + i) -= *(vector2 + i); // In this case alloc == 'n'
        return NULL;
    }
}

static long double norm2(int size, long double* vector)
{
    long double norm = 0;
    for (int i = 0; i < size; i++) norm += (long double) pow(*(vector + i), 2);
    norm = (long double) sqrt(norm);

    return norm;
}

static void multiplyVectorByScalar(int size, long double scalar, long double* vector)
{
    for (int i = 0; i < size; i++) *(vector + i) *= scalar;
}

/* This function initiate Graham-Schmidt algorithm on given matrix's *columns*, and stores the result in Q's columns
If the matrix is "tall" (more rows than columns) - the function adds to the input some random vectors until it can output 'numRow' vectors
(this is done so that Q can be square)

if on the way we run into a vector that is dependent on the ones before it - we replace it with a random vector and carry on */
static void grahamSchmidtForQR(int numRow, int numColumn, long double* matrix, long double* Q)
{
    /* Instead of taking a random vector when we run into a dependent vector in the process - we take a standard one, to make the numbers prettier
    this counter tracks which one to take now. First we take (1, 0, 0,...), then (0, 1, 0,...) and so on
    We are *assured* that we'll finish the process before we've run out of vectors;
    since it's not possible for *all* of the standard vectors of 'numRow' size to be dependent on a group of size less than 'numRow',
    since that means that said group is spanning a vector space of higher dimension that the group's size */
    int dependent_counter = 0;
    bool get_standard = 0;

    for (int j = 0; j < numRow; j++)
    {
        long double* w = NULL;
        if (j < numColumn && !get_standard) w = separateRowOrColumn(numRow, numColumn, j, matrix, 'c');
        else // We're getting here only if it's a tall matrix, so that we generate more input vectors
        {
            w = getStandardVector(numRow, dependent_counter);
            dependent_counter++;
            get_standard = 0;
        }

        for (int k = 0; k < j; k++) // "peeling" the vector
        {
            long double* u = separateRowOrColumn(numRow, numRow, k, Q, 'c');
            long double scalar = dotProduct(numRow, u, w);

            multiplyVectorByScalar(numRow, scalar, u);
            vectorSubtraction(numRow, w, u, 'n');

            free(u);
        }

        /* We find out if 'w' is 0 through the theorem: 'a matrix is 0 iff its rank is 0'
        this tells us it's dependent on the vectors before it, so we replace it with a random vector and do the current iteration again */
        if (!negativeZerosAndFindRank(numRow, 1, w))
        {
            j--;
            get_standard = 1;
        }
        else
        {
            long double norm = norm2(numRow, w);
            multiplyVectorByScalar(numRow, 1 / norm, w);
            insertColumn(numRow, numRow, j, Q, w, 'o');
        }

        free(w);
    }
}

/* This function finds the QR decomposition of any matrix and puts it in 'Q' and 'R'
If the matrix is "tall", it even finds the economy QR decomposition and puts it in 'economy_Q' and 'economy_R'.
'economy' version basically means that we stop Graham-Schmidt after we get 'numColumn' vectors, and not try to make 'Q' a square matrix
Its columns are still orthonormal, the matrix itself is just not orthonormal, since it's not square */
void QRDecomposition(int numRow, int numColumn, long double* inputMatrix, long double* Q, long double* R, long double* economy_Q, long double* economy_R)
{
    grahamSchmidtForQR(numRow, numColumn, inputMatrix, Q); // Getting Q
    negativeZerosAndFindRank(numRow, numRow, Q);
    long double* Q_T = transpose(numRow, numRow, Q);
    long double* temp = matrixMultiplier(numRow, numRow, numColumn, Q_T, inputMatrix); // Now we got R, since Q^T * A = R
    matrixCopy(numRow, numColumn, temp, R);
    negativeZerosAndFindRank(numRow, numColumn, R);
    free(temp);

    // The matrix is "tall", so there is also the economy version.
    // We also produce it for square matrices, although it's identical to the regular 'Q' and 'R', to fit better with 'leastSquares()'
    if (numColumn <= numRow)
    {
        truncateMatrix(numRow, Q, numRow, numColumn, economy_Q); // Getting economy Q
        negativeZerosAndFindRank(numRow, numColumn, economy_Q);
        truncateMatrix(numColumn, R, numColumn, numColumn, economy_R); // Getting economy R
        negativeZerosAndFindRank(numColumn, numColumn, economy_R);
    }

    free(Q_T);
}

/* THIS FUNCTION ALLOCATES MEMORY FOR 'x' AND 'p_x'
This function solves a system of equations through the least squares algorithm,
only instead of solving the problem the "regular" way, we utilize the QR decomposition:
Ax=b -> QRx=b -> Rx=Q^Tb -> x=R^-1Q^Tb

*We can be sure 'R' has an inverse, since we only call 'QRDecomposition()' only when A's columns are independent;
and A's columns are independent iff R's columns are independent

If there's a *single* minimum point, it stores it in 'x' (note it's a long double pointer), and returns the remainder vector
(When the minimum point is a true solution, the remainder is 0, obviously)
If there's infinite solutions, it returns NULL */
long double* leastSquares(int numRow, int numColumn, long double* A, long double** x, long double* b)
{
    if (!areColumnsIndependent(numRow, numColumn, A)) return NULL;

    long double* Q = (long double*) malloc(numRow * numRow * sizeof(long double));
    long double* R = (long double*) malloc(numRow * numColumn * sizeof(long double));
    long double* economy_Q = (long double*) malloc(numRow * numColumn * sizeof(long double));
    long double* economy_R = (long double*) malloc(numColumn * numColumn * sizeof(long double));
    QRDecomposition(numRow, numColumn, A, Q, R, economy_Q, economy_R);
    free(Q);
    free(R);

    long double* Q_T = transpose(numRow, numColumn, economy_Q); // Now we have Q^T
    long double* R_inverse = (long double*) malloc(numColumn * numColumn * sizeof(long double));
    setIdentityMatrix(numColumn, R_inverse);
    gaussJordanAndFindInverse(numColumn, numColumn, economy_R, R_inverse, NULL, NULL, NULL, NULL); // Now we have R^-1

    long double* Q_T_times_b = matrixMultiplier(numColumn, numRow, 1, Q_T, b);
    *x = matrixMultiplier(numColumn, numColumn, 1, R_inverse, Q_T_times_b); // This is the minimal point

    long double* A_times_x = matrixMultiplier(numRow, numColumn, 1, A, *x);
    long double* p_x = vectorSubtraction(numRow, A_times_x, b, 'y'); // This is the function p(x) (the remainder function) at the minimal point

    free(economy_R);
    free(economy_Q);
    free(Q_T);
    free(R_inverse);
    free(Q_T_times_b);
    free(A_times_x);

    return p_x;
}



