#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "auxFuncs.h"
#include "coreFuncs.h"

int main()
{
    enableColor();
    int mode = printOpeningAndGetNum();

    //printText(mode);

    switch(mode)
    {
    case 1: // Matrix calculator
    {
        int numRow, numColumn;

        getRowAndColumnNum(&numRow, &numColumn);
        bool isSquare = (numRow == numColumn);

        long double* inputMatrix = (long double*) malloc(numRow * numColumn * sizeof(long double));
        doubleVerify(numRow, numColumn, inputMatrix, 1);

        long double* echelonFormMatrix = (long double*) malloc(numRow * numColumn * sizeof(long double));
        matrixCopy(numRow, numColumn, inputMatrix, echelonFormMatrix);


        long double* inverseMatrix = NULL;
        if (isSquare)
        {
            inverseMatrix = (long double*) malloc(numRow * numColumn * sizeof(long double));
            setIdentityMatrix(numRow, inverseMatrix); // We set it to identity matrix, and in the end of the process get the inverse
        }

        long double* permutationMatrix = (long double*) malloc(numRow * numRow * sizeof(long double));
        setIdentityMatrix(numRow, permutationMatrix); // We set it to identity matrix, and in the end of the process get the permutation
        bool isPermutationIdentity = 1;
        long double* L = (long double*) malloc(numRow * numRow * sizeof(long double));
        setIdentityMatrix(numRow, L);
        long double* U = (long double*) malloc(numRow * numColumn * sizeof(long double));

        long double* Q = (long double*) malloc(numRow * numRow * sizeof(long double));
        long double* R = (long double*) malloc(numRow * numColumn * sizeof(long double));
        long double* economy_Q = (long double*) malloc(numRow * numColumn * sizeof(long double));
        long double* economy_R = (long double*) malloc(numColumn * numColumn * sizeof(long double));
        QRDecomposition(numRow, numColumn, inputMatrix, Q, R, economy_Q, economy_R);

        printf("Your \x1b[96mmatrix\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, inputMatrix);


        GaussJordanAndFindInverse(numRow, numColumn, echelonFormMatrix, inverseMatrix, permutationMatrix, &isPermutationIdentity, NULL, NULL);

        long double* inputMatrixCopy = (long double*) malloc(numRow * numColumn * sizeof(long double));
        matrixCopy(numRow, numColumn, inputMatrix, inputMatrixCopy);
        LUDecomposition(numRow, numColumn, inputMatrixCopy, permutationMatrix, L, U);
        free(inputMatrixCopy);

        int rank = negativeZerosAndFindRank(numRow, numColumn, echelonFormMatrix);
        negativeZerosAndFindRank(numRow, numColumn, inverseMatrix);

        printf("The matrix in \x1b[92mreduced row echelon form\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, echelonFormMatrix);

        printf("The \x1b[93mrank\x1b[0m is: %d\n\n", rank);

        long double deter;
        if (isSquare)
        {
            deter = deterCalc(numColumn, inputMatrix);

            if (!areEqual(deter, 0.0)) // If deter != 0, then the matrix has an inverse
            {
                printf("The \x1b[94minverse matrix\x1b[0m is:\n");
                displayMatrix(numRow, numColumn, inverseMatrix);
            }
            else printf("This matrix is \x1b[94msingular\x1b[0m\n\n");

            printf("The \x1b[91mdeterminant\x1b[0m is: %Lf\n\n", deter);

            printf("The \x1b[95madjoint matrix\x1b[0m is:\n\n");
            findAndPrintAdjoint(numColumn, inputMatrix);
       }
        else
            printf("Non-square matrices don't have an \x1b[94minverse matrix\x1b[0m, "
            "nor a \x1b[91mdeterminant\x1b[0m, nor an \x1b[95madjoint\x1b[0m\n\n\n");


        printf("The LU decomposition is:\n");
        printf("L:\n");
        displayMatrix(numRow, numRow, L);
        printf("U:\n");
        displayMatrix(numRow, numColumn, U);
        if (isSquare && !areEqual(deter, 0.0) && isPermutationIdentity)
            printf("This matrix is regular (square, invertible, and has an LU decomposition without a permutation)\n");
        else
        {
            printf("With a permutation:\n");
            displayMatrix(numRow, numRow, permutationMatrix);
        }

        if (isSquare && !areEqual(deter, 0.0))
        {
            printf("The LDV decomposition is:\n");
            printf("L:\n");
            displayMatrix(numRow, numRow, L);
            findDVFromUAndPrint(numRow, U);
        }
        else printf("Non-square or singular matrices don't have an LDV decomposition\n\n");


        printf("The QR decomposition is:\n");
        printf("Q:\n");
        displayMatrix(numRow, numRow, Q);
        printf("R:\n");
        displayMatrix(numRow, numColumn, R);

        if (numColumn < numRow)
        {
            printf("This is a \"tall\" matrix, so it has an \"economy\" QR decomposition as well:\n");
            printf("Economy Q:\n");
            displayMatrix(numRow, numColumn, economy_Q);
            printf("Economy R:\n");
            displayMatrix(numColumn, numColumn, economy_R);
        }

        free(inputMatrix);
        free(echelonFormMatrix);
        free(inverseMatrix);
        free(permutationMatrix);
        free(L);
        free(U);
        free(Q);
        free(economy_Q);
        free(R);
        free(economy_R);
        break;
    }

    case 2: // System of equations calculator.
    {
        int numRow, numColumn;

        printf("Enter variable coefficients matrix's (\x1b[31mA\x1b[0m) dimensions:\n");
        getRowAndColumnNum(&numRow, &numColumn);
        //while (getchar() != '\n');
        printf("\n");

        long double* A = (long double*) malloc(numRow * numColumn * sizeof(long double));
        doubleVerify(numRow, numColumn, A, 1);

        printf("Enter the free coefficients (\x1b[34mb\x1b[0m):\n");
        long double* b = (long double*) malloc(numRow * sizeof(long double));
        doubleVerify(numRow, 1, b, 0);

        long double* A_star = insertColumn(numRow, numColumn, numColumn, A, b, 'e');
        printf("The scalars of your system are:\n");
        displayMatrix(numRow, numColumn + 1, A_star);

        long double* A_copy = (long double*) malloc(numRow * numColumn * sizeof(long double)); // We want 'A' to remain the same, so we copy to find rank
        matrixCopy(numRow, numColumn, A, A_copy);
        GaussJordanAndFindInverse(numRow, numColumn, A_copy, NULL, NULL, NULL, NULL, NULL);
        int A_rank = negativeZerosAndFindRank(numRow, numColumn, A_copy);
        free(A_copy);

        GaussJordanAndFindInverse(numRow, numColumn + 1, A_star, NULL, NULL, NULL, NULL, NULL); // We don't care about A_star besides its rank
        int A_star_rank = negativeZerosAndFindRank(numRow, numColumn + 1, A_star);
        free(A_star);

        long double* x = NULL;
        long double* remainder = leastSquares(numRow, numColumn, A, &x, b);


        /* If 'leastSquares' returned something (the columns are independent) - no problem.
        Either there's a single true solution, or a single minimal point, in least squares sense

        But if it returned 'NULL' (the columns are dependent) - that's where our 'A' and 'A_star' ranks come in place.
        If the ranks are equal, there are infinite *true* solutions
        If the ranks are not equal, there are infinite minimal points (not true solutions), in least squares sense */
        if (remainder)
        {
            negativeZerosAndFindRank(numColumn, 1, x);
            negativeZerosAndFindRank(numRow, 1, remainder);

            // We find out if the remainder is non-0 through the theorem: 'a matrix is 0 iff its rank is 0'
            if (!negativeZerosAndFindRank(numRow, 1, remainder))
            {
                printf("The solution vector \x1b[92mx\x1b[0m is: ");
                displaySolution(numColumn, x);
            }
            else
            {
                printf("There isn't a solution. The minimal solution vector \x1b[92mx\x1b[0m in 'Least Square' sense is:\n");
                displaySolution(numColumn, x);
                printf("With a \x1b[96mremainder\x1b[0m vector:\n");
                displaySolution(numRow, remainder);
            }
        }
        else
        {
            if (A_rank == A_star_rank) printf("There are infinite solutions\n");
            else printf("There isn't a solution, but there are infinite minimal vectors in 'Least Square' sense\n");
        }

        free(A);
        free(b);
        free(x);
        free(remainder);
        break;
    }

    case 3: // Matrix multiplier
    {
        int numRow1, numColumn1, numRow2, numColumn2;

        printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
        getRowAndColumnNum(&numRow1, &numColumn1);

        printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
        getRowAndColumnNum(&numRow2, &numColumn2);

        while (numColumn1 != numRow2)
        {
            printf("Number of columns of \x1b[91mmatrix 1\x1b[0m must be equal to number of rows of \x1b[94mmatrix 2\x1b[0m. Please try again:\n");

            printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
            getRowAndColumnNum(&numRow1, &numColumn1);
            printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
            getRowAndColumnNum(&numRow2, &numColumn2);
        }

        printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
        long double* matrix1 = (long double*) malloc(numRow1 * numColumn1 * sizeof(long double));
        doubleVerify(numRow1, numColumn1, matrix1, 1);

        printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
        long double* matrix2 = (long double*) malloc(numRow2 * numColumn2 * sizeof(long double));
        doubleVerify(numRow2, numColumn2, matrix2, 1);

        long double* product = matrixMultiplier(numRow1, numColumn1, numColumn2, matrix1, matrix2);
        negativeZerosAndFindRank(numRow1, numColumn2, product);
        printf("The \x1b[92mproduct\x1b[0m is:\n");
        displayMatrix(numRow1, numColumn2, product);

        free(matrix1);
        free(matrix2);
        free(product);
        break;
    }
    }

    system("pause");

    return 0;
}