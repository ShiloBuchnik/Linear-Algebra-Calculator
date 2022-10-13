#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "coreFuncs.h"
#include "auxFuncs.h"

int main()
{
    //setbuf(stdout, 0);
    enableColor();
    int mode = printOpeningAndGetNum();

    printText(mode);

    switch(mode)
    {
    case 1: // Matrix calculator
    {
        int numRow, numColumn;

        getRowAndColumnNum(&numRow, &numColumn);
        bool isSquare = (numRow == numColumn);

        double* inputMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
        doubleVerify(numRow, numColumn, inputMatrix, 1);

        double* echelonFormMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
        matrixCopy(numRow, numColumn, inputMatrix, echelonFormMatrix);


        double* inverseMatrix = NULL;
        if (isSquare)
        {
            inverseMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
            setIdentityMatrix(numRow, inverseMatrix); // We set it to identity matrix, and in the end of the process get the inverse
        }

        double* permutationMatrix = (double*) malloc(numRow * numRow * sizeof(double));
        setIdentityMatrix(numRow, permutationMatrix); // We set it to identity matrix, and in the end of the process get the permutation
        bool isPermutationIdentity = 1;
        double* L = (double*) malloc(numRow * numRow * sizeof(double));
        setIdentityMatrix(numRow, L);
        double* U = (double*) malloc(numRow * numColumn * sizeof(double));

        double* Q = (double*) malloc(numRow * numRow * sizeof(double));
        double* R = (double*) malloc(numRow * numColumn * sizeof(double));
        double* economy_Q = (double*) malloc(numRow * numColumn * sizeof(double));
        double* economy_R = (double*) malloc(numColumn * numColumn * sizeof(double));
        QRDecomposition(numRow, numColumn, inputMatrix, Q, R, economy_Q, economy_R);

        printf("Your \x1b[96mmatrix\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, inputMatrix);


        gaussJordanAndFindInverse(numRow, numColumn, echelonFormMatrix, inverseMatrix, permutationMatrix, &isPermutationIdentity, NULL, NULL);

        double* inputMatrixCopy = (double*) malloc(numRow * numColumn * sizeof(double));
        matrixCopy(numRow, numColumn, inputMatrix, inputMatrixCopy);
        LUDecomposition(numRow, numColumn, inputMatrixCopy, permutationMatrix, L, U);
        free(inputMatrixCopy);

        int rank = negativeZerosAndFindRank(numRow, numColumn, echelonFormMatrix);
        negativeZerosAndFindRank(numRow, numColumn, inverseMatrix);

        printf("The matrix in \x1b[92mreduced row echelon form\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, echelonFormMatrix);

        printf("The \x1b[93mrank\x1b[0m is: %d\n\n", rank);

        double deter;
        if (isSquare)
        {
            deter = deterCalc(numColumn, inputMatrix);

            if (!areEqual(deter, 0.0)) // If deter != 0, then the matrix has an inverse
            {
                printf("The \x1b[94minverse matrix\x1b[0m is:\n");
                displayMatrix(numRow, numColumn, inverseMatrix);
            }
            else printf("This matrix is \x1b[94msingular\x1b[0m\n\n");

            printf("The \x1b[91mdeterminant\x1b[0m is: %lf\n\n", deter);

            printf("The \x1b[95madjoint matrix\x1b[0m is:\n\n");
            findAndPrintAdjoint(numColumn, inputMatrix);
       }
        else
            printf("Non-square matrices don't have an \x1b[94minverse matrix\x1b[0m, "
            "nor a \x1b[91mdeterminant\x1b[0m, nor an \x1b[95madjoint\x1b[0m\n\n\n");


        printf("The \x1b[32mLU\x1b[0m decomposition is:\n");
        printf("L:\n");
        displayMatrix(numRow, numRow, L);
        printf("U:\n");
        displayMatrix(numRow, numColumn, U);
        if (isSquare && !areEqual(deter, 0.0) && isPermutationIdentity)
            printf("This matrix is regular (square, invertible, and has an \x1b[32mLU\x1b[0m decomposition without a permutation)\n\n");
        else
        {
            printf("With a permutation:\n\n");
            displayMatrix(numRow, numRow, permutationMatrix);
        }

        if (isSquare && !areEqual(deter, 0.0))
        {
            printf("The \x1b[33mLDV\x1b[0m decomposition is:\n");
            printf("L:\n");
            displayMatrix(numRow, numRow, L);
            findDVFromUAndPrint(numRow, U);
        }
        else printf("Non-square or singular matrices don't have an \x1b[33mLDV\x1b[0m decomposition\n\n");


        printf("The \x1b[31mQR\x1b[0m decomposition is:\n");
        printf("Q:\n");
        displayMatrix(numRow, numRow, Q);
        printf("R:\n");
        displayMatrix(numRow, numColumn, R);

        if (numColumn < numRow)
        {
            printf("This is a \"tall\" matrix, so it has an \"economy\" \x1b[31mQR\x1b[0m decomposition as well:\n");
            printf("Economy Q:\n");
            displayMatrix(numRow, numColumn, economy_Q);
            printf("Economy R:\n");
            displayMatrix(numColumn, numColumn, economy_R);
        }


        printf("The \x1b[36mFrobenius norm\x1b[0m is: %0.10lf\n\n", frobeniusNorm(numRow, numColumn, inputMatrix));

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

        double* A = (double*) malloc(numRow * numColumn * sizeof(double));
        doubleVerify(numRow, numColumn, A, 1);

        printf("Enter the free coefficients (\x1b[34mb\x1b[0m):\n");
        double* b = (double*) malloc(numRow * sizeof(double));
        doubleVerify(numRow, 1, b, 0);

        double* A_star = insertColumn(numRow, numColumn, numColumn, A, b, 'e');
        printf("The scalars of your system are:\n");
        displayMatrix(numRow, numColumn + 1, A_star);

        double* A_copy = (double*) malloc(numRow * numColumn * sizeof(double)); // We want 'A' to remain the same, so we copy to find rank
        matrixCopy(numRow, numColumn, A, A_copy);
        gaussJordanAndFindInverse(numRow, numColumn, A_copy, NULL, NULL, NULL, NULL, NULL);
        int A_rank = negativeZerosAndFindRank(numRow, numColumn, A_copy);
        free(A_copy);

        gaussJordanAndFindInverse(numRow, numColumn + 1, A_star, NULL, NULL, NULL, NULL, NULL); // We don't care about A_star besides its rank
        int A_star_rank = negativeZerosAndFindRank(numRow, numColumn + 1, A_star);
        free(A_star);

        double* x = NULL;
        double* remainder = leastSquares(numRow, numColumn, A, &x, b);


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

        while (numColumn1 != numRow2) // Verify the dimensions agree
        {
            printf("Number of columns of \x1b[91mmatrix 1\x1b[0m must be equal to number of rows of \x1b[94mmatrix 2\x1b[0m. Please try again:\n");

            printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
            getRowAndColumnNum(&numRow1, &numColumn1);
            printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
            getRowAndColumnNum(&numRow2, &numColumn2);
        }

        printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
        double* matrix1 = (double*) malloc(numRow1 * numColumn1 * sizeof(double));
        doubleVerify(numRow1, numColumn1, matrix1, 1);

        printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
        double* matrix2 = (double*) malloc(numRow2 * numColumn2 * sizeof(double));
        doubleVerify(numRow2, numColumn2, matrix2, 1);

        double* product = matrixMultiplier(numRow1, numColumn1, numColumn2, matrix1, matrix2);
        negativeZerosAndFindRank(numRow1, numColumn2, product);
        printf("The \x1b[92mproduct\x1b[0m is:\n");
        displayMatrix(numRow1, numColumn2, product);

        free(matrix1);
        free(matrix2);
        free(product);
        break;
    }

    case 4: // Convolution calculator
    {
        int size1, size2;

        printf("For \x1b[91mvector 1\x1b[0m:\n\n");
        getVectorSize(&size1);

        printf("For \x1b[94mvector 2\x1b[0m:\n\n");
        getVectorSize(&size2);

        printf("For \x1b[91mvector 1\x1b[0m:\n");
        printf("Enter vector entries: ");
        double* vector1 = (double*) malloc(size1 * sizeof(double));
        doubleVerify(1, size1, vector1, 0);

        printf("For \x1b[94mvector 2\x1b[0m:\n");
        printf("Enter vector entries: ");
        double* vector2 = (double*) malloc(size2 * sizeof(double));
        doubleVerify(1, size2, vector2, 0);


        double* z = linearConvolution(size1, vector1, size2, vector2);
        printf("The \x1b[92mlinear\x1b[0m convolution is:\n");
        displaySolution(size1 + size2 - 1, z);
        free(z);

        z = circularConvolution(size1, vector1, size2, vector2);
        printf("The \x1b[32mcircular\x1b[0m convolution is:\n");
        displaySolution(max(size1, size2), z);
        free(z);

        free(vector1);
        free(vector2);
        break;
    }
    }

    system("pause");

    return 0;
}