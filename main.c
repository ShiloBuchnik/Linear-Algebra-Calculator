#include <stdio.h>
#include <stdlib.h>
#include "auxFuncsDec.h"
#include "coreFuncsDec.h"

int main()
{
    enableColor();
    int mode = printOpeningAndGetNum();

    //printText(mode);

    switch(mode)
    {
    case 1: // Matrix calculator
    {
        int numRow, numColumn, rank;
        double deter;

        getRowAndColumnNum(&numRow, &numColumn);
        int isSquare = (numRow == numColumn);

        double* inputMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
        doubleVerify(numRow, numColumn, inputMatrix, "row");

        double* echelonFormMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
        arrayCopy(numRow, numColumn, inputMatrix, echelonFormMatrix);


        double* inverseMatrix = NULL;
        if (isSquare)
        {
            inverseMatrix = (double*) malloc(numRow * numColumn * sizeof(double));
            setIdentityMatrix(numRow, inverseMatrix); // We set it to identity matrix, and in the end of the process get the inverse
        }

        double* permutationMatrix = (double*) malloc(numRow * numRow * sizeof(double));
        setIdentityMatrix(numRow, permutationMatrix); // We set it to identity matrix, and in the end of the process get the permutation
        int isPermutationIdentity = 1;
        double* L = (double*) malloc(numRow * numRow * sizeof(double));
        setIdentityMatrix(numRow, L);
        double* U = (double*) malloc(numRow * numColumn * sizeof(double));

        printf("Your \x1b[96mmatrix\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, inputMatrix);

        GaussJordanAndFindInverse(numRow, numColumn, echelonFormMatrix, inverseMatrix, permutationMatrix, &isPermutationIdentity, NULL, NULL);

        double* inputMatrixCopy = (double*) malloc(numRow * numColumn * sizeof(double));
        arrayCopy(numRow, numColumn, inputMatrix, inputMatrixCopy);
        LUDecomposition(numRow, numColumn, inputMatrixCopy, permutationMatrix, L, U);
        free(inputMatrixCopy);


        rank = negativeZerosAndFindRank(numRow, numColumn, echelonFormMatrix);
        negativeZerosAndFindRank(numRow, numColumn, inverseMatrix);

        printf("The matrix in \x1b[92mreduced row echelon form\x1b[0m is:\n");
        displayMatrix(numRow, numColumn, echelonFormMatrix);

        printf("The \x1b[93mrank\x1b[0m is: %d\n\n", rank);

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

            printf("The \x1b[95madjoint matrix\x1b[0m is:\n");
            findAndPrintAdjoint(numColumn, inputMatrix);
       }
        else
            printf("Non-square matrices don't have an \x1b[94minverse matrix\x1b[0m, "
            "nor a \x1b[91mdeterminant\x1b[0m, nor an \x1b[95madjoint\x1b[0m\n\n");


        printf("The LU decomposition is:\n");
        printf("L:\n");
        displayMatrix(numRow, numRow, L);
        printf("U:\n");
        displayMatrix(numRow, numColumn, U);
        printf("With a permutation:\n");
        displayMatrix(numRow, numRow, permutationMatrix);
        if (isSquare && !areEqual(deter, 0.0) && isPermutationIdentity)
            printf("This matrix is regular (square, invertible, and has an LU decomposition without a permutation)\n");

        if (isSquare && !areEqual(deter, 0.0))
        {
            printf("The LDV decomposition is:\n");
            printf("L:\n");
            displayMatrix(numRow, numRow, L);
            findDVFromUAndPrint(numRow, U);
        }
        else printf("Non-square or singular matrices don't have an LDV decomposition\n");


        free(inputMatrix);
        free(echelonFormMatrix);
        free(inverseMatrix);
        free(permutationMatrix);
        free(L);
        free(U);
        break;
    }

    case 2: // Set of equations calculator
    {
        int size;

        printf("Enter size of system: ");
        positiveIntVerify(&size, NULL, 0);
        while (getchar() != '\n');
        printf("\n");

        double* solution = (double*) malloc(size * sizeof(double));
        double* inputMatrix = (double*) malloc(size * (size + 1) * sizeof(double)); // The extra column in 'inputMatrix' is for the b vector in Ax=b
        doubleVerify(size, size + 1, inputMatrix, "equation");

        printf("The \x1b[96mscalars\x1b[0m of your system are:\n");
        displayMatrix(size, size + 1, inputMatrix);

        double* bColumn = (double *) malloc(size * sizeof(double)); // stores the b vector from Ax=b
        double* scalarMatrix = separateColumn(size, size + 1, size, inputMatrix, bColumn); // Stores the A matrix from Ax=b
        free(inputMatrix);

        if (CramersRule(size, solution, scalarMatrix, bColumn))
        {
            negativeZerosAndFindRank(1, size, solution);
            printf("The \x1b[92msolution vector\x1b[0m is: ");
            displaySolution(size, solution);
        }
        else printf("The system doesn't have a single solution (might have infinite or none)\n\n");

        free(scalarMatrix);
        free(solution);
        free(bColumn);
        break;
    }

    case 3: // Matrix multiplier
    {
        int numRow1, numColumn1, numRow2, numColumn2;

        printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
        getRowAndColumnNum(&numRow1, &numColumn1);

        printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
        getRowAndColumnNum(&numRow2, &numColumn2);

        if (numColumn1 != numRow2)
        {
            printf("Number of columns of \x1b[91mmatrix 1\x1b[0m must be equal to number of rows of \x1b[94mmatrix 2\x1b[0m. Terminating . . .\n\n");
            system("pause");
            exit(0);
        }

        printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
        double* matrix1 = (double*) malloc(numRow1 * numColumn1 * sizeof(double));
        doubleVerify(numRow1, numColumn1, matrix1, "row");

        printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
        double* matrix2 = (double*) malloc(numRow2 * numColumn2 * sizeof(double));
        doubleVerify(numRow2, numColumn2, matrix2, "row");

        double* product = matrixMultiplier(numRow1, numColumn1, numColumn2, matrix1, matrix2);
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