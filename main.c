#include <stdio.h>
#include <stdlib.h>
#include "mainDec.h"

// Tip: *(matrix + i*numColumn + j) == matrix[i][j]

int main()
{
   enableColor();
   int mode = printOpeningAndGetNum();

   printText(mode);

   switch(mode) {
   case 1: // Matrix calculator
   {
       int numRow, numColumn, rank;
       double deter;

       getRowAndColumnNum(&numRow, &numColumn);

       double *inputMatrix = (double *) malloc(numRow * numColumn * sizeof(double));
       doubleVerify(numColumn, numRow, inputMatrix, "row");

       double *echelonFormMatrix = (double *) malloc(numRow * numColumn * sizeof(double));
       arrayCopy(numRow, numColumn, inputMatrix, echelonFormMatrix);

       double *inverseMatrix = (double *) malloc(numRow * numColumn * sizeof(double));
       setIdentityMatrix(numRow, numColumn, inverseMatrix); // We set it to identity matrix, and in the end of the process get the inverse

       printf("Your \x1b[96mmatrix\x1b[0m is:\n");
       displayMatrix(numColumn, numRow, inputMatrix);

       GaussJordanAndFindInverse(numRow, numColumn, echelonFormMatrix, inverseMatrix);

       rank = negativeZerosAndFindRank(numColumn, numRow, echelonFormMatrix);
       negativeZerosAndFindRank(numColumn, numRow, inverseMatrix);

       printf("The matrix in \x1b[92mreduced row echelon form\x1b[0m is:\n");
       displayMatrix(numColumn, numRow, echelonFormMatrix);

       printf("The \x1b[93mrank\x1b[0m is: %d\n\n", rank);

       if (numRow == numColumn) {
           deter = deterCalc(numColumn, inputMatrix);

           if (deter != 0) // If deter != 0, then the matrix has an inverse
           {
               printf("The \x1b[94minverse matrix\x1b[0m is:\n");
               displayMatrix(numColumn, numRow, inverseMatrix);
           } else printf("This matrix is \x1b[94msingular\x1b[0m\n\n");

           printf("The \x1b[91mdeterminant\x1b[0m is: %lf\n\n", deter);

           printf("The \x1b[95madjoint matrix\x1b[0m is:\n");
           findAndPrintAdjoint(numColumn, inputMatrix);
       } else
           printf("Non-square matrices don't have an \x1b[94minverse matrix\x1b[0m, nor a \x1b[91mdeterminant\x1b[0m, nor an \x1b[95madjoint\x1b[0m\n\n");

       free(inputMatrix);
       free(echelonFormMatrix);
       free(inverseMatrix);
       break;
   }

   case 2: // Set of equations calculator
   {
       int size;

       printf("Enter size of system: ");
       positiveIntVerify(&size, NULL, 0);
       while (getchar() != '\n');
       printf("\n");

       double *solution = (double *) malloc(size * sizeof(double));
       double *inputMatrix = (double *) malloc(
               size * (size + 1) * sizeof(double)); // The extra column in 'inputMatrix' is for the b vector in Ax=b
       doubleVerify(size + 1, size, inputMatrix, "equation");

       printf("The \x1b[96mscalars\x1b[0m of your system are:\n");
       displayMatrix(size + 1, size, inputMatrix);

       double *bColumn = (double *) malloc(size * sizeof(double)); // stores the b vector from Ax=b
       double *scalarMatrix = separateColumn(size, size + 1, size, inputMatrix, bColumn); // Stores the A matrix from Ax=b
       free(inputMatrix);

       if (CramersRule(size, solution, scalarMatrix, bColumn)) {
           negativeZerosAndFindRank(size, 1, solution);
           printf("The \x1b[92msolution vector\x1b[0m is: ");
           displaySolution(size, solution);
       } else printf("The system doesn't have a single solution (might have infinite or none)\n\n");

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

       if (numColumn1 != numRow2) {
           printf("Number of columns of \x1b[91mmatrix 1\x1b[0m must be equal to number of rows of \x1b[94mmatrix 2\x1b[0m. Terminating . . .\n\n");
           system("pause");
           exit(0);
       }

       printf("For \x1b[91mmatrix 1\x1b[0m:\n\n");
       double *matrix1 = (double *) malloc(numRow1 * numColumn1 * sizeof(double));
       doubleVerify(numColumn1, numRow1, matrix1, "row");

       printf("For \x1b[94mmatrix 2\x1b[0m:\n\n");
       double *matrix2 = (double *) malloc(numRow2 * numColumn2 * sizeof(double));
       doubleVerify(numColumn2, numRow2, matrix2, "row");

       double *product = matrixMultiplier(numRow1, numColumn1, numColumn2, matrix1, matrix2);
       negativeZerosAndFindRank(numColumn2, numRow1, product);
       printf("The \x1b[92mproduct\x1b[0m is:\n");
       displayMatrix(numColumn2, numRow1, product);

       free(matrix1);
       free(matrix2);
       free(product);
       break;
   }
   }

   system("pause");

   return 0;
}