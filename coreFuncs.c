#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coreFuncsDec.h"

// Getting rid of those annoying negative zeros, and returning rank
int negativeZerosAndFindRank(int numColumn, int numRow, double *matrix)
{
   int flag = 0, nonZeroRows = 0;
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
void getMinor(int n, int rowValue, int columnValue, double *originMatrix, double *dupMatrix)
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

/* Calculating determinant recursively, input array is 2d (a matrix), 'n' is number of columns
Yes, I know we can multiply the numbers on the diagonal of an echelon form matrix, but isn't recursion fun? */
double deterCalc(int n, double *originMatrix)
{
   if (n == 1) return *(originMatrix); // Base case, Det(a) = a

   double *dupMatrix = (double*) malloc( (n - 1) * (n - 1) * sizeof(double)); // This array represents a minor
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

void findAndPrintAdjoint(int n, double *matrix)
{
   if (n == 1)
   {
      printf("%lf\n\n", 1.0);
      return;
   }

   double *dupArray = (double*) malloc( (n - 1) * (n - 1) * sizeof(double)); // This array stores the minor
   double *adjoint = (double*) malloc(n * n * sizeof(double));

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
}

/* This function takes the 'A' and 'b' in Ax=b, and stores the 'x', if there is a solution, in 'solution'
The calculation is done with Cramer's rule */
int CramersRule(int size, double *solution, double *scalarMatrix, double *bColumn)
{
   double matrixDeter = deterCalc(size, scalarMatrix);

   if (matrixDeter == 0) return 0; // No single solution

   for (int i = 0; i < size; i++)
   {
      switchColumns(size, size, i, scalarMatrix, bColumn); // Switching to calculate determinant ('size' represents the b vector in Ax=b)

      *(solution + i) = deterCalc(size, scalarMatrix) / matrixDeter; // The actual formula

      switchColumns(size, size, i, scalarMatrix, bColumn); // Switching back
   }

   return 1; // A single solution found
}

/* 'input' is passed to perform elimination on, and an identity matrix is passed as 'inverse'
Through elimination, as 'input' gets echelon form, 'inverse' becomes the inverse matrix (if 'input' is invertible)
If we don't want to find the inverse, we pass NULL
*Remember that 'matrix[i][j] = matrix + i*numColumn + j' */
void GaussJordanAndFindInverse(int numRow, int numColumn, double *inputMatrix, double *inverseMatrix)
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

         switchRows(numColumn, i, z, inputMatrix);
         if (inverseMatrix != NULL) switchRows(numColumn, i, z, inverseMatrix);

         pivot = *(inputMatrix + z * numColumn + j);
         for (k = 0; k < numColumn; k++) // Dividing the row by leading entry (pivot)
         {
            *(inputMatrix + z*numColumn + k) /= pivot;
            if (inverseMatrix != NULL) *(inverseMatrix + z*numColumn + k) /= pivot;
         }

         for (k = z + 1; k < numRow; k++) // Subtracting pivot's row ('z's row) from the rows below it
         {
            scalar = *(inputMatrix + k * numColumn + j);
            for (p = 0; p < numColumn; p++)
            {
               *(inputMatrix + k*numColumn + p) -= scalar * *(inputMatrix + z * numColumn + p);
               if (inverseMatrix != NULL) *(inverseMatrix + k*numColumn + p) -= scalar * *(inverseMatrix + z * numColumn + p);

               if (fabs(*(inputMatrix + k*numColumn + p)) < 1E-10) *(inputMatrix + k*numColumn + p) = 0; // Handling rounding errors
               if (inverseMatrix != NULL) if (fabs(*(inverseMatrix + k*numColumn + p)) < 1E-10) *(inverseMatrix + k*numColumn + p) = 0;
            }
         }
         z++;
         break; // We've zeroed all elements below pivot, and we're done with this column
      }
   }


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
                  if (inverseMatrix != NULL) *(inverseMatrix + k*numColumn + p) -= scalar * *(inverseMatrix + i * numColumn + p);
               }
            }
         break;
      }
   }
}

/* This function multiplies 'matrix1' and 'matrix2' (in this order)
Note that 'numColumn1 == numRow2', so no need for another variable */
double* matrixMultiplier(int numRow1, int numColumn1, int numColumn2, double *matrix1, double *matrix2)
{
   double *product = (double*) malloc(numRow1 * numColumn2 * sizeof(double));

   for (int i = 0; i < numRow1; i++)
   {
      for (int j = 0; j < numColumn2; j++)
      {
         *(product + i*numColumn2 + j) = vectorMultiplier(numColumn1, numColumn2, matrix1 + numColumn1*i, matrix2 + j);
      }
   }

   return product;
}
