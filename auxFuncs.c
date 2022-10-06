#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

// Sending out error message and cleaning buffer
void invalidInput(void)
{
   printf("Invalid input, please try again: ");
   while (getchar() != '\n');
}

// Using Windows' API to enable color in program
int enableColor(void)
{
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE); // We need to get the windows handle of our stdout
    if (hOut == INVALID_HANDLE_VALUE) return -1; // If something went wrong - return -1

    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode)) return -1; // checking the output mode, idk, anyway if something went wrong - return -1

    dwMode |= 0x0004; // ENABLE_VIRTUAL_TERMINAL_PROCESSING macro, this is what we tell to setConsoleMode() to enable ansi escape sequences
    if (!SetConsoleMode(hOut, dwMode)) return -1;
    // setting the console mode to enable ansi escape sequences, if something went wrong - return -1

    return 1; // Return 1 on success
}

void swap(double *num1, double *num2)
{
   double temp = *num1;
   *num1 = *num2;
   *num2 = temp;
}

/* This function takes an input and verify that it's a *positive* int; and if it's in the array 'arr'.
if we only want to check whether the input is a positive int, we can pass NULL as 'arr', and an arbitrary int as 'arrLen' */
void positiveIntVerify(int *input, int *arr, int arrLen)
{
    while (1)
    {
        while (scanf("%d", input) != 1 || *input < 1) invalidInput(); // We want it to be an int and a positive one

        if (arr == NULL) return;
        else // checking if input is in the array
        {
            for (int i = 0; i < arrLen; i++)
            {
                if (*input == arr[i]) return;
            }
            invalidInput();
        }
    }
}

/* For getting matrix elements & verifying them. Verification part is a tad spaghetti, I just wanted to verify perfectly
parameter 'str' can be 'row' or 'equation', based on what mode the user chose */
void doubleVerify(int numColumn, int numRow, double *matrix, char *str)
{
   int i, j, symbol, flag = 0;
   double temp1, temp2;

   for (i = 0; i < numRow; i++)
   {
      printf("Enter %s #%d\n", str, i+1);

      for (j = 0; j < numColumn; j++)
      {
         while (flag == 1 || scanf("%lf", &temp1) != 1)
         {
            invalidInput();
            flag = 0;
         }

         symbol = getchar();

         if (symbol == '/' || symbol == '\\')
         {
            if (scanf("%lf", &temp2) != 1)
            {
               j--;
               flag = 1;
               continue;
            }
            *(matrix + i*numColumn + j) = temp1 / temp2;
         }

         else if (symbol == ' ' || symbol == '\n') *(matrix + i*numColumn + j) = temp1;

         else
         {
            j--;
            flag = 1;
            continue;
         }
      }
      printf("\n");
   }
}

// This function gets number of rows and number of columns from user
void getRowAndColumnNum(int *numRow, int *numColumn)
{
   printf("Enter row size: "); // Number of arrays
   positiveIntVerify (numRow, NULL, 0);
   while (getchar() != '\n');
   printf("Enter column size: "); // Number of elements in array
   positiveIntVerify (numColumn, NULL, 0);
   while (getchar() != '\n');
   printf("\n");
}

void displayMatrix(int numColumn, int numRow, double *matrix)
{
   for (int i = 0; i < numRow; i++)
   {
      for (int j = 0; j < numColumn; j++)
      {
         printf("%lf  ", *(matrix + i*numColumn + j));
      }
      printf("\n");
   }

   printf("\n");
}

// 'displaySolution' displays vectors in a slightly different way than 'displayMatrix'
void displaySolution(int size, double *solution) // For displaying the solution
{
   int i;

   printf("(");
   for (i = 0; i < size - 1; i++) printf("%lf, ", solution[i]);
   printf("%lf)\n\n", solution[i]); // I don't want a comma in the end :>
}

// Setting given matrix to identity matrix
void setIdentityMatrix(int numRow, int numColumn, double *matrix)
{
   for (int i = 0; i < numRow; i++)
   {
      for (int j = 0; j < numColumn; j++)
      {
         if (j == i) *(matrix + i*numColumn + j) = 1;
         else *(matrix + i*numColumn + j) = 0;
      }
   }
}

void arrayCopy(int numRow, int numColumn, double *srcMatrix, double *destMatrix)
{
   for (int i = 0; i < numRow; i++)
   {
      for (int j = 0; j < numColumn; j++)
      {
         *(destMatrix + i*(numColumn) + j) = *(srcMatrix + i*(numColumn) + j);
      }
   }
}

/* The function takes a matrix, a column vector and a number 'column'.
it swaps between the column vector and the column in the matrix whose number is stored in 'column' */
void switchColumns(int numRow, int numColumn, int column, double *matrix, double *arrColumn)
{
   if (column < 0 || column >= numColumn) return;

   for (int i = 0; i < numRow; i++) swap(matrix + i*numColumn + column, arrColumn + i);
}

void switchRows(int numColumn, int row1, int row2, double *matrix)
{
   for (int j = 0; j < numColumn; j++) swap(matrix + row1*numColumn + j, matrix + row2*numColumn + j);
}

/* This function takes a matrix and separate a column (listed as 'column') from it.
It returns a pointer to a new matrix without that column, and it stores that column in 'arrColumn' */
double* separateColumn(int numRow, int numColumn, int column, double *matrix, double *arrColumn)
{
   if (column < 0 || column >= numColumn) return NULL;

   int i, j, x = 0;
   double *scalarMatrix = (double*) malloc(numRow * (numColumn - 1) * sizeof(double));

   for (i = 0; i < numRow; i++)
   {
      for (j = 0; j < numColumn; j++)
      {
         if (j == column) *(arrColumn + i) = *(matrix + i*numColumn + j);
         else
         {
            *(scalarMatrix + x) = *(matrix + i*numColumn + j);
            x++;
         }
      }
   }

   return scalarMatrix;
}

double vectorMultiplier(int column1, int column2, double *vector1, double *vector2)
{
   double sum = 0;
   for (int i = 0; i < column1; i++) sum += *(vector1 + i) * *(vector2 + i*column2);

   return sum;
}

// Print introduction and get mode from user
int printOpeningAndGetNum(void)
{
   printf("Hello! This program offers: \n");
   printf("\x1b[92m1)\x1b[0m Various matrix-related calculations\n");
   printf("\x1b[91m2)\x1b[0m Square system of equations (number of variables = number of equations) calculator\n");
   printf("\x1b[94m3)\x1b[0m Matrix multiplication calculator\n\n");
   printf("Please input your choice[\x1b[92m1\x1b[0m/\x1b[91m2\x1b[0m/\x1b[94m3\x1b[0m]: ");

   int num, *arr = (int*) malloc(3 * sizeof(int));
   arr[0] = 1; arr[1] = 2; arr[2] = 3;

   positiveIntVerify(&num, arr, 3);
   free(arr);

   system("cls");

   return num;
}

// Printing instructions for user's choice of mode
// Did you know? Separating strings in 'printf' with whitespace - concatenate them
void printText(int num)
{
   if (num == 1)
   {
      printf("In this program you input a \x1b[96mmatrix\x1b[0m\n"
             "The output is the \x1b[92mreduced row echelon form\x1b[0m,\n"
             "as well as the \x1b[93mrank\x1b[0m, the \x1b[94minverse matrix\x1b[0m, the \x1b[91mdeterminant\x1b[0m and the \x1b[95madjoint matrix\x1b[0m\n"
             "Please input each row as a string of numbers, with spaces between them. For example: 0.3 1/2 5 8 23");
   }
   else if (num == 2)
   {
      printf("In this program you input, for each equation, the \x1b[96mscalars\x1b[0m of its variables by order\n"
             "So, for example, for equation #1, '2x+3y=1', the input should be '2 3 1'\n"
             "The output is \x1b[92msolution vector\x1b[0m for the set of equations");
   }
   else
   {
      printf("In this program you input two matrices: \x1b[91mmatrix 1\x1b[0m and \x1b[94mmatrix 2\x1b[0m\n"
             "The output is the \x1b[92mproduct\x1b[0m: \x1b[91mmatrix 1\x1b[0m * \x1b[94mmatrix 2\x1b[0m");
   }

    printf("\n\nThings to take heed of:\n"
            "1) If you type a non-number character in a row,\n"
            "the program would require you to retype every entry \x1b[4mfrom said non-number\x1b[0m to the end of the row\n"
            "2) If you wish to input a fraction, \x1b[4mdo not\x1b[0m include spaces before and after the fraction line\n"
            "This is valid: a/b\n"
            "This is not valid: a / b; a /b; a/ b\n"
            "3) This calculator is accurate up to 10 digits after the decimal\n"
            "4) Enjoy ^~^\n\n\n");
}
