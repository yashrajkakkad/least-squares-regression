/*
    Written as per ANSI C Standards. Can be verified by '-std=c89'
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Macros */
#define _Bool int
#define true 1
#define false 0

/* Uncomment the below line for input from files */
/* #define FILE_IN 1 */

/* Preprocessor directives for debugging */
/* #define DEBUG_MAIN 1 */
/* #define DEBUG_RREF 1 */
/* #define DEBUG_PS 1 */
/* #define DEBUG_NS 1 */

/* Encapsulating the matrix in a structure */
typedef struct
{
    int rows, cols;
    long double **data;
} matrix;

/* Function declarations */
int min(int a, int b);
matrix *newMatrix(int rows, int cols);
void readMatrix(matrix *m, FILE *fp);
void printMatrix(matrix *m);
matrix *transpose(matrix *m);
matrix *mulMatrix(matrix *m1, matrix *m2);
void switchRows(matrix *m, int row1, int row2);
void rref(matrix *A, matrix *b);
matrix *particularSolution(matrix *R, matrix *b);
_Bool isZero(long double a, long double threshold);

int min(int a, int b)
{
    return (a < b) ? a : b;
}

/* Return a pointer to a matrix struct after allocating appropriate space for data */
matrix *newMatrix(int rows, int cols)
{
    int i;
    matrix *m = (matrix *)malloc(sizeof(matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = (long double **)calloc(rows, sizeof(long double *));
    for (i = 0; i < rows; i++)
    {
        m->data[i] = (long double *)calloc(cols, sizeof(long double));
    }
    return m;
}

/* Read data to a matrix from a file pointer (stdin for console) */
void readMatrix(matrix *m, FILE *fp)
{
    int row, col;
    for (row = 0; row < m->rows; row++)
    {
        for (col = 0; col < m->cols; col++)
        {
            fscanf(fp, "%Lf", &(m->data[row][col]));
        }
    }
}

/* Print content of a matrix */
void printMatrix(matrix *m)
{
    int row, col;
    for (row = 0; row < m->rows; row++)
    {
        for (col = 0; col < m->cols; col++)
        {
            printf("%d %d %Lf ", row, col, m->data[row][col]);
        }
        /*        printf("%lf", m->data[row][m->cols - 1]); */
        printf("\n");
    }
}

matrix *transpose(matrix *m)
{
    int row, col;
    matrix *mt = newMatrix(m->cols, m->rows);
    for (row = 0; row < mt->rows; row++)
    {
        for (col = 0; col < mt->cols; col++)
        {
            mt->data[row][col] = m->data[col][row];
        }
    }
    return mt;
}

matrix *mulMatrix(matrix *m1, matrix *m2)
{
    int i, j, k;
    matrix *m = newMatrix(m1->rows, m2->cols);
    long double sum = 0;
    if (m1->cols != m2->rows)
    {
        return NULL;
    }
    for (i = 0; i < m->rows; i++)
    {
        for (j = 0; j < m->cols; j++)
        {
            m->data[i][j] = 0;
        }
    }
    for (i = 0; i < m1->rows; i++)
    {
        for (j = 0; j < m2->cols; j++)
        {
            for (k = 0; k < m1->cols; k++)
            {
                sum += (m1->data[i][k] * m2->data[k][j]);
            }
            m->data[i][j] = sum;
            sum = 0;
        }
    }
    return m;
}

/* Switch two rows - row1 and row 2 */
void switchRows(matrix *m, int row1, int row2)
{
    double *arr = calloc(m->cols, sizeof(long double));
    int i;
    if (row1 >= m->rows || row2 >= m->rows)
    {
        printf("Invalid row shift!\n");
        return;
    }
    for (i = 0; i < m->cols; i++)
    {
        arr[i] = m->data[row1][i];
    }
    for (i = 0; i < m->cols; i++)
    {
        m->data[row1][i] = m->data[row2][i];
    }
    for (i = 0; i < m->cols; i++)
    {
        m->data[row2][i] = arr[i];
    }
    free(arr);
}

/* Convert a matrix to rref form */
void rref(matrix *A, matrix *b)
{
    int pivotRow = 0;
    int pivotCol = 0;
    int j;
    int switchRow;
    long double r;
    int k;
    long double quotient;
    /* Converting to upper triangular */
    while (pivotRow < A->rows && pivotCol < A->cols)
    {
        switchRow = -1;
        if (isZero(A->data[pivotRow][pivotCol], 1e-5))
        {
#ifdef DEBUG_RREF
            printf("Zero at: %d %d\n", pivotRow, pivotCol);
#endif
            for (j = pivotRow + 1; j < A->rows; j++)
            {
                if (!isZero(A->data[j][pivotCol], 1e-5))
                {
                    switchRow = j;
#ifdef DEBUG_RREF
                    printf("Switch %d %d\n", pivotRow, j);
#endif
                    switchRows(A, pivotRow, switchRow);
                    switchRows(b, pivotRow, switchRow);
                    break;
                }
            }
            if (switchRow == -1)
            {
#ifdef DEBUG_RREF
                printf("No switch for row: %d\n", pivotRow);
#endif
                pivotCol++;
                continue;
            }
        }
        for (j = pivotRow + 1; j < A->rows; j++)
        {
            if (!isZero(A->data[j][pivotCol], 1e-5))
            {
                r = A->data[j][pivotCol] / A->data[pivotRow][pivotCol];
#ifdef DEBUG_RREF
                printf("%d %d %d %lf %lf r : %lf\n\n", pivotRow, pivotCol, j, A->data[j][pivotCol], A->data[pivotRow][pivotCol], r);
#endif
                for (k = pivotCol; k < A->cols; k++)
                {
                    A->data[j][k] -= r * A->data[pivotRow][k];
                }
                b->data[j][0] -= r * b->data[pivotRow][0];
            }
#ifdef DEBUG_RREF
            printMatrix(A);
            printf("\n");
#endif
        }
        pivotRow++;
        pivotCol++;
    }
#ifdef DEBUG_RREF
    printf("A converted to upper triangular: ");
    printMatrix(A);
    printMatrix(b);
#endif

    /*  Converting to rref from upper triangular
        Select a pivot */
    pivotRow = 0;
    pivotCol = 0;
    while (pivotRow < A->rows && pivotCol < A->cols)
    {
        if (isZero(A->data[pivotRow][pivotCol], 1e-5))
        {
            pivotCol++;
            continue;
        }
        else
        {
            for (j = pivotRow - 1; j >= 0; j--)
            {
                r = A->data[j][pivotCol] / A->data[pivotRow][pivotCol];
#ifdef DEBUG_RREF
                printMatrix(A);
                printf("\n");
                printf("%d %d %d %lf %lf r : %lf\n\n", pivotRow, pivotCol, j, A->data[j][pivotCol], A->data[pivotRow][pivotCol], r);
                // printf("%d %d %lf %lf r : %lf\n", , j, A->data[j][i], A->data[i][i], r);
#endif
                /* Row operations */
                for (k = pivotCol; k < A->cols; k++)
                {
                    A->data[j][k] -= r * A->data[pivotRow][k];
                }
                b->data[j][0] -= r * b->data[pivotRow][0];
            }
        }
        pivotRow++;
        pivotCol++;
    }
    /* Make all pivot elements '1' */
    pivotRow = 0;
    pivotCol = 0;
    while (pivotRow < A->rows && pivotCol < A->cols)
    {
        if (isZero(A->data[pivotRow][pivotCol], 1e-5))
        {
            pivotCol++;
            continue;
        }
        quotient = A->data[pivotRow][pivotCol];
        for (j = pivotCol; j < A->cols; j++)
        {
            if (A->data[pivotRow][j])
                A->data[pivotRow][j] /= quotient;
        }
        if (b->data[pivotRow][0])
            b->data[pivotRow][0] /= quotient;
        pivotRow++;
        pivotCol++;
    }
}

/* Check if there are any inconsistencies in Rx = b */
_Bool solutionExists(matrix *R, matrix *b)
{
    int i, j;
    _Bool zeroFlag;
    for (i = 0; i < R->rows; i++)
    {
        zeroFlag = true;
        for (j = 0; j < R->cols; j++)
        {
            if (!isZero(R->data[i][j], 1e-3))
            {
                zeroFlag = false;
                break;
            }
        }
        if (zeroFlag)
        {
            if (!isZero(b->data[i][0], 1e-3))
            {
                return false;
            }
        }
    }
    return true;
}

/* Deduce particular solution from reduced row echelon matrix with pivot elements 1 */
matrix *particularSolution(matrix *R, matrix *b)
{
    matrix *xp = newMatrix(R->cols, 1);
    int pivotRow = 0, pivotCol = 0;
    while (pivotRow < R->rows && pivotCol < R->cols)
    {
        if (!R->data[pivotRow][pivotCol])
        {
            pivotCol++;
            continue;
        }
        xp->data[pivotCol][0] = b->data[pivotRow][0] / R->data[pivotRow][pivotCol];
#ifdef DEBUG_PS
        printf("%d %d\n", pivotRow, pivotCol);
#endif
        pivotRow++;
        pivotCol++;
    }
#ifdef DEBUG_PS
    printf("Particular solution calculated\n");
#endif
    return xp;
}

/* Check if a number is approximately zero, upto a threshold. Used to deal with precision issues */
_Bool isZero(long double a, long double threshold)
{
    if (a >= -threshold && a <= threshold)
    {
        return true;
    }
    return false;
}
