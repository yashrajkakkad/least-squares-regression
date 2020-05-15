#include "matrix.h"
#include <stdio.h>

int main()
{
    int DEGREE = 2;
    int H[] = {160, 155, 143, 162, 170, 175, 167, 163, 156, 172, 180, 159, 169, 157, 161, 171, 181, 177, 140};
    int W[] = {68, 78, 56, 70, 72, 78, 68, 54, 43, 71, 65, 59, 63, 62, 90, 45, 62, 57, 80};
    int nItems = sizeof(H) / sizeof(int);
    int i;
    matrix *A = newMatrix(nItems, DEGREE + 1);
    matrix *b = newMatrix(nItems, 1);
    matrix *P, *q, *x;
    printf("Height Weight\n");
    for (i = 0; i < nItems; i++)
    {
        printf("%d\t%d\n", H[i], W[i]);
    }
    for (i = 0; i < nItems; i++)
    {
        A->data[i][0] = H[i] * H[i];
        A->data[i][1] = H[i];
        A->data[i][2] = 1;
        b->data[i][0] = W[i];
    }
    printf("We want the least squares solution for Ax = b.");
    printf("Matrix A: \n");
    printMatrix(A);
    printf("Vector b: \n");
    printMatrix(b);
    printf("We solve for At*A = At*b\n");
    P = mulMatrix(transpose(A), A);
    printf("At*A: \n");
    printMatrix(P);
    q = mulMatrix(transpose(A), b);
    printf("At*b:\n");
    printMatrix(q);
    rref(P, q);
    printf("Using Gaussian elimination, we obtain the solution: \n");
    x = particularSolution(P, q);
    printMatrix(x);
    printf("Therefore, in the second degree equation ah^2 + bh + c, a = %Lf, b = %Lf, c = %Lf \n", x->data[0][0], x->data[1][0], x->data[2][0]);
    printf("Predicted weight for h = 185cm : %Lf \n", x->data[0][0] * 185 * 185 + x->data[1][0] * 185 + x->data[2][0]);
    return 0;
}