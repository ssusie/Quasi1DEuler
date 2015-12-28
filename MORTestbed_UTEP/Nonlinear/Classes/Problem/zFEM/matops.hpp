#ifndef __MATOPS
#define __MATOPS

//Determinant of 3x3 matrix
double determinant3x3(double** A);

//Inverse of 3x3 matrix
void inverse3x3(double** A, double** Ai);

void eig3x3sym(double **A, double *lam);
#endif
