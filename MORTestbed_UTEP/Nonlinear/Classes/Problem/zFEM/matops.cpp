#include "matops.hpp"
#include <math.h>

#define   PI   4.0*atan(1.0)

double determinant3x3(double** A) {

    return   ( A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1]  );
}

void inverse3x3(double** A, double** Ai) {

    double Ji = 1.0/(determinant3x3(A));

    Ai[0][0] = Ji*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
    Ai[0][1] = Ji*(A[0][2]*A[2][1]-A[0][1]*A[2][2]);
    Ai[0][2] = Ji*(A[0][1]*A[1][2]-A[0][2]*A[1][1]);

    Ai[1][0] = Ji*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
    Ai[1][1] = Ji*(A[0][0]*A[2][2]-A[0][2]*A[2][0]);
    Ai[1][2] = Ji*(A[0][2]*A[1][0]-A[0][0]*A[1][2]);

    Ai[2][0] = Ji*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    Ai[2][1] = Ji*(A[0][1]*A[2][0]-A[0][0]*A[2][1]);
    Ai[2][2] = Ji*(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
}

void eig3x3sym(double **A, double *lam) {
    //This function computes the eigenvalues of a SYMMETRIC 3x3 matrix
    //Algorithm from Wikipedia page "Eigenvalue Algorithm"

    double p = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];

    if (p == 0){
      //A is diagonal
      lam[0]=A[0][0];
      lam[1]=A[1][1];
      lam[2]=A[2][2];
    } else {
      double q = (1.0/3.0)*(A[0][0]+A[1][1]+A[2][2]);
      p*=2;
      p += (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q) + (A[2][2]-q)*(A[2][2]-q);
      p*=(1.0/6.0);
      p = pow(p,0.5);

      A[0][0]-=q; A[1][1]-=q; A[2][2]-=q;
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          A[i][j]*=(1.0/q);

      double r = 0.5*determinant3x3(A);

      //In exact arithmetic for a symmetric matrix -1 <= r <= 1
      //but computation error can leave it slightly outside this range.
      double phi;
      if (r<=-1)
        phi = PI / (3.0);
      else if (r >= 1)
        phi = 0.0;
      else
        phi = acos(r) / 3.0;

      //the eigenvalues satisfiy eig1 >= eig2 >= eig3
      lam[0] = q + 2*p*cos(phi);
      lam[1] = q + 2*p*cos(phi + PI*(2.0/3.0));
      lam[2] = 3*q - lam[0] - lam[1];
    }
    return;
}
