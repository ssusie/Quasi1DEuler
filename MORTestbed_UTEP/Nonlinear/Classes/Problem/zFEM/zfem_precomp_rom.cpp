#include "zfem_forces.hpp"
#include "zfem_elems.hpp"
#include "matops.hpp"
#include "array.hpp"

#include <cmath>
#include <assert.h>
#include <stdio.h>

void  PrecompRomQuant(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha, array<double> beta, array<double> gamma, array<double> omega) {

/*    int cnt=0;
    for (int r1=0; r1 < nY; ++r1)
       alpha(r1)=++cnt;

    cnt=0;
    for (int r1=0; r1 < nY; ++r1)
       for (int r2=0; r2 < nY; ++r2)
          beta(r1,r2)=++cnt;
    cnt=0;
    for (int r1=0; r1 < nY; ++r1)
       for (int r2=0; r2 < nY; ++r2)
          for (int r3=0; r3 < nY; ++r3)
             gamma(r1,r2,r3)=++cnt;
    cnt=0;
    for (int r1=0; r1 < nY; ++r1)
       for (int r2=0; r2 < nY; ++r2)
          for (int r3=0; r3 < nY; ++r3)
             for (int r4=0; r4 < nY; ++r4)
                omega(r1,r2,r3,r4)=++cnt;

    return; */

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];

    size_t dims[4];

    //Entries of ROB on current element
    dims[0]=msh->nsd; dims[1]=msh->nen; dims[2]=nY;
    array<double> ve(3,dims);

    //Pre-computed FEM matrix ("A"); see notes
    dims[0]=msh->nsd; dims[1]=msh->nsd; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> A(4,dims); 
    array<double> Ab(4,dims); 

    //Pre-computed FEM matrix ("B"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> Bb(2,dims); 

    //Pre-computed FEM matrix ("C"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nsd;
    array<double> C(4,dims); 
    array<double> C1b(4,dims); 
    array<double> C2b(4,dims); 

    //Pre-computed FEM matrix ("D"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> D(4,dims); 
    array<double> Db(4,dims); 

    //Set vector to contain ubar_iI*ubar_iJ
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> u2(2,dims);

    //Set vector to contain ve2
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=nY; dims[3]=nY;
    array<double> ve2(4,dims);
    array<double> dNdN(2,dims);

    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    double** dN = new double* [msh->nen];
    double** dNdX = new double* [msh->nen];
    for (int i=0; i < msh->nen; ++i) 
     { dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    double J;

    alpha.zeroout();
    beta.zeroout();
    gamma.zeroout();
    omega.zeroout();

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              for (int k=0; k < nY; ++k) {
                 if ((msh->LM)(e, i+msh->nsd*j) > -1)
                    ve(i,j,k)=V( (msh->LM)( e, i+msh->nsd*j ) , k );
                 else
                    ve(i,j,k)=0.0;
              }
           }
        }

        /*if (e == 0) {
           for (int k=0; k < nY; ++k){
              for (int i=0; i < msh->nsd; ++i){
                 for (int j=0; j < msh->nen; ++j){
                    std::cout << ve(i,j,k) << "  ";
                 }
               std::cout << std::endl;
               }
           std::cout << std::endl;
           std::cout << std::endl;
           }
        }*/

        //Pre-compute terms used in computation of alpha, beta, gamma, and omega
        u2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1){
           for (int n2=0; n2<msh->nen; ++n2){
              for (int i=0; i < msh->nsd; ++i){ 
                 if ((msh->LM)(e, i+msh->nsd*n1) > -1 && (msh->LM)(e, i+msh->nsd*n2) > -1)
                    u2(n1,n2)+=ub( (msh->LM)( e, i+msh->nsd*n1) )*ub( (msh->LM)( e, i+msh->nsd*n2) );
              }
           }
        }
 
        ve2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1)
           for (int n2=0; n2<msh->nen; ++n2)
              for (int r1=0; r1 < nY; ++r1)
                 for (int r2=0; r2 < nY; ++r2)
                    for (int s=0; s < msh->nsd; ++s)
                       ve2(n1,n2,r1,r2)+=ve(s,n1,r1)*ve(s,n2,r2);

        A.zeroout();
        C.zeroout();
        D.zeroout();

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 dNdX[i][j]=0.0;
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];

           //Compute dNdX*dNdX'
           dNdN.zeroout();
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int s=0; s < msh->nsd; ++s)
                    dNdN(n1,n2)+=dNdX[n1][s]*dNdX[n2][s];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int s1=0; s1 < msh->nsd; ++s1)
              for (int s2=0; s2 < msh->nsd; ++s2)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       A(s1,s2,n1,n2) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s1]*dNdX[n2][s2];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int s=0; s < msh->nsd; ++s)
                       C(n1,n2,n3,s) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s]*dNdN(n2,n3);
         
           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
               for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int n4=0; n4 < msh->nen; ++n4)
                       D(n1,n2,n3,n4) += (1.0/6.0)*J*(quad->W(k))*dNdN(n1,n2)*dNdN(n3,n4);
        }

       /* double sum = 0.0;
        for (int n1=0; n1 < msh->nen; ++n1)
            for (int n2=0; n2 < msh->nen; ++n2)
              for (int n3=0; n3 < msh->nen; ++n3)
                 for (int n4=0; n4 < msh->nen; ++n4)
                    sum += pow(pow(D(n1,n2,n3,n4) - D(n3,n4,n1,n2),2.0),0.5);
        if (sum > 1e-3)
           std::cout << "D Major Symmetry error of element " << e << " = " << sum << std::endl;*/

        /*for (int n2=0; n2 < msh->nen; ++n2){
           for (int n1=0; n1 < msh->nen; ++n1){
              for (int s1=0; s1 < msh->nsd; ++s1){
                 for (int s2=0; s2 < msh->nsd; ++s2){
                    std::cout << A(s1,s2,n1,n2) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }*/

        /*for (int s=0; s < msh->nsd; ++s){
           for (int n3=0; n3 < msh->nen; ++n3){
              for (int n1=0; n1 < msh->nen; ++n1){
                 for (int n2=0; n2 < msh->nen; ++n2){
                    std::cout << C(n1,n2,n3,s) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }*/

        /*for (int n4=0; n4 < msh->nen; ++n4){
           for (int n3=0; n3 < msh->nen; ++n3){
              for (int n1=0; n1 < msh->nen; ++n1){
                 for (int n2=0; n2 < msh->nen; ++n2){
                    std::cout << D(n1,n2,n3,n4) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }*/
         

        //std::cout << (phys[msh->mat(e)].mu) << std::endl;
        //Compute Abar from definition of A
        Ab.zeroout();
        for (int s1=0; s1 < msh->nsd; ++s1)
           for (int s2=0; s2 < msh->nsd; ++s2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    Ab(s1,s2,n1,n2) = (phys[msh->mat(e)].lam)*A(s1,s2,n1,n2) + (phys[msh->mat(e)].mu)*A(s1,s2,n2,n1);

        Bb.zeroout();
        //Compute Bbar from definition of A
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int s=0; s < msh->nsd; ++s)
                 Bb(n1,n2) += (phys[msh->mat(e)].mu)*A(s,s,n1,n2);

        //Compute C1bar from definition of C
        C1b.zeroout();
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int n3=0; n3 < msh->nen; ++n3)
                 for (int s=0; s < msh->nsd; ++s)
                    C1b(n1,n2,n3,s) = 0.5*(phys[msh->mat(e)].lam)*C(n1,n2,n3,s) + (phys[msh->mat(e)].mu)*C(n3,n2,n1,s);
        
        //Compute C2bar from definition of C
        C2b.zeroout();
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int n3=0; n3 < msh->nen; ++n3)
                 for (int s=0; s < msh->nsd; ++s)
                    C2b(n1,n2,n3,s) = (phys[msh->mat(e)].lam)*C(n1,n2,n3,s) + (phys[msh->mat(e)].mu)*(C(n2,n1,n3,s)+C(n3,n1,n2,s));

        //Compute Dbar from definition of D
        Db.zeroout();
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int n3=0; n3 < msh->nen; ++n3)
                 for (int n4=0; n4 < msh->nen; ++n4)
                    Db(n1,n2,n3,n4) = 0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4) + (phys[msh->mat(e)].mu)*D(n1,n3,n2,n4);
                    //Db(n1,n2,n3,n4) = 0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4) + (phys[msh->mat(e)].mu)*D(n2,n3,n1,n4);
                    //Db(n1,n2,n3,n4) = 0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4) + (phys[msh->mat(e)].mu)*D(n1,n3,n2,n4);

        /*sum = 0.0;
        for (int n1=0; n1 < msh->nen; ++n1)
            for (int n2=0; n2 < msh->nen; ++n2)
              for (int n3=0; n3 < msh->nen; ++n3)
                 for (int n4=0; n4 < msh->nen; ++n4)
                    sum += pow(pow(Db(n1,n2,n3,n4) - Db(n3,n4,n1,n2),2.0),0.5)/pow(pow(Db(n1,n2,n3,n4),2.0),0.5);
        if (sum > 1e-3)
           std::cout << "Db Major Symmetry error of element " << e << " = " << sum << std::endl;*/

        /*for (int n2=0; n2 < msh->nen; ++n2){
           for (int n1=0; n1 < msh->nen; ++n1){
              for (int s1=0; s1 < msh->nsd; ++s1){
                 for (int s2=0; s2 < msh->nsd; ++s2){
                    std::cout << A(s1,s2,n1,n2) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }*/

        /*for (int s=0; s < msh->nsd; ++s){
           for (int n3=0; n3 < msh->nen; ++n3){
              for (int n1=0; n1 < msh->nen; ++n1){
                 for (int n2=0; n2 < msh->nen; ++n2){
                    std::cout << C(n1,n2,n3,s) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }*/

        /*for (int n4=0; n4 < msh->nen; ++n4){
           for (int n3=0; n3 < msh->nen; ++n3){
              for (int n1=0; n1 < msh->nen; ++n1){
                 for (int n2=0; n2 < msh->nen; ++n2){
                    std::cout << Db(n1,n2,n3,n4) << "  ";
                  }
                  std::cout << std::endl;
                }
             std::cout << std::endl << std::endl;
             }
         std::cout << std::endl << std::endl << std::endl;
         }
       
         Bb.printFormatted();*/
        
        //Compute alpha
        for (int r1=0; r1 < nY; ++r1)
           alpha(r1)=0.0;
 
        //Compute beta
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    for (int s1=0; s1 < msh->nsd; ++s1)
                       for (int s2=0; s2 < msh->nsd; ++s2)
                          beta(r1,r2)+=Ab(s1,s2,n1,n2)*ve(s1,n1,r1)*ve(s2,n2,r2);

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    beta(r1,r2)+=Bb(n1,n2)*ve2(n1,n2,r1,r2);
             
        //Compute gamma  
        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int r3=0; r3 < nY; ++r3){
                 for (int n1=0; n1 < msh->nen; ++n1){
                    for (int n2=0; n2 < msh->nen; ++n2){
                       for (int n3=0; n3 < msh->nen; ++n3){
                          for (int s=0; s < msh->nsd; ++s){
                             //for (int t=0; t < msh->nsd; ++t){
                             //   gamma(r1,r2,r3)+=0.5*(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*ve(t,n2,r2)*ve(t,n3,r3)*ve(s,n1,r1);
                             //   gamma(r1,r2,r3)+=(phys[msh->mat(e)].lam)*C(n3,n1,n2,s)*ve(t,n1,r1)*ve(t,n2,r2)*ve(s,n3,r3);
                             //   gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n3,n1,n2,s)*ve(t,n2,r2)*ve(t,n3,r3)*ve(s,n1,r1);
                             //   gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n1,n2,n3,s)*ve(t,n1,r1)*ve(t,n2,r2)*ve(s,n3,r3);
                             //   gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n2,n1,n3,s)*ve(t,n1,r1)*ve(t,n2,r2)*ve(s,n3,r3);
                             //}
                             gamma(r1,r2,r3)+=0.5*(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*(ve2(n2,n3,r2,r3)*ve(s,n1,r1));
                             gamma(r1,r2,r3)+=(phys[msh->mat(e)].lam)*C(n3,n1,n2,s)*(ve2(n1,n2,r1,r2)*ve(s,n3,r3));
                             gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n3,n1,n2,s)*(ve2(n2,n3,r2,r3)*ve(s,n1,r1));
                             gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n1,n2,n3,s)*(ve2(n1,n2,r1,r2)*ve(s,n3,r3));
                             gamma(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n2,n1,n3,s)*(ve2(n1,n2,r1,r2)*ve(s,n3,r3));
                           }
                        }
                     }
                  }
               }
            }
         }

        /*for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma(r1,r2,r3)+=C1b(n1,n2,n3,s)*ve2(n2,n3,r2,r3)*ve(s,n1,r1);

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma(r1,r2,r3)+=C2b(n1,n2,n3,s)*ve2(n2,n3,r1,r2)*ve(s,n1,r3);*/

        if (e == 0){
           std::cout << "D = " << std::endl;
           for (int l=0; l<nY; ++l){
              for (int k=0; k<nY; ++k){
                for (int i=0; i<nY; ++i){
                   for (int j=0; j<nY; ++j){
                      std::cout << D(i,j,k,l) << "  ";
                    }
                 std::cout << std::endl;
                 }
              std::cout << std::endl;
              }
           std::cout << std::endl;
           }

           std::cout << "ve2 = " << std::endl;
           for (int l=0; l<nY; ++l){
              for (int k=0; k<nY; ++k){
                for (int i=0; i< msh->nen; ++i){
                   for (int j=0; j< msh->nen; ++j){
                      std::cout << ve2(i,j,k,l) << "  ";
                    }
                 std::cout << std::endl;
                 }
              std::cout << std::endl;
              }
           std::cout << std::endl;
           }
        }        

        //Compute omega
        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int r3=0; r3 < nY; ++r3){
                 for (int r4=0; r4 < nY; ++r4){
                    for (int n1=0; n1 < msh->nen; ++n1){
                       for (int n2=0; n2 < msh->nen; ++n2){
                          for (int n3=0; n3 < msh->nen; ++n3){
                             for (int n4=0; n4 < msh->nen; ++n4){
                                //for (int s=0; s < msh->nsd; ++s){
                                //   for (int t=0; t < msh->nsd; ++t){
                                //      omega(r1,r2,r3,r4)+=0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4)*ve(s,n1,r1)*ve(s,n2,r2)*ve(t,n3,r3)*ve(t,n4,r4);
                                //      omega(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n1,n3,n2,n4)*ve(s,n1,r1)*ve(s,n2,r2)*ve(t,n3,r3)*ve(t,n4,r4);
                                //    }
                                //}
                                omega(r1,r2,r3,r4)+=0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4)*(ve2(n1,n2,r1,r2)*ve2(n3,n4,r3,r4));
                                omega(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n1,n3,n2,n4)*(ve2(n1,n2,r1,r2)*ve2(n3,n4,r3,r4));
                             }
                          }
                       }
                    }
                 }
              }
           }
        }

        //Compute omega
        /*for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int r4=0; r4 < nY; ++r4)
                    for (int n1=0; n1 < msh->nen; ++n1)
                       for (int n2=0; n2 < msh->nen; ++n2)
                          for (int n3=0; n3 < msh->nen; ++n3)
                             for (int n4=0; n4 < msh->nen; ++n4)
                                omega(r1,r2,r3,r4)+=Db(n1,n2,n3,n4)*ve2(n1,n2,r1,r2)*ve2(n3,n4,r3,r4);
                                //omega(r1,r2,r3,r4)+=Db(n1,n2,n3,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
       */
     }


    delete[] N;

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    delete[] Xe;
    delete[] dXdZ;
    delete[] dZdX;

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
    delete[] dN;
    delete[] dNdX;
}

void  PrecompRomQuant(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_lam, array<double> beta_lam, array<double> gamma_lam, array<double> omega_lam, array<double> alpha_mu, array<double> beta_mu, array<double> gamma_mu, array<double> omega_mu) {


    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];
    size_t dims[4];

    //Entries of ROB on current element
    dims[0]=msh->nsd; dims[1]=msh->nen; dims[2]=nY;
    array<double> ve(3,dims);

    //Pre-computed FEM matrix ("A"); see notes
    dims[0]=msh->nsd; dims[1]=msh->nsd; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> A(4,dims); 

    //Pre-computed FEM matrix ("B"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> Bb(2,dims); 

    //Pre-computed FEM matrix ("C"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nsd;
    array<double> C(4,dims); 

    //Pre-computed FEM matrix ("D"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> D(4,dims); 

    //Set vector to contain ubar_iI*ubar_iJ
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> u2(2,dims);

    //Set vector to contain ve2
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=nY; dims[3]=nY;
    array<double> ve2(4,dims);
    array<double> dNdN(2,dims);

    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    double** dN = new double* [msh->nen];
    double** dNdX = new double* [msh->nen];
    for (int i=0; i < msh->nen; ++i) 
     { dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    double J;

    alpha_lam.zeroout();  alpha_mu.zeroout();
    beta_lam.zeroout();   beta_mu.zeroout();
    gamma_lam.zeroout();  gamma_mu.zeroout();
    omega_lam.zeroout();  omega_mu.zeroout();

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              for (int k=0; k < nY; ++k) {
                 if ((msh->LM)(e, i+msh->nsd*j) > -1)
                    ve(i,j,k)=V( (msh->LM)( e, i+msh->nsd*j ) , k );
                 else
                    ve(i,j,k)=0.0;
              }
           }
        }
        //Pre-compute terms used in computation of alpha, beta, gamma, and omega
        u2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1){
           for (int n2=0; n2<msh->nen; ++n2){
              for (int i=0; i < msh->nsd; ++i){ 
                 if ((msh->LM)(e, i+msh->nsd*n1) > -1 && (msh->LM)(e, i+msh->nsd*n2) > -1)
                    u2(n1,n2)+=ub( (msh->LM)( e, i+msh->nsd*n1) )*ub( (msh->LM)( e, i+msh->nsd*n2) );
              }
           }
        }
 
        ve2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1)
           for (int n2=0; n2<msh->nen; ++n2)
              for (int r1=0; r1 < nY; ++r1)
                 for (int r2=0; r2 < nY; ++r2)
                    for (int s=0; s < msh->nsd; ++s)
                       ve2(n1,n2,r1,r2)+=ve(s,n1,r1)*ve(s,n2,r2);

        A.zeroout();
        C.zeroout();
        D.zeroout();

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 dNdX[i][j]=0.0;
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];

           //Compute dNdX*dNdX'
           dNdN.zeroout();
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int s=0; s < msh->nsd; ++s)
                    dNdN(n1,n2)+=dNdX[n1][s]*dNdX[n2][s];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int s1=0; s1 < msh->nsd; ++s1)
              for (int s2=0; s2 < msh->nsd; ++s2)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       A(s1,s2,n1,n2) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s1]*dNdX[n2][s2];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int s=0; s < msh->nsd; ++s)
                       C(n1,n2,n3,s) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s]*dNdN(n2,n3);
         
           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
               for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int n4=0; n4 < msh->nen; ++n4)
                       D(n1,n2,n3,n4) += (1.0/6.0)*J*(quad->W(k))*dNdN(n1,n2)*dNdN(n3,n4);
        }

        Bb.zeroout();
        //Compute Bbar from definition of A
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int s=0; s < msh->nsd; ++s)
                 Bb(n1,n2) += (phys[msh->mat(e)].mu)*A(s,s,n1,n2);

        //Compute alpha
        for (int r1=0; r1 < nY; ++r1){
           alpha_lam(r1)=0.0;
           alpha_mu(r1)=0.0;
        }
 
        //Compute beta
        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int n1=0; n1 < msh->nen; ++n1){
                 for (int n2=0; n2 < msh->nen; ++n2){
                    for (int s1=0; s1 < msh->nsd; ++s1){
                       for (int s2=0; s2 < msh->nsd; ++s2){
                          beta_lam(r1,r2)+=(phys[msh->mat(e)].lam)*A(s1,s2,n1,n2)*ve(s1,n1,r1)*ve(s2,n2,r2);
                          beta_mu(r1,r2)+=(phys[msh->mat(e)].mu)*A(s1,s2,n2,n1)*ve(s1,n1,r1)*ve(s2,n2,r2);
                       }
                    }
                 }
              }
           }
        }

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    beta_mu(r1,r2)+=Bb(n1,n2)*ve2(n1,n2,r1,r2);
             
        //Compute gamma  
        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int r3=0; r3 < nY; ++r3){
                 for (int n1=0; n1 < msh->nen; ++n1){
                    for (int n2=0; n2 < msh->nen; ++n2){
                       for (int n3=0; n3 < msh->nen; ++n3){
                          for (int s=0; s < msh->nsd; ++s){
                             gamma_lam(r1,r2,r3)+=0.5*(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*ve2(n2,n3,r2,r3)*ve(s,n1,r1);
                             gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n3,n2,n1,s)*ve2(n2,n3,r2,r3)*ve(s,n1,r1);
                             //gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n1,n2,n3,s)*ve2(n2,n1,r2,r3)*ve(s,n3,r1);
                          }
                       }
                    }
                 }
              }
           }
        }

        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int r3=0; r3 < nY; ++r3){
                 for (int n1=0; n1 < msh->nen; ++n1){
                    for (int n2=0; n2 < msh->nen; ++n2){
                       for (int n3=0; n3 < msh->nen; ++n3){
                          for (int s=0; s < msh->nsd; ++s){
                             gamma_lam(r1,r2,r3)+=(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*ve2(n2,n3,r1,r2)*ve(s,n1,r3);
                             gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*(C(n2,n1,n3,s)+C(n3,n1,n2,s))*ve2(n2,n3,r1,r2)*ve(s,n1,r3);
                             //gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*(C(n1,n2,n3,s)+C(n3,n2,n1,s))*ve2(n1,n3,r1,r2)*ve(s,n2,r3);
                          }
                       }
                    }
                 }
              }
           }
        } 
 
        //Compute omega
        for (int r1=0; r1 < nY; ++r1){
           for (int r2=0; r2 < nY; ++r2){
              for (int r3=0; r3 < nY; ++r3){
                 for (int r4=0; r4 < nY; ++r4){
                    for (int n1=0; n1 < msh->nen; ++n1){
                       for (int n2=0; n2 < msh->nen; ++n2){
                          for (int n3=0; n3 < msh->nen; ++n3){
                             for (int n4=0; n4 < msh->nen; ++n4){
                                omega_lam(r1,r2,r3,r4)+=0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
                                omega_mu(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n2,n3,n1,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
                                //omega_mu(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n2,n3,n1,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
                             }
                          }
                       }
                    }
                 }
              }
           }
        }
     }


    delete[] N;

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    delete[] Xe;
    delete[] dXdZ;
    delete[] dZdX;

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
    delete[] dN;
    delete[] dNdX;
}

void  PrecompRomQuant_lam(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_lam, array<double> beta_lam, array<double> gamma_lam, array<double> omega_lam) {

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];

    size_t dims[4];

    //Entries of ROB on current element
    dims[0]=msh->nsd; dims[1]=msh->nen; dims[2]=nY;
    array<double> ve(3,dims);

    //Pre-computed FEM matrix ("A"); see notes
    dims[0]=msh->nsd; dims[1]=msh->nsd; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> A(4,dims); 

    //Pre-computed FEM matrix ("C"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nsd;
    array<double> C(4,dims); 

    //Pre-computed FEM matrix ("D"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> D(4,dims); 

    //Set vector to contain ubar_iI*ubar_iJ
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> u2(2,dims);

    //Set vector to contain ve2
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=nY; dims[3]=nY;
    array<double> ve2(4,dims);
    array<double> dNdN(2,dims);

    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    double** dN = new double* [msh->nen];
    double** dNdX = new double* [msh->nen];
    for (int i=0; i < msh->nen; ++i) 
     { dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    double J;

    alpha_lam.zeroout();
    beta_lam.zeroout();
    gamma_lam.zeroout();
    omega_lam.zeroout();

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              for (int k=0; k < nY; ++k) {
                 if ((msh->LM)(e, i+msh->nsd*j) > -1)
                    ve(i,j,k)=V( (msh->LM)( e, i+msh->nsd*j ) , k );
                 else
                    ve(i,j,k)=0.0;
              }
           }
        }
        //Pre-compute terms used in computation of alpha, beta, gamma, and omega
        u2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1){
           for (int n2=0; n2<msh->nen; ++n2){
              for (int i=0; i < msh->nsd; ++i){ 
                 if ((msh->LM)(e, i+msh->nsd*n1) > -1 && (msh->LM)(e, i+msh->nsd*n2) > -1)
                    u2(n1,n2)+=ub( (msh->LM)( e, i+msh->nsd*n1) )*ub( (msh->LM)( e, i+msh->nsd*n2) );
              }
           }
        }
 
        ve2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1)
           for (int n2=0; n2<msh->nen; ++n2)
              for (int r1=0; r1 < nY; ++r1)
                 for (int r2=0; r2 < nY; ++r2)
                    for (int s=0; s < msh->nsd; ++s)
                       ve2(n1,n2,r1,r2)+=ve(s,n1,r1)*ve(s,n2,r2);

        A.zeroout();
        C.zeroout();
        D.zeroout();

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 dNdX[i][j]=0.0;
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];

           //Compute dNdX*dNdX'
           dNdN.zeroout();
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int s=0; s < msh->nsd; ++s)
                    dNdN(n1,n2)+=dNdX[n1][s]*dNdX[n2][s];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int s1=0; s1 < msh->nsd; ++s1)
              for (int s2=0; s2 < msh->nsd; ++s2)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       A(s1,s2,n1,n2) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s1]*dNdX[n2][s2];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int s=0; s < msh->nsd; ++s)
                       C(n1,n2,n3,s) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s]*dNdN(n2,n3);
         
           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
               for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int n4=0; n4 < msh->nen; ++n4)
                       D(n1,n2,n3,n4) += (1.0/6.0)*J*(quad->W(k))*dNdN(n1,n2)*dNdN(n3,n4);
        }

        //Compute alpha
        for (int r1=0; r1 < nY; ++r1)
           alpha_lam(r1)=0.0;
 
        //Compute beta
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    for (int s1=0; s1 < msh->nsd; ++s1)
                       for (int s2=0; s2 < msh->nsd; ++s2)
                          beta_lam(r1,r2)+=(phys[msh->mat(e)].lam)*A(s1,s2,n1,n2)*ve(s1,n1,r1)*ve(s2,n2,r2);

        //Compute gamma  
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma_lam(r1,r2,r3)+=0.5*(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*ve2(n2,n3,r2,r3)*ve(s,n1,r1);

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma_lam(r1,r2,r3)+=(phys[msh->mat(e)].lam)*C(n1,n2,n3,s)*ve2(n2,n3,r1,r2)*ve(s,n1,r3);
 
        //Compute omega
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int r4=0; r4 < nY; ++r4)
                    for (int n1=0; n1 < msh->nen; ++n1)
                       for (int n2=0; n2 < msh->nen; ++n2)
                          for (int n3=0; n3 < msh->nen; ++n3)
                             for (int n4=0; n4 < msh->nen; ++n4)
                                omega_lam(r1,r2,r3,r4)+=0.5*(phys[msh->mat(e)].lam)*D(n1,n2,n3,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
     }

    delete[] N;

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    delete[] Xe;
    delete[] dXdZ;
    delete[] dZdX;

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
    delete[] dN;
    delete[] dNdX;
}

void  PrecompRomQuant_mu(array<double> V, array<double> ub, int nY, mesh* msh, physics* phys, DirichletBCs* dbc, integration* quad, array<double> alpha_mu, array<double> beta_mu, array<double> gamma_mu, array<double> omega_mu) {

    //Make arrays for elemental positions (reference) and displacements
    double* N   = new double [msh->nen];

    size_t dims[4];

    //Entries of ROB on current element
    dims[0]=msh->nsd; dims[1]=msh->nen; dims[2]=nY;
    array<double> ve(3,dims);

    //Pre-computed FEM matrix ("A"); see notes
    dims[0]=msh->nsd; dims[1]=msh->nsd; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> A(4,dims); 

    //Pre-computed FEM matrix ("B"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> Bb(2,dims); 

    //Pre-computed FEM matrix ("C"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nsd;
    array<double> C(4,dims); 

    //Pre-computed FEM matrix ("D"); see notes
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=msh->nen; dims[3]=msh->nen;
    array<double> D(4,dims); 

    //Set vector to contain ubar_iI*ubar_iJ
    dims[0]=msh->nen; dims[1]=msh->nen;
    array<double> u2(2,dims);

    //Set vector to contain ve2
    dims[0]=msh->nen; dims[1]=msh->nen; dims[2]=nY; dims[3]=nY;
    array<double> ve2(4,dims);
    array<double> dNdN(2,dims);

    double** Xe = new double* [msh->nsd];
    double** dXdZ = new double* [msh->nsd];
    double** dZdX = new double* [msh->nsd];
    for (int i=0; i < msh->nsd; ++i)
       { Xe[i] = new double [msh->nen];
         dXdZ[i] =  new double [msh->nsd];
         dZdX[i] =  new double [msh->nsd]; }
         
    double** dN = new double* [msh->nen];
    double** dNdX = new double* [msh->nen];
    for (int i=0; i < msh->nen; ++i) 
     { dN[i] = new double [msh->nsd];
       dNdX[i] = new double [msh->nsd]; }

    double J;

    alpha_mu.zeroout();
    beta_mu.zeroout();
    gamma_mu.zeroout();
    omega_mu.zeroout();

    for (int e=0; e < msh->nel; ++e){
        //Extract element positions (reference) and displacements,
        //taking DBCs into account
        for (int i=0; i < msh->nsd; ++i) {
           for (int j=0; j < msh->nen; ++j) {
              Xe[i][j]=(msh->X)( msh->IEN(e, j), i);
              for (int k=0; k < nY; ++k) {
                 if ((msh->LM)(e, i+msh->nsd*j) > -1)
                    ve(i,j,k)=V( (msh->LM)( e, i+msh->nsd*j ) , k );
                 else
                    ve(i,j,k)=0.0;
              }
           }
        }
        //Pre-compute terms used in computation of alpha, beta, gamma, and omega
        u2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1){
           for (int n2=0; n2<msh->nen; ++n2){
              for (int i=0; i < msh->nsd; ++i){ 
                 if ((msh->LM)(e, i+msh->nsd*n1) > -1 && (msh->LM)(e, i+msh->nsd*n2) > -1)
                    u2(n1,n2)+=ub( (msh->LM)( e, i+msh->nsd*n1) )*ub( (msh->LM)( e, i+msh->nsd*n2) );
              }
           }
        }
 
        ve2.zeroout();
        for (int n1=0; n1<msh->nen; ++n1)
           for (int n2=0; n2<msh->nen; ++n2)
              for (int r1=0; r1 < nY; ++r1)
                 for (int r2=0; r2 < nY; ++r2)
                    for (int s=0; s < msh->nsd; ++s)
                       ve2(n1,n2,r1,r2)+=ve(s,n1,r1)*ve(s,n2,r2);

        A.zeroout();
        C.zeroout();
        D.zeroout();

        for (int k=0; k< quad->nint; ++k){
           //Compute shape functions
           ShapeFunctions(quad, k, N, dN);

           //Compute jacobian of map      
           for (int i=0; i < msh->nsd; ++i) {
              for (int j=0; j < msh->nsd; ++j) {
                 dXdZ[i][j] = 0.0;
                 for (int p=0; p < msh->nen; ++p)
                    dXdZ[i][j]+=Xe[i][p]*dN[p][j];
              }
           }
           
           //Compute dZdX = inv(dXdZ)
           inverse3x3(dXdZ,dZdX);
           //Compute J = det(dXdZ)
           J=determinant3x3(dXdZ);

           //Compute dNdX
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 dNdX[i][j]=0.0;
           for (int i=0; i < msh->nen; ++i)
              for (int j=0; j < msh->nsd; ++j)
                 for (int p=0; p < msh->nsd; ++p)
                    dNdX[i][j]+=dN[i][p]*dZdX[p][j];

           //Compute dNdX*dNdX'
           dNdN.zeroout();
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int s=0; s < msh->nsd; ++s)
                    dNdN(n1,n2)+=dNdX[n1][s]*dNdX[n2][s];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int s1=0; s1 < msh->nsd; ++s1)
              for (int s2=0; s2 < msh->nsd; ++s2)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       A(s1,s2,n1,n2) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s1]*dNdX[n2][s2];

           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
              for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int s=0; s < msh->nsd; ++s)
                       C(n1,n2,n3,s) += (1.0/6.0)*J*(quad->W(k))*dNdX[n1][s]*dNdN(n2,n3);
         
           //Pre-compute terms used in computation of alpha, beta, gamma, and omega
           for (int n1=0; n1 < msh->nen; ++n1)
               for (int n2=0; n2 < msh->nen; ++n2)
                 for (int n3=0; n3 < msh->nen; ++n3)
                    for (int n4=0; n4 < msh->nen; ++n4)
                       D(n1,n2,n3,n4) += (1.0/6.0)*J*(quad->W(k))*dNdN(n1,n2)*dNdN(n3,n4);
        }

        Bb.zeroout();
        //Compute Bbar from definition of A
        for (int n1=0; n1 < msh->nen; ++n1)
           for (int n2=0; n2 < msh->nen; ++n2)
              for (int s=0; s < msh->nsd; ++s)
                 Bb(n1,n2) += (phys[msh->mat(e)].mu)*A(s,s,n1,n2);

        //Compute alpha
        for (int r1=0; r1 < nY; ++r1)
           alpha_mu(r1)=0.0;
 
        //Compute beta
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    for (int s1=0; s1 < msh->nsd; ++s1)
                       for (int s2=0; s2 < msh->nsd; ++s2)
                          beta_mu(r1,r2)+=(phys[msh->mat(e)].mu)*A(s1,s2,n1,n2)*ve(s1,n1,r1)*ve(s2,n2,r2);

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int n1=0; n1 < msh->nen; ++n1)
                 for (int n2=0; n2 < msh->nen; ++n2)
                    beta_mu(r1,r2)+=Bb(n1,n2)*ve2(n1,n2,r1,r2);
             
        //Compute gamma  
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*C(n1,n2,n3,s)*ve2(n2,n1,r2,r3)*ve(s,n3,r1);

        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int n1=0; n1 < msh->nen; ++n1)
                    for (int n2=0; n2 < msh->nen; ++n2)
                       for (int n3=0; n3 < msh->nen; ++n3)
                          for (int s=0; s < msh->nsd; ++s)
                             gamma_mu(r1,r2,r3)+=(phys[msh->mat(e)].mu)*(C(n1,n2,n3,s)+C(n3,n2,n1,s))*ve2(n1,n3,r1,r2)*ve(s,n2,r3);
 
        //Compute omega
        for (int r1=0; r1 < nY; ++r1)
           for (int r2=0; r2 < nY; ++r2)
              for (int r3=0; r3 < nY; ++r3)
                 for (int r4=0; r4 < nY; ++r4)
                    for (int n1=0; n1 < msh->nen; ++n1)
                       for (int n2=0; n2 < msh->nen; ++n2)
                          for (int n3=0; n3 < msh->nen; ++n3)
                             for (int n4=0; n4 < msh->nen; ++n4)
                                omega_mu(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n2,n3,n1,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
                                //omega_mu(r1,r2,r3,r4)+=(phys[msh->mat(e)].mu)*D(n1,n2,n3,n4)*ve2(n4,n3,r1,r4)*ve2(n1,n2,r2,r3);
     }

    delete[] N;

    for (int i=0; i < msh->nsd; ++i)
       { delete[] Xe[i];
         delete[] dXdZ[i];
         delete[] dZdX[i]; }
    delete[] Xe;
    delete[] dXdZ;
    delete[] dZdX;

    for (int i=0; i < msh->nen; ++i)
       { delete[] dN[i];
         delete[] dNdX[i]; }
    delete[] dN;
    delete[] dNdX;
}
