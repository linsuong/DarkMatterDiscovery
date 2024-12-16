
#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

extern REAL MassArray(int id,  int i);
extern REAL MixMatrix(int id, int i,int j);


int main(void)
{
   REAL  MTestSym[6] = {0., 1., 1.,
                            0., 1.,
                                0.};
printf("Diagonalization of symmetric  matrix:\n"
"          0  1  1\n"
"          1  0  1\n"
"          1  1  0\n" );

   int dim=3;

   printf("\n  For   func*.mdl\n");
   REAL m[3][3];
   m[0][0]=MTestSym[0];
   m[0][1]=MTestSym[1];
   m[0][2]=MTestSym[2];
   m[1][1]=MTestSym[3];
   m[1][2]=MTestSym[4];
   m[2][2]=MTestSym[5];

   initDiagonal();
   int id=rDiagonal(dim,m[0][0],m[0][1],m[0][2],m[1][1],m[1][2],m[2][2]);

   printf("Result: masses"); for(int i=1;i<=dim;i++)  printf(" %E ",(double)MassArray(id,i)); printf("\n");
   printf("Rotation matrix:\n");
   for(int i=1;i<=3;i++)
   { for(int j=1;j<=3;j++) printf("%12.4E  ",(double)MixMatrix(id,i,j)); printf("\n");}

   printf("\n Restoration of initial matrix\n");

   for(int i=1;i<=dim;i++) for(int j=1; j<=dim;j++)
   { REAL mm=0;
     for(int k=1;k<=dim;k++) mm+=MixMatrix(id,k,j)*MassArray(id,k)*MixMatrix(id,k, i);
      printf(" %12.4E ", (double) mm); if(j==dim) printf("\n");
   }
// Below is an example of the function to diagonalize the mass matric MTestSym
// when not working with a matrix defined in a model file func*.mdl
printf("\n  When not using model files :\n");

   REAL E[dim],  V[dim*dim];
   rJacobi(MTestSym,dim, E , V);
   printf("Result  Masses: "); for(int k=0;k<dim;k++) printf("%E ",(double)E[k]); printf("\n");
    printf("Rotation matrix:\n");

    for(int i=0;i<3;i++)
    { for(int j=0;j<3;j++) printf("%12.4E  ",(double)V[i*dim+j]); printf("\n");}

}
