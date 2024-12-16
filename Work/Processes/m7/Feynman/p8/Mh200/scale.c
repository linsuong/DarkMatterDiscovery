#include<math.h>
#include<stdlib.h>
#include"/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/nType.h"
#define min(x,y) (x<y? x:y)
#define max(x,y) (x>y? x:y)
extern void ScaleCC(int nsub, REAL*, double (*calcPV)(char,char*,double*), double*,double*,double*,double*,double *);
void ScaleCC(int nsub,REAL*modelVal, double (*calcPV)(char,char*,double*), double *pvect,double *Qren, double *Qpdf1,double *Qpdf2, double*Qshow)
{ double X[9];          
  double ss,s1; int i;
X[0] = modelVal[42];  /* Mhc  */ 
X[1]=X[0];
 *Qren= X[1];
X[2] = modelVal[42];  /* Mhc  */ 
X[3]=X[2];
 *Qpdf1= X[3];
X[4] = modelVal[42];  /* Mhc  */ 
X[5]=X[4];
 *Qpdf2= X[5];
X[6] = modelVal[42];  /* Mhc  */ 
X[7]=X[6];
 *Qshow= X[7];
}
