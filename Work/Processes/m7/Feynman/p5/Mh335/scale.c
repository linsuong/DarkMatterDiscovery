#include<math.h>
#include<stdlib.h>
#include"/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/nType.h"
#define min(x,y) (x<y? x:y)
#define max(x,y) (x>y? x:y)
extern void ScaleCC(int nsub, REAL*, double (*calcPV)(char,char*,double*), double*,double*,double*,double*,double *);
void ScaleCC(int nsub,REAL*modelVal, double (*calcPV)(char,char*,double*), double *pvect,double *Qren, double *Qpdf1,double *Qpdf2, double*Qshow)
{ double                
  double ss,s1; int i;
