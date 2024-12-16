#include <stddef.h>

#include"include/LHA2C.h"
void * initLHApdf(char * name, int num, int*maxNum, int*ID,  double*xMin,double*xMax,double*qMin,double *qMax)
{return NULL;
}

void * initpdfLHA(char * name, int num){ return 0;}
void   pdfInfoLHA(void *p, int*ID,  int*maxNum,  double*xMin,double*xMax,double*qMin,double *qMax, int *pdg)
{ 
  if(ID) *ID=0;
  if(maxNum) *maxNum=0;
  if(xMin) *xMin=0;
  if(xMax) *xMax=0;
  if(qMin) *qMin=0;
  if(qMax) *qMax=0;
}  
int    pdgInLHA( void *p,int pdg){ return 0;}
double pdfFuncLHA(void*p, int pdg, double x, double q){ return 0;}
double alphaLHA(void*p,double q){ return 0;};
void   deletePdfLHA(void *p){;}
char** listLHA(void){return NULL;}


  