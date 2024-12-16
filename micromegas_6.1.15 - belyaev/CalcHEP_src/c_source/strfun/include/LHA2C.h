
#ifndef LHA2C
#define LHA2C
  extern void * initpdfLHA(char * name, int num);
  extern void   pdfInfoLHA(void *p,  int*ID, int*maxNum, double*xMin,double*xMax,double*qMin,double *qMax, int *pdg);
  extern int    pdgInLHA( void *p,int pdg);
  extern double pdfFuncLHA(void*p, int pdg, double x, double q);
  extern double alphaLHA(void*p,double q);
  extern void   deletePdfLHA(void *p);
  extern char** listLHA(void);
#endif
