
#include"micromegas.h"
#include"micromegas_aux.h"
#include"../CalcHEP_src/c_source/strfun/include/LHA2C.h" 

static void*pdf=NULL;
static char ** setArr=NULL;


int  LHAPDFList(void)
{ if(!setArr) setArr=listLHA();   
  if(!setArr  || !setArr[0] ) { printf("LHAPDF is not available\n"); return 1;}
  for(int n=0;setArr[n];n++) 
  { void *p=initpdfLHA(setArr[n],0);
    int ID,maxNum;
    pdfInfoLHA(p, &ID, &maxNum, NULL,NULL,NULL,NULL,NULL);
    printf(" setName: %-25.25s  ID=%-6d  nIitems=%d\n", setArr[n], ID,maxNum);
    deletePdfLHA(p);
  }  
}    

static double lhaPartonPdf(int pdg, double x,double q) { return pdfFuncLHA(pdf, pdg, x, q);}

static double lhaAlpha(double q) { return alphaLHA(pdf,q);}

int  setLHAPDF(char *name, int num )
{  
   if(pdf){  deletePdfLHA(pdf); pdf=NULL;}
   if(!setArr) setArr=listLHA();
   if(!setArr  || !setArr[0] ) return 1;
   for(int n=0;setArr[n];n++) if(strcmp(name, setArr[n])==0)
   {  pdf=initpdfLHA(setArr[n],num); 
      parton_distr=lhaPartonPdf;
      parton_alpha=lhaAlpha;
      return 0;
   }
   return 1;
}
   

//  FORTRAN 

extern int  setlhapdf_(char *fname, int* num, int len);

int  setlhapdf_(char *fname, int*num, int len)
{ 
  char cname[50];
  fName2c(fname,cname,len);
  return setLHAPDF(cname, *num);
}  


extern int lhapdflist_(void);

int  lhapdflist_(void) { return LHAPDFList(); }
