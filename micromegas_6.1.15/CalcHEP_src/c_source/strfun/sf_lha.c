/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <pthread.h>

#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"

#include "strfun.h"
#include "lha.h"
#include "alphas2.h"
#include "sf_lha.h"
#include "LHA2C.h" 


char ** pdfSetArr=NULL; //  array of LHAPDF sets filled by listLHA

static void *pdf[2]={NULL,NULL};  // initiated   LHAPDF  PDF* structure    
static double xMin[2]={0,0},xMax[2]={1,1},qMin[2]={1,1},qMax[2]={1E10,1E10};  // limits

static char  setName[2][50]={"", ""};                      // input parameters dor LHAPDF  
static int setNum[2]={0,0};
static int  sgn[2]={1,1}, pnum[2]={0,0};      // sign to incoming proton and pdg number for incoming quark

int mc_lha(int i) { return 2212*sgn[i-1];}

//static pthread_mutex_t alpha_key=PTHREAD_MUTEX_INITIALIZER;

static double alpha_lha1(double q )  { return alphaLHA(pdf[0],q);}  //!
static double alpha_lha2(double q )  { return alphaLHA(pdf[1],q);}  //!

int p_lha(int * pNum) 
{  
  if(!pdfSetArr) pdfSetArr=listLHA();
  if(!pdfSetArr) return 0;
  for(int n=0;pdfSetArr[n];n++)
  {  void*p=initpdfLHA( pdfSetArr[n],0);
     int i;
     for(i=0;pNum[i];i++) if(!pdgInLHA(p,pNum[i])) break ;
     deletePdfLHA(p);
     if(pNum[i]==0) return 1;
  }
  return 0;
}


void n_lha(int i, char *name) 
{  int i1=i-1;
   if(strlen(setName[i1]))  sprintf(name,"LHA:%s:%d:%d",setName[i1],setNum[i1],sgn[i1]);
   else  strcpy(name,"LHA:"); 
}


int init_lha(int i,double * be, double * mass) 
{ 
   *mass=0.9383;
   *be=1;

   int k;
   int N;
   pinf_int(Nsub,i,NULL,&N);
   
   if(i==1) sf_alpha[0]=alpha_lha1; else sf_alpha[1]=alpha_lha2;   
   pnum[i-1]=N; 
   return 1;
}

 
static char * lhaMenu(int*pdg)
{  
  if(!pdfSetArr) pdfSetArr=listLHA();
  int len=0,k=0,n;
  for(n=0; pdfSetArr[n];n++) { int l=strlen(pdfSetArr[n]); if(l>len) len=l;}
  len+=2;
  char*menuTxt=malloc(n*len+2);
  
  for(n=0; pdfSetArr[n];n++) 
  {  void*p=initpdfLHA( pdfSetArr[n],0);
     int i;
     for(i=0;pdg[i];i++) if(!pdgInLHA(p,pdg[i])) break ;
     deletePdfLHA(p);
     if(pdg[i]==0)
     { 
       for(int j=0;j<len;j++) 
       if(j<1 || j>strlen(pdfSetArr[n])) menuTxt[1+k*len+j]=' '; else  menuTxt[1+k*len+j]=pdfSetArr[n][j-1];
       k++;
     }      
  }
  if(k==0) { free(menuTxt); return NULL;}  
  menuTxt[1+k*len]=0;
  menuTxt[0]=len;
  return menuTxt; 
}

int r_lha(int i, char *name)
{ int i1=i-1;
  char txt[50];
  char*men,*c;
  int max;

  if(3!=sscanf(name,"LHA:%[^:]:%d:%d",setName[i1],setNum+i1,sgn+i1)) return 0;  
  if(abs(sgn[i1])!=1) return 0;

  if(!pdfSetArr) pdfSetArr=listLHA();
  if(!pdfSetArr) return 0;
  for(int n=0;pdfSetArr[n];n++) 
  if(strcmp(pdfSetArr[n],setName[i1])==0) 
  {   deletePdfLHA(pdf[i1]);
      pdf[i1]=initpdfLHA(pdfSetArr[n],setNum[i1]);
      if(pdf[i1]) 
      {
        pdfInfoLHA(pdf[i1], NULL, NULL, xMin+i1,xMax+i1,qMin+i1,qMax+i1,NULL);
        return 1;
      }
  } 
  return 0;
}

                        
int m_lha(int i,int*pString)
{ 
  void *pscr=NULL;
  void *pscr0=NULL;
  static int n1=0;
  char * strmen=lhaMenu(pString); 
  int i1=i-1;

  if(!strmen) return 0;
  

  int n0=1,k,l;

  if(n1==0 && strlen(setName[i1]))
  { char *ch=strstr(strmen,setName[i1]);
    if(ch) n1= 1+(ch-strmen)/strmen[0];
  }
  menu1(5,10,"LHAlib menu",strmen,"",&pscr,&n1);
  if(n1)
  { char buff[50];
    sscanf(strmen+1+strmen[0]*(n1-1),"%s",buff);
    if(strcmp(buff,setName[i1]))
    { strcpy(setName[i1],buff);
      setNum[i1]=0;
      deletePdfLHA(pdf[i1]);
      pdf[i1]=initpdfLHA(buff,setNum[i1]);
    }
  }  
  else 
  { setName[i1][0]=0;
    setNum[i1]=0;
    sgn[i1]=1;
    return 0;
  } 


  for(;n0!=0 && n0!=3;)
  {  char buff[50];
     int nMax;
     char strmen0[]="\030"
                    " Set = 0                "   
                    " Proton                 "
                    " OK                     ";
     deletePdfLHA(pdf[i1]);
     pdf[i1]=initpdfLHA(setName[i1],0);
     pdfInfoLHA(pdf[i1], NULL, &nMax, xMin+i1,xMax+i1,qMin+i1,qMax+i1,NULL);

     if(nMax>1) improveStr(strmen0,"Set = 0","Set = %d [0,%d]",setNum[i1],nMax-1);
     else       improveStr(strmen0,"Set = 0","Set = 0 (only)");
     
     if(sgn[i1]<0) improveStr(strmen0,"Proton","%s","antiProton");
     
     menu1(5,10,"",strmen0,"",&pscr0,&n0);
     switch(n0) 
     { 
       case 1: if(nMax>1) 
               { correctInt(50,12,"Enter new value ",setNum+i1,1);
                 if(setNum[i1]<0) setNum[i1]=0;
                 if(setNum[i1]>=nMax) setNum[i1]=nMax-1;
                 deletePdfLHA(pdf[i1]);
                 pdf[i1]=initpdfLHA(setName[i1],setNum[i1]);
               }   
               break;
       case 2: sgn[i1]=-sgn[i1]; break;
       case 3: put_text(&pscr0); break;
     }
  }
  
  put_text(&pscr); 
  free(strmen);
  return 1;
}

static pthread_mutex_t strfun_key=PTHREAD_MUTEX_INITIALIZER;

            
double c_lha(int i, double x, double q)
{
  int i1=i-1;
  int p=pnum[i1];
  double z;
  static int nqMax[2]={0,0};
  

//  if(x<xMin[i1]) x=xMin[i1]; else if(x>xMax[i1]) x=xMax[i1];
//  if(q<qMin[i1]) q=qMin[i1]; else 
  if(q>qMax[i1]) 
  {  nqMax[i1]++;
     if(nqMax[i1]==1) printf(" Call of LHA:%s with q=%.2E > Qmax=%.2E\n", setName[i1],q,qMax[i1]);
     else if(nqMax[i1]==100) printf(" More than 100 calls  of LHA:%s with q > Qmax=%.2E\n", setName[i1],qMax[i1]);  
     q=qMax[i1];
  }

  if(sgn[i1]<0 && abs(p)<=6 ) p=-p;
  
//  if(nPROCSS>1 ) pthread_mutex_lock(&strfun_key);
    z= pdfFuncLHA(pdf[i1], p,x,q);   // !!
//  if(nPROCSS>1 ) pthread_mutex_unlock(&strfun_key);  

  if(z<0) z=0;
//  if(z<=0) printf("x=%E q=%E z=%E sc2=%E   \n",x,q,z,sc2); 
//printf("i1=%d pdf=%p pnum=%d x=%E q=%e  %E\n",i1, pdf[i1], pnum[i1],x,q,z);  


  return z;  
}
