#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/extern.h"
#include "/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=17;
static ModelPrtclsStr ModelPrtcls_[17]=
{
  {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","wZ",2,1,3,0}
, {"G","G",1, 21, "0","0",2,8,16,0}
, {"W+","W-",0, 24, "MW","wW",2,1,3,3}
, {"ne","Ne",0, 12, "0","0",1,1,1,0}
, {"e","E",0, 11, "0","0",1,1,2,-3}
, {"nm","Nm",0, 14, "0","0",1,1,1,0}
, {"m","M",0, 13, "0","0",1,1,2,-3}
, {"nl","Nl",0, 16, "0","0",1,1,1,0}
, {"l","L",0, 15, "Ml","0",1,1,2,-3}
, {"u","U",0, 2, "0","0",1,3,6,2}
, {"d","D",0, 1, "0","0",1,3,6,-1}
, {"c","C",0, 4, "Mc","0",1,3,6,2}
, {"s","S",0, 3, "Ms","0",1,3,6,-1}
, {"t","T",0, 6, "Mtp","wt",1,3,6,2}
, {"b","B",0, 5, "Mb","0",1,3,6,-1}
, {"h","h",1, 25, "Mh","wh",0,1,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=10;
int nModelFunc=10;
static int nCurrentVars=9;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[20]={
 "EE","SW","Q","MW","Mtp","McMc","MbMb","alphaSMZ","Ml","Mh"
,"CW","GF","MZ","LamQCD","Mb","Mc","Ms","LAAh","LGGh","aQCDh"
};
char**varNames=varNames_;
static REAL varValues_[20]={
   3.133300E-01,  4.740000E-01,  1.000000E+02,  8.038500E+01,  1.725000E+02,  1.230000E+00,  4.250000E+00,  1.184000E-01,  1.777000E+00,  1.250000E+02
};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   nCurrentVars=10;
   V[10]=Sqrt(1-Pow(V[1],2));
   if(!isfinite(V[10]) || FError) return 10;
   nCurrentVars=11;
   V[11]=Pow(V[0],2)/(Pow(2*V[1]*V[3],2))/(M_SQRT2);
   if(!isfinite(V[11]) || FError) return 11;
   nCurrentVars=12;
   V[12]=V[3]/(V[10]);
   if(!isfinite(V[12]) || FError) return 12;
   nCurrentVars=13;
   V[13]=initQCD5(V[7],V[5],V[6],V[4]);
   if(!isfinite(V[13]) || FError) return 13;
 FirstQ:
 cErr=1;
   nCurrentVars=14;
   V[14]=MbEff(V[2]);
   if(!isfinite(V[14]) || FError) return 14;
   nCurrentVars=15;
   V[15]=McEff(V[2]);
   if(!isfinite(V[15]) || FError) return 15;
   nCurrentVars=16;
   V[16]=MqEff(0.096,V[2]);
   if(!isfinite(V[16]) || FError) return 16;
   nCurrentVars=17;
   V[17]=-Cabs(lAAhiggs(V[9],"h"));
   if(!isfinite(V[17]) || FError) return 17;
   nCurrentVars=18;
   V[18]=-Cabs(lGGhiggs(V[9],"h"));
   if(!isfinite(V[18]) || FError) return 18;
   nCurrentVars=19;
   V[19]=alphaQCD(V[9])/(Acos(-1));
   if(!isfinite(V[19]) || FError) return 19;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
