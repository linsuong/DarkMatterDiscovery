#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/extern.h"
#include "/Users/linusong/Documents/DMDiscovery/calchep_3.9.1/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=20;
static ModelPrtclsStr ModelPrtcls_[20]=
{
  {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","wZ",2,1,3,0}
, {"G","G",1, 21, "0","0",2,8,16,0}
, {"W+","W-",0, 24, "MW","wW",2,1,3,3}
, {"n1","N1",0, 12, "0","0",1,1,1,0}
, {"e1","E1",0, 11, "0","0",1,1,2,-3}
, {"n2","N2",0, 14, "0","0",1,1,1,0}
, {"e2","E2",0, 13, "Mm","0",1,1,2,-3}
, {"n3","N3",0, 16, "0","0",1,1,1,0}
, {"e3","E3",0, 15, "Mtau","0",1,1,2,-3}
, {"u","U",0, 2, "Mu","0",1,3,6,2}
, {"d","D",0, 1, "Md","0",1,3,6,-1}
, {"c","C",0, 4, "Mcp","0",1,3,6,2}
, {"s","S",0, 3, "Ms","0",1,3,6,-1}
, {"t","T",0, 6, "Mtop","wtop",1,3,6,2}
, {"b","B",0, 5, "Mbp","0",1,3,6,-1}
, {"H","H",1, 25, "MH","wH",0,1,1,0}
, {"~h1","~h1",1, 3000022, "Mh1","wh1",0,1,1,0}
, {"~h2","~h2",1, 3000027, "Mh2","wh2",0,1,1,0}
, {"~h+","~h-",0, 3000025, "Mhc","wHC",0,1,1,3}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=27;
int nModelFunc=45;
static int nCurrentVars=26;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[72]={
 "alphaSMZ","EE","Q","Mm","Mtau","Ms","McMc","MbMb","Mtop","MH"
,"MZ","MW","wtop","wZ","wW","s12","s23","s13","Maux","ld345"
,"MD1","DMP","DM3","ld","Mu","Md","wh1","Mcp","Mbp","alphaE0"
,"CW","SW","GF","vv","c12","c23","c13","LamQCD","Mb","Mt"
,"Mc","Mh1","Mhc","Mh2","md2","ld3","ld4","ld5","lam","aQCD"
,"ahF_c","ahF_b","ahF_t","ahF_e3","a_hV_W","a_hS_Hc","aQCD_h","Rqcd_h","Quq","Qdq"
,"LGGH","LAAH","B00000","B00001","B00002","B00003","B00004","B00005","B00006","B00007"
,"B00008","B00009"};
char**varNames=varNames_;
static REAL varValues_[72]={
   1.184000E-01,  3.134300E-01,  1.000000E+02,  1.057000E-01,  1.777000E+00,  2.000000E-01,  1.230000E+00,  4.250000E+00,  1.725000E+02,  1.250000E+02
,  9.118800E+01,  8.038500E+01,  1.590000E+00,  2.494440E+00,  2.088950E+00,  2.210000E-01,  4.000000E-02,  3.500000E-03,  1.000000E+00,  0.000000E+00
,  1.000000E+02,  1.000000E+01,  1.000000E+01,  1.000000E+00,  3.000000E-03,  5.000000E-03,  0.000000E+00};
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
   nCurrentVars=27;
   V[27]=V[6]*(1+4/(double)((3))*alphaQCD(V[6])/(M_PI));
   if(!isfinite(V[27]) || FError) return 27;
   nCurrentVars=28;
   V[28]=V[7]*(1+4/(double)((3))*alphaQCD(V[7])/(M_PI));
   if(!isfinite(V[28]) || FError) return 28;
   nCurrentVars=29;
   V[29]=1/(137.036);
   if(!isfinite(V[29]) || FError) return 29;
   nCurrentVars=30;
   V[30]=V[11]/(V[10]);
   if(!isfinite(V[30]) || FError) return 30;
   nCurrentVars=31;
   V[31]=Sqrt(1-Pow(V[30],2));
   if(!isfinite(V[31]) || FError) return 31;
   nCurrentVars=32;
   V[32]=Pow(V[1],2)/(Pow(2*V[31]*V[11],2))/(M_SQRT2);
   if(!isfinite(V[32]) || FError) return 32;
   nCurrentVars=33;
   V[33]=2*V[11]/(V[1])*V[31];
   if(!isfinite(V[33]) || FError) return 33;
   nCurrentVars=34;
   V[34]=Sqrt(1-Pow(V[15],2));
   if(!isfinite(V[34]) || FError) return 34;
   nCurrentVars=35;
   V[35]=Sqrt(1-Pow(V[16],2));
   if(!isfinite(V[35]) || FError) return 35;
   nCurrentVars=36;
   V[36]=Sqrt(1-Pow(V[17],2));
   if(!isfinite(V[36]) || FError) return 36;
   nCurrentVars=37;
   V[37]=initQCD5(V[0],V[6],V[7],V[8]);
   if(!isfinite(V[37]) || FError) return 37;
 FirstQ:
 cErr=1;
   nCurrentVars=38;
   V[38]=MbEff(V[2]);
   if(!isfinite(V[38]) || FError) return 38;
   nCurrentVars=39;
   V[39]=MtEff(V[2]);
   if(!isfinite(V[39]) || FError) return 39;
   nCurrentVars=40;
   V[40]=McEff(V[2]);
   if(!isfinite(V[40]) || FError) return 40;
   nCurrentVars=41;
   V[41]=V[20];

   nCurrentVars=42;
   V[42]=V[41]+V[21];
   if(!isfinite(V[42]) || FError) return 42;
   nCurrentVars=43;
   V[43]=V[42]+V[22];
   if(!isfinite(V[43]) || FError) return 43;
   nCurrentVars=44;
   V[44]=1/(double)((2))*V[19]*Pow(V[33],2)-Pow(V[41],2);
   if(!isfinite(V[44]) || FError) return 44;
   nCurrentVars=45;
   V[45]=2/(Pow(V[33],2))*(Pow(V[42],2)+1/(double)((2))*V[19]*Pow(V[33],2)-Pow(V[41],2));
   if(!isfinite(V[45]) || FError) return 45;
   nCurrentVars=46;
   V[46]=1/(Pow(V[33],2))*(Pow(V[43],2)+Pow(V[41],2)-2*Pow(V[42],2));
   if(!isfinite(V[46]) || FError) return 46;
   nCurrentVars=47;
   V[47]=1/(Pow(V[33],2))*(Pow(V[41],2)-Pow(V[43],2));
   if(!isfinite(V[47]) || FError) return 47;
   nCurrentVars=48;
   V[48]=Pow(V[1]/(V[31])*V[9]/(V[11]),2)/(8);
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=alphaQCD(V[9])/(M_PI);
   if(!isfinite(V[49]) || FError) return 49;
   nCurrentVars=50;
   V[50]=-V[1]/(V[11])*V[40]/(V[31])/(2)/(V[27]);
   if(!isfinite(V[50]) || FError) return 50;
   nCurrentVars=51;
   V[51]=-V[1]/(V[11])*V[38]/(V[31])/(2)/(V[28]);
   if(!isfinite(V[51]) || FError) return 51;
   nCurrentVars=52;
   V[52]=-V[1]/(V[11])*V[8]/(V[31])/(2)/(V[8]);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=-V[1]/(V[11])*V[4]/(V[31])/(2)/(V[4]);
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=V[1]*V[11]/(V[31])/(Pow(V[11],2));
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=-2/(V[1])*V[11]*V[31]*V[45]/(Pow(V[42],2));
   if(!isfinite(V[55]) || FError) return 55;
   nCurrentVars=56;
   V[56]=alphaQCD(V[9])/(M_PI);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=Sqrt(1+V[56]*(149/(double)((12))+V[56]*(68.6482-V[56]*212.447)));
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=4/(double)((9));
   if(!isfinite(V[58]) || FError) return 58;
   nCurrentVars=59;
   V[59]=1/(double)((9));
   if(!isfinite(V[59]) || FError) return 59;
   nCurrentVars=60;
   V[60]=-Cabs(hGGeven(V[9],V[56],3,1,3,V[28],V[51],1,3,V[27],V[50],1,3,V[8],V[52]));
   if(!isfinite(V[60]) || FError) return 60;
   nCurrentVars=61;
   V[61]=-Cabs(V[59]*hAAeven(V[9],V[56],1,1,3,V[28],V[51])+V[58]*hAAeven(V[9],V[56],2,1,3,V[8],V[52],1,3,V[27],V[50])+hAAeven(V[9],V[56],3,1,1,V[4],V[53],2,1,V[11],V[54],0,1,V[42],V[55]));
   if(!isfinite(V[61]) || FError) return 61;
   nCurrentVars=62;
   V[62]=1+2*Pow(V[30],2);
   if(!isfinite(V[62]) || FError) return 62;
   nCurrentVars=63;
   V[63]=1-4*Pow(V[30],2);
   if(!isfinite(V[63]) || FError) return 63;
   nCurrentVars=64;
   V[64]=1-2*Pow(V[30],2);
   if(!isfinite(V[64]) || FError) return 64;
   nCurrentVars=65;
   V[65]=Pow(V[1],2)*V[33]-8*V[60]*Pow(V[11],2)*V[57]*Pow(V[31],2);
   if(!isfinite(V[65]) || FError) return 65;
   nCurrentVars=66;
   V[66]=V[45]+V[46]+V[47];
   if(!isfinite(V[66]) || FError) return 66;
   nCurrentVars=67;
   V[67]=V[45]+V[46]-V[47];
   if(!isfinite(V[67]) || FError) return 67;
   nCurrentVars=68;
   V[68]=V[46]+V[47];
   if(!isfinite(V[68]) || FError) return 68;
   nCurrentVars=69;
   V[69]=V[46]-V[47];
   if(!isfinite(V[69]) || FError) return 69;
   nCurrentVars=70;
   V[70]=1-4*Pow(V[30],2)+4*Pow(V[30],4);
   if(!isfinite(V[70]) || FError) return 70;
   nCurrentVars=71;
   V[71]=V[45]+V[46];
   if(!isfinite(V[71]) || FError) return 71;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
