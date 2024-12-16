#include<stdlib.h>

#ifndef __Variables_and_Particles__
#define __Variables_and_Particles__

#include "nType.h"

extern int nModelParticles;
 
typedef struct
{ 
   char* name; char* aname;  int selfC; int NPDG; char* mass; char* width; int spin2; int cdim;  int g;  int q3;
}  ModelPrtclsStr;

extern ModelPrtclsStr*ModelPrtcls;
extern int nModelVars;
extern int nModelFunc;
extern int*currentVarPtr;
extern char**varNames;
extern REAL *varValues;
extern int calcMainFunc(void);

/*  VandPgate */

extern double usrFF_(int n_in, int n_out,double * pvect,char**pnames,int*pdg);
extern double usrfun_(char * name,int n_in, int n_out, double * pvect,char**pnames,int*pdg);
extern double usrFF(int n_in, int n_out,double * pvect,char**pnames,int*pdg);
extern double usrfun(char * name,int n_in, int n_out, double * pvect,char**pnames,int*pdg);

extern double alpha_lha(double q );

#endif
