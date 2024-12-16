#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"


int  main(void)
{
  double vSigma=2.9979E-26; //  [cm^3/s] = 1pbc
  double  SpA[NZ];

vSigma/=3.586373E+00;
SpectraFlag=2;   
  Mcdm=8;
  basicSpectra(Mcdm,5,0,SpA);  
  double res=DwarfSignal(vSigma,SpA);
  printf("res=%E\n",res);
  printf("DwarfSignal> 1  indicates that the model is excluded.\n");
  
}