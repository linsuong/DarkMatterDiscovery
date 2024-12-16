#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int main(int argc,char** argv)
{ int err,n,i;
  SpectraFlag=0;      
  double pi0[NZ],pich[NZ],kl[NZ],ks[NZ],kch[NZ];
  double M=0.8;
  basicSpectra(M,111,0,pi0);
  basicSpectra(M,211,0,pich);
  basicSpectra(M,130,0,kl);
  basicSpectra(M,310,0,ks);
  basicSpectra(M,321,0,kch);
  displayPlot(" Edgamma/dE","E",0.001,M ,1,5,"p0",0, eSpectdNdE,pi0
                                            ,"pich",0, eSpectdNdE,pich
                                            ,"KL",0, eSpectdNdE,kl
                                            ,"Ks",0, eSpectdNdE,ks
                                            ,"K+/-",0, eSpectdNdE,kch
                                          );
exit(0);
                                           
 M=100;
 int pdg=5, out=0;
 double sp0[NZ],sp1[NZ], sp2[NZ],sp3[NZ];
 SpectraFlag=0; basicSpectra(M,pdg,out,sp0);
 SpectraFlag=1; basicSpectra(M,pdg,out,sp1);
 SpectraFlag=2; basicSpectra(M,pdg,out,sp2);
 SpectraFlag=3; basicSpectra(M,pdg,out,sp3);
 displayPlot("spectra", "E", 1E-3*M,M,0,4,"Pythia6",0,eSpectdNdE,sp0 
                                        ,"Ausrtalian",0,eSpectdNdE,sp1
                                        ,"PPPC",0,eSpectdNdE,sp2
                                        ,"CosmiXs",0,eSpectdNdE,sp3);
                                          
}                                           

