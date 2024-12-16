#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

double electronFSR(double E, void*Mvoid)
{
   double M=*(double*)(Mvoid);
   double x=E/M;
   double me=5.11E-4;
   return 1./137./M_PI*(1+(1-x)*(1-x))*(log( 4*M*M*(1-x)/me/me)-1);
}

   char * name[4]={"Pythia_6", "Pythia_8", "PPPC", "CosmiXs"};
   char * out[7]={"A","e","p","ne","nm","ml","D"};  

double energyTest(double M, char*pdgTxt)
{  double Etot=0;
   double sp[NZ];
   int pdg; 
   sscanf(pdgTxt,"%d",&pdg);   
   for(int Nout=0;Nout<7;Nout++) 
   { basicSpectra(M,pdg,Nout,sp);
     int err=0;
     double e=simpson_arg(eSpectdNdE,sp,M*1E-6,M,1E-3,&err);
     if(Nout==2) e+=0.939*simpson_arg(SpectdNdE,sp,M*1E-6,M,1E-3,&err);
     if(Nout==6) e+=2*0.939*simpson_arg(SpectdNdE,sp,M*1E-6,M,1E-3,&err);
     if(err)
     {  char txt[200]; 
        sprintf(txt,"Problem in integration err=%d: SpectraFlag=%s M=%.2E pdg=%d Nout=%s",err, name[SpectraFlag],M, pdg,out[Nout]);
        printf("%s\n",txt);
//        displayPlot(txt,"E",M*1E-6,M,1,1,"EdNdE",0,eSpectdNdE, sp);  exit(1);
     }  
     if(Nout) Etot+=2*e; else Etot=e;
   }
   return Etot/(2*M)-1;
}                                         

double eTestPhyth6(double M)  { SpectraFlag=0; return energyTest(M,"11");}
double eTestPhyth8(double M)  { SpectraFlag=1; return energyTest(M,"11");} 
double eTestPPPC(double M)    { SpectraFlag=2; return energyTest(M,"11");} 
double eTestCosmiXs(double M) { SpectraFlag=3; return energyTest(M,"11");}

int main(int argc,char** argv)
{ 
   double M=100;      //  DM mass
   int pdg=5,         //  b-quark  
   Nout=0;            //  photons
/*
M=2;
pdg=130 ;
Nout=1;
*/
   char title[200];   
   double sp[7][NZ];  //  foe spectra

   for(SpectraFlag=0;SpectraFlag<4;SpectraFlag++) basicSpectra(M,pdg,Nout,sp[SpectraFlag]); 
   sprintf(title,"EdNdE for  %d,%d -> %s. M=%.3E",pdg,-pdg,out[Nout],M);  
   displayPlot(title, "E", 1.E-9*M,M,1,4 
//                                     ,"Theory", 0,electronFSR,&M
                                       ,name[0],0,eSpectdNdE,sp[0] 
                                       ,name[1],0,eSpectdNdE,sp[1]
                                       ,name[2],0,eSpectdNdE,sp[2]
                                       ,name[3],0,eSpectdNdE,sp[3]
                                       );    
                                                           
   spectUncert=1;
   spectraUncertainty(M, pdg, Nout, sp[4]);
   sprintf(title,"Uncertainty  EdNdE for  %d,%d -> %s. M=%.3E",pdg,-pdg,out[Nout],M);    
   displayPlot(title, "E", 1.E-9*M,M,1,1,"delta(EdNdE)",0,eSpectdNdE,sp[4]); 

   spectUncert=0;
                                                      
   SpectraFlag=1;
   sprintf(title,"Energy conservation E/(2Mdm)-1 for %s spectra",name[SpectraFlag]);
   displayPlot(title, "Mdm",5,1E5,1,8
                                ,"e",0,energyTest,"11"
                                ,"tay",0,energyTest,"13"
                                ,"W",0,energyTest,"24"
                                ,"h",0,energyTest,"25"
                                ,"g",0,energyTest,"21"
                                ,"s",0,energyTest,"3"
                                ,"b",0,energyTest,"5"
                                ,"t",0,energyTest,"6"  
      );

}                                           

