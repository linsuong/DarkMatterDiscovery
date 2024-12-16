#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
using namespace LHAPDF;
using namespace std;

#include <cstring>


extern "C" {

#include"include/LHA2C.h"

void * initpdfLHA(char * name, int num)
{  setVerbosity(0);
   PDF*pdf=mkPDF(name,num);
   if(pdf==NULL) return NULL;
   return pdf;
}   

void pdfInfoLHA(void *p,  int*ID, int*maxNum, double*xMin,double*xMax,double*qMin,double *qMax, int *pdg)
{  PDF*pdf =(PDF*)p;
   if(maxNum) *maxNum= pdf->info().get_entry_as<int>("NumMembers");
   if(ID) *ID=pdf->lhapdfID();
   if(xMin) *xMin=pdf->info().get_entry_as<double>("XMin");
   if(xMax) *xMax=pdf->info().get_entry_as<double>("XMax");
   if(qMin) *qMin=pdf->info().get_entry_as<double>("QMin");
   if(qMax) *qMax=pdf->info().get_entry_as<double>("QMax"); 
   if(pdg)
   {  vector<int> pids = pdf->flavors();
      int i=0;
      for (int pid : pids) pdg[i++]=pid;
      pdg[i]=0;
   }       
}

int pdgInLHA( void *p,int pdg)
{    
   PDF*pdf =(PDF*)p;
   vector<int> pids = pdf->flavors();
   for (int pid : pids) if(pid==pdg) return 1;
   return 0;
}

double pdfFuncLHA(void*p, int pdg, double x, double q) { PDF*pdf=(PDF*)p; return pdf->xfxQ(pdg,x,q)/x;}

double alphaLHA(void*p,double q) { if(!p) return NAN;    PDF*pdf=(PDF*)p; return pdf->alphasQ(q);}

void deletePdfLHA(void *p) { if(p) { PDF*pdf =(PDF*)p; delete(pdf); } }


char** listLHA(void)
{
   char**list=NULL;
   int n=0;
   for (const string& setname : LHAPDF::availablePDFSets()) 
   { list=(char**) realloc(list,(n+1)*sizeof(char*));
     list[n]=(char*)malloc(setname.length()+1);
      strcpy(list[n],setname.c_str());
      n++;
   } 
   list=(char**) realloc(list,(n+1)*sizeof(char*));
   list[n]=NULL;
   return list;
}

}

