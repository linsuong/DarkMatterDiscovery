//Computation of the Sommerfeld enhancement factor for both the Yukawa potential and the hulthen potential
//These assume fermionic dark matter with a scalar or vector mediator

/*====== Modules ===============
   Keys to switch on
   various modules of micrOMEGAs
================================*/

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"



typedef struct { double mu;    //  mass of mediator divided by reduced mass of colliding particles.
                 double alpha; //  Yukawa coupling  g1*g2/(4 pi)
                 double vr;    //  relative velocity
               } SommPar;


static double SY_v(double v, void * par_)  // Sommerfeld factor for Yukawa potential as a funtion of velocity
{
  SommPar *par=par_;
  double a=par->alpha/v, b=par->mu/v;
  return SYukawa(a, b);
}


double SY_mu(double mu, void * par_)      // Sommerfeld factor for Yukawa potential  as a funtion of mu
{
  SommPar *par=par_;
  double a=par->alpha/par->vr, b=mu/par->vr;
  return SYukawa(a, b);
}


double SH_v(double v, void*par_)   // Sommerfeld factor for Hulthen  potential as a function of v
{   SommPar *par=par_;
    double a=par->alpha/v, b=par->mu/v;
    return SHulthen(a,b);
}



double SH_mu(double mu, void*par_)   // Sommerfeld factor for Hulthen  potential as a function of mu
{   SommPar *par=par_;
    double a=par->alpha/par->vr, b=mu/par->vr;
    return SHulthen(a,b);
}

double SC_v(double v,void*par_)          // Sommerfeld factor for Coulomb potential   as a funtion of v
{  SommPar *par=par_;
   double a=2*par->alpha/v;
   return M_PI*a/(1-exp(-M_PI*a));
}



int main(int argc,char** argv)
{  int err,n,i;
//                     mu   alpha  v
   SommPar  par[5]={ {0.01, 0.01, 1E-3}
                    ,{0.01, 0.10, 1E-3}
                    ,{0.10, 0.01, 1E-3}
                    ,{0.10, 0.10, 1E-3}
                    ,{0.01,-0.01, 1E-3}
                   };

   char chpar[5][20];
   for(int i=0;i<4;i++) sprintf(chpar[i],"mu=%.2f,al=%.2f",par[i].mu,par[i].alpha);

// For four representative values of mu and alpha, here we plot the Sommerfeld enhancement factor
// as function of the relative velocity for the Yukawa potential

   displayPlot(" Sommerfeld enhancement for Yukawa","v",1E-4,1, 1,4,chpar[0],0,SY_v,par+0
                                                                   ,chpar[1],0,SY_v,par+1
                                                                   ,chpar[2],0,SY_v,par+2
                                                                   ,chpar[3],0,SY_v,par+3
                                                       );

// For the same values of mu and alpha as above, here we compare the Sommerfeld enhancement factor
//  for the Yukawa (same color code as above) and Hulthen potentials
   displayPlot(" Yukawa/Hulthen comparison","v",1E-4,1, 1,8
                                                          ,"Y",0,SY_v,par+0
                                                          ,"Y",0,SY_v,par+1
                                                          ,"Y",0,SY_v,par+2
                                                          ,"Y",0,SY_v,par+3
                                                           ,"H",0,SH_v,par+0
                                                           ,"H",0,SH_v,par+1
                                                           ,"H",0,SH_v,par+2
                                                           ,"H",0,SH_v,par+3   
                                                         );
  char mess[100];

  // Comparison of  the Sommerfeld enhancement factor
  //  for the Yukawa and Hulthen potentials for negative values of alpha

  sprintf(mess,"Yukawa/Hulthen for repulsion  mu=%.2f al=%.2f",par[4].mu,par[4].alpha);
  displayPlot(mess, "v",1E-4,1,1,2, "Yukawa",0,SY_v,par+4
                                  , "Hulthen",0,SH_v,par+4
                                  );

   SommPar  clmb={ 1E-4, 0.01,1E-3};

  //Comparison of the Sommerfeld enhancement factor for Yukawa, Hulthen and Coulomb potential, here mu=1e-4 and alpha=0.01
  //Here we assume a vector mediator

  sprintf(mess,"Yukawa/Hulthen/Coulomb mu=%.2E, alpha=%.2E",clmb.mu,clmb.alpha);
  displayPlot(mess,"v",1E-5,1,1,3,"Yukawa",0, SY_v,&clmb
                                 ,"Hulthen",0,SH_v,&clmb
                                 ,"Coulomb",0,SC_v,&clmb
                               );

   SommPar pole={0.1,0.01, 1E-3};  //  mu, alpha, vr

   sprintf(mess,"alpha=%.2E, v=%.2E  Resonance is expected at mu=%.2E",pole.alpha, pole.vr,

   (6*pole.alpha + sqrt( 36*pole.alpha*pole.alpha - 144*pole.vr*pole.vr))/M_PI/M_PI
   );

   displayPlot(mess,"mu", 0.00001,1, 1,2,"Yukawa",  0,  SY_mu,       &pole
                                        ,"Hulthen", 0, SH_mu,&pole);

   killPlots();
   return 0;
}
