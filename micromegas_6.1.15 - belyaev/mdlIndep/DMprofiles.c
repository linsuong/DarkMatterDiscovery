/*====== Modules ===============
     Display DM-profiles in Milky Way galaxy
================================*/


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"


int main(int argc,char** argv)
{ int err,n,i;


  setProfileZhao(1,3,1,8.1);
  displayPlot("1906.08419, rhoDM=0.51","R[kpc]", 0,50,0,1,"NFW",0,hProfileZhao,NULL);
   
  setProfileEinasto(2.33,10.7);
  displayPlot("2309.00048","R[kpc]", 0,50,0,1,"Einasto",0,hProfileEinasto,NULL);


  setProfileEinasto(0.91,6.65);
  setProfileZhao(1,3,0.0258,5.26);
  displayPlot("2303.12838","R[kpc]", 0,50,0,2,"Einasto",0,hProfileEinasto,NULL,"Zhao",0,hProfileZhao,NULL);

  killPlots();

}
