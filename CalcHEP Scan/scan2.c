#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
#include <dlfcn.h>
#include <sys/wait.h> 
#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h" 
#include <time.h>

int main(void)
{ int err,i;

	/* INTPUT PARAMETERS (to scan over) */

	/* OUTPUT PARAMETERS */
    // Higgs decay branching ratios
    double  wh,braa;
txtList branchings_MP,branchings_HD;

//set model dir here
char mdldir[] = "Work/models";

 // Set model number and number of points to collect, mdlnr is your model number
int mdlnr=6;

//a model to switch between to reset values when reloading
 setModel(mdldir , mdlnr ); 

/*****************************************************************************/
 srand (time(NULL)); //this is used to seed the random number by the system time

 if (remove("scan2.dat") == -1)
	perror("Error in deleting a file");

 FILE *file;
 file = fopen("scan2.dat","a+"); /* apend file (add text to
					a file or create a file if it does not exist.*/

 // Writing parameter names at first line to keep track of columns:
 //input parameters (1)
 //output parameters (3)
 fprintf(file,"MHD\t\t Br(MP->V' mu) \t Br(Mp->HD mu)\t  Br(HD->Mu mu) \n");		
 fclose(file); /*done with header of file*/

 /*** Starting randomizing loop ***/
int npoints=50;
double MD1_min = 10, MD1_max = 1e3;
double step_size = 1;

/*double wMP,wHD,BrMP__Vp_mu,BrMP__HD_mu,BrHD__mu_mu;
double MMP=1100,MMD=1050,MV=10;
double MD1_min = 10, MD1_max = 1e3;*/
 for (i = 0; i <= npoints; i++){

 /********** generate random values for variables **********/
 /*Mh     = Mhmin+(double) random()/RAND_MAX*(Mhmax-Mhmin);*/
 double MHD,LogMHDmin=0,LogMHDmax=4;
 
 MHD  = pow(10,(LogMHDmin+ i*(LogMHDmax-LogMHDmin)/npoints));

 /* Have to reset model every time, otherwise widths are not recalculated */
 setModel(mdldir , mdlnr ); 

 /******* assign variable values ********/
 /* the string is the calchep var name */
 /*err=assignValW("MMD", MMD);
  err=assignValW("MMP", MMP);
  err=assignValW("MHD", MHD);
  err=assignValW("MV", MV); (example)*/

  err=assignValW('%Mh1', MD1);
  err=assignValW('%Mh2', MD2);
  err=assignValW('%Mhc', MDP);

 // Calculation of public constraints  
 err=calcMainFunc();

 if(err!=0) { 
	  printf("Can not calculate constrained parameter %s\n",varNames[err]);i--;
 }
 else {
		// if the point survives the constraints collect more output values:
		// width and branchings of a particle
		wMP    = pWidth("mp",&branchings_MP);
		wHD    = pWidth("HD",&branchings_HD);
		
		BrMP__Vp_mu= findBr(branchings_MP,"V0,m");
		BrMP__HD_mu= findBr(branchings_MP,"HD,m");
		BrHD__mu_mu= findBr(branchings_HD,"m,M");
	
		// write values to file
  		file  = fopen("scan2.dat","a+");
		//input parameters
  		fprintf(file,"%f\t",MHD);
		//output parameters
  		fprintf(file,"%e\t%e\t%e\n",BrMP__Vp_mu, BrMP__HD_mu,BrHD__mu_mu);
  		fclose(file); 
 }
  
 }// *** end of rand loop ***

  return 0;
}

