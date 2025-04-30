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

double rnd_slg(double xmin, double xmax) {
    //randomly scans with a log scale as well as a negative side
    double sgn=1.;
    if(xmin<0) { 
    xmin=1E-10;
    sgn= (rand() % 2 == 0) ? 1 : -1; // Randomly return +1 or -1 
    }
    
    double log_min = log(xmin);
    double log_max = log(xmax);
    
    // Generate a random value on the log scale
    double random_log = log_min + ((double)rand() / RAND_MAX) * (log_max - log_min);
   
    return sgn*exp(random_log);
      }

double rnd_lin(double xmin, double xmax) {
  //not using this!
    return xmin + (xmax - xmin) * ((double)rand() / RAND_MAX);
}

int main(void)
{ int err,i;

	/* INTPUT PARAMETERS (to scan over) */

	/* OUTPUT PARAMETERS */
    // Higgs decay branching ratios
    double  wh,braa;
txtList branchings_W,branchings_HD;

//set model dir here
char mdldir[] = "../models";

 // Set model number and number of points to collect, mdlnr is your model number
int mdlnr=6;

//a model to switch between to reset values when reloading
 setModel(mdldir , mdlnr); 

 
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
 fprintf(file,"MD1 \t MD2\t DMP \t DM2 \t Br(h-->e-nu) \t Br(h-->mu-nu) \t Br(h-->tau-nu) \t Br(h-->W-h1) \t Br(W-->e-nu) \tBr(W-->mu-nu) \t Br(W-->tau-nu)\n");		
 fclose(file); /*done with header of file*/

 /*** Starting randomizing loop ***/
 int npoints=50000;

double MD1_min = 10, MD1_max = 1e3;
double DMP_min = 10, DMP_max = 1e2;
double DM3_min = 10, DM3_max = 1e2;

double step_size = 1;

 for (i = 0; i <= npoints; i++){

 /********** generate random values for variables **********/
 //double MHD,LogMHDmin=10,LogMHDmax=100;
 double MD1, DMP, DM3, MD2;
 
 //MHD  = pow(10,(LogMHDmin+ i*(LogMHDmax-LogMHDmin)/npoints));
 MD1 = rnd_slg(MD1_min,MD1_max); //MD1 = mass of DM particle h1
 DMP = rnd_slg(DMP_min,DMP_max); //DMP = MDP-MD1                                      
 DM3 = rnd_slg(DM3_min,DM3_max); //DM3 = MD2-MDP    
 MD2 = DM3 + DMP + MD1;

 /* Have to reset model every time, otherwise widths are not recalculated */
 setModel(mdldir , mdlnr);

 /******* assign variable values ********/
 // Debugging: Print values before assignment
printf("MD1 = %f, Mh1+DMP = %f, Mhc+DM3 = %f\n", MD1, DMP, DM3);

// Assigning values to the model parameters
// string =  CalcHEP  vars.mdl expression, variable can be whatever
err = assignValW("MD1", MD1);
err = assignValW("DMP", DMP);
err = assignValW("DM3", DM3);
//err = assignValW('~h2', MD2);

 // Calculation of public constraints
 err=calcMainFunc();

 if(err!=0) {
	  printf("Can not calculate constrained parameter %s\n",varNames[err]);i--;
 }
 else {
		// if the point survives the constraints collect more output values:
		// width and branchings of a particle
		double widthW, wHD;
		double BrHD__electron, BrHD__muon, BrHD__tau, BrHD__neutrino, BrHD__W;
		double BrW__electron, BrW__muon, BrW__tau, BrW__neutrino;
		wHD    = pWidth("~h-",&branchings_HD);
		widthW = pWidth("W-", &branchings_W);
		
		BrHD__electron= findBr(branchings_HD,"N1,e1,~h1");
		BrHD__muon= findBr(branchings_HD,"N2,e2,~h1"); 
		BrHD__tau= findBr(branchings_HD,"N3,e3,~h1");

		BrHD__W = findBr(branchings_HD, "W-, ~h1");

		BrW__electron= findBr(branchings_W,"N1,e1");
		BrW__muon= findBr(branchings_W,"N2,e2"); 
		BrW__tau= findBr(branchings_W,"N3,e3");

		//cSec = findValW()
		// write values to file
  		file  = fopen("scan2.dat","a+");
  		fprintf(file,"%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", MD1, MD2, DMP, DM3, 
		BrHD__electron, BrHD__muon, BrHD__tau, BrHD__W, 
		BrW__electron,BrW__muon, BrW__tau);
  		fclose(file);	
 }

//findValW() to get value of something.
 }// *** end of rand loop ***

  return 0;
}
