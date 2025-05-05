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
txtList branchings_Z,branchings_HD;

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
 fprintf(file,"MD1 \t MD2\t DMP \t DM3 \t Br(h2->e+e-h1) \t Br(h2->mu+mu-h1) \t Br(h2->tau+tau-h1) \t Br(h2->n+n-h1) \t Br(h2->Zh1) \t Br(Z->e+e-) \t Br(Z->mu+mu-) \t Br(Z->tau+tau-) \t Br(Z->n+n-)\n");		
 fclose(file); /*done with header of file*/

 /*** Starting fixed values loop ***/

double MD1_values[] = {1, 10, 20, 30, 40, 50, 60, 70, 80};
double DMP_values[] = {5, 20, 40, 60, 80, 120};
double DM3_values[] = {1, 10, 100};

int n_MD1 = sizeof(MD1_values) / sizeof(MD1_values[0]);
int n_DMP = sizeof(DMP_values) / sizeof(DMP_values[0]);
int n_DM3 = sizeof(DM3_values) / sizeof(DM3_values[0]);

for (int i_MD1 = 0; i_MD1 < n_MD1; i_MD1++) {
    for (int i_DMP = 0; i_DMP < n_DMP; i_DMP++) {
        for (int i_DM3 = 0; i_DM3 < n_DM3; i_DM3++) {
            double MD1 = MD1_values[i_MD1];
            double DMP = DMP_values[i_DMP];
            double DM3 = DM3_values[i_DM3];
            double MD2 = DM3 + DMP + MD1;

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
                printf("Can not calculate constrained parameter %s\n",varNames[err]);
            }
            else {
                // if the point survives the constraints collect more output values:
                // width and branchings of a particle
                double wZ, wHD;
                double BrHD__electron, BrHD__muon, BrHD__tau, BrHD__neutrino, BrHD__Z;
                double BrZ__electron, BrZ__muon, BrZ__tau, BrZ__neutrino;
                wHD = pWidth("~h2",&branchings_HD);
                wZ = pWidth("Z", &branchings_Z);

                
                BrHD__electron= findBr(branchings_HD,"e1,E1,~h1");   // h2 -> e+ e- ~h1
                BrHD__muon= findBr(branchings_HD,"e2,E2,~h1"); // h2 -> mu+ mu- ~h1
                BrHD__tau= findBr(branchings_HD,"e3,E3,~h1"); // h2 -> tau+ tau- ~h1
                BrHD__neutrino = findBr(branchings_HD, "n1, N1, ~h1"); 

                BrHD__Z = findBr(branchings_HD, "Z, ~h1");

                BrZ__electron = findBr(branchings_Z,"e1,E1");
                BrZ__muon = findBr(branchings_Z,"e2,E2" );
                BrZ__tau = findBr(branchings_Z, "e3,E3");
                BrZ__neutrino = findBr(branchings_Z, "n1, N1");

                
                //cSec = findValW()
                // write values to file
                file  = fopen("scan2.dat","a+");
                //fprintf(file,"MD1 \t MD2\t DMP \t DM3 \t Br(h2->e+e-h1) \t Br(h2->mu+mu-h1) 
                //	\t Br(h2->tau+tau-h1) \t Br(h2->n+n-h1) \t Br(h2->Zh1) 
                //	\t Br(Z->e+e-h1) \t Br(Z->mu+mu-h1) \t Br(Z->tau+tau-h1) \t Br(Z->n+n-h1)\n");		

                fprintf(file,"%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n", MD1, MD2, DMP, DM3,
                    BrHD__electron, BrHD__muon, BrHD__tau, BrHD__neutrino, BrHD__Z,
                    BrZ__electron, BrZ__muon, BrZ__tau, BrZ__neutrino
                     );
                fclose(file);
            }
        }
    }
}

  return 0;
}

//zeros cause leading branching is 