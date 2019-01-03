/**
 * resolved mass spring model
 * using the leap frog integrator. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "spring.h"


int NS; 
struct spring* springs;
void reb_springs();// to pass springs to display

double gamma_all; // for gamma  of all springs
double t_damp;    // end faster damping, relaxation
double t_print;   // for table printout 
char froot[30];   // output files
int npert=0;

int icentral=-1; // central mass location

#define NPMAX 10  // maximum number of point masses
double itaua[NPMAX],itaue[NPMAX]; // inverse of migration timescales
double itmig[NPMAX];  // inverse timescale to get rid of migration

void heartbeat(struct reb_simulation* const r);

void additional_forces(struct reb_simulation* r){
   spring_forces(r); // spring forces
}


int main(int argc, char* argv[]){
	struct reb_simulation* const r = reb_create_simulation();
        struct spring spring_mush; // spring parameters for mush
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_BASIC;
	r->boundary	= REB_BOUNDARY_NONE;
	r->G 		= 1.00;	//gravitational constant, should be 1 
        r->additional_forces = additional_forces;  // setup callback function for additional forces
        double mball = 1.0;          // total mass of ball
        double rball = 1.0;          // radius of a ball
        double tmax = 0.0;  // if 0 integrate forever

// things to set!  can be read in with parameter file
        double dt; 
        double b_distance,omegax,omegay,omegaz,ks,mush_fac,gamma_fac;
        double ratio1,ratio2,obliquity_deg;
        int lattice_type;
        double rad[NPMAX],mp[NPMAX];
        double aa[NPMAX],ee[NPMAX],ii[NPMAX];
        double longnode[NPMAX],argperi[NPMAX],meananom[NPMAX];
        int npointmass=0;

    if (argc ==1){
        strcpy(froot,"t1");   // to make output files
	dt	   = 1e-3;    // Timestep
        tmax       = 0.0;     // max integration time
        t_print    = 1.0;     // printouts for table
	lattice_type  = 0;    // 0=rand 1=hcp
        b_distance = 0.15;    // min separation between particles
        mush_fac    = 2.3;    // ratio of smallest spring distance to minimum interparticle dist
        ks          = 8e-2;   // spring constant
        // spring damping
        gamma_all   = 1.0;    // final damping coeff
        gamma_fac   = 5.0;    // factor initial gamma is higher that gamma_all
        t_damp      = 1.0;    // gamma from initial gamma 
                              // to gamma_all for all springs at this time
        ratio1      = 0.7;    // axis ratio resolved body  y/x  b/a
        ratio2      = 0.5;    // axis ratio c/a
        omegaz      = 0.2;    // initial spin
        omegax      = 0.0;    // initial spin
        omegay      = 0.0;    // initial spin
        obliquity_deg = 0.0;  // obliquity

        npointmass=1;  // add one point mass
        int ip=0;   // index
        mp[ip]       = 1.0;   // mass
        rad[ip]      = 0.0;   // display radius
        itaua[ip]    = 0.0;   // inverse drift rate in a
        itaue[ip]    = 0.0;   // inverse drift rate in e
        itmig[ip]    = 0.0;   // get rid of drift rate in inverse of this time
	// orbit
        aa[ip]       =7.0;    // distance of m1 from resolved body, semi-major orbital
	ee[ip]       =0.0;    // initial eccentricity
	ii[ip]       =0.0;    // initial inclination 
	argperi[ip]  =0.0;    // initial orbtal elements
	longnode[ip] =0.0;    // initial 
	meananom[ip] =0.0;    // initial 

     }
     else{
        FILE *fpi;
        fpi = fopen(argv[1],"r");
        char line[300];
        fgets(line,300,fpi);  sscanf(line,"%s",froot);   // fileroot for outputs
        fgets(line,300,fpi);  sscanf(line,"%lf",&dt);    // timestep
        fgets(line,300,fpi);  sscanf(line,"%lf",&tmax);  // integrate to this time
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_print); // output timestep
        fgets(line,300,fpi);  sscanf(line,"%d",&lattice_type); // particle lattice type
        fgets(line,300,fpi);  sscanf(line,"%lf",&b_distance); // min interparticle distance
        fgets(line,300,fpi);  sscanf(line,"%lf",&mush_fac);   // sets max spring length
        fgets(line,300,fpi);  sscanf(line,"%lf",&ks);         // spring constant
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_fac);  // factor initial gamma is higher
        fgets(line,300,fpi);  sscanf(line,"%lf",&gamma_all);  // damping final
        fgets(line,300,fpi);  sscanf(line,"%lf",&t_damp);     // time to switch
        fgets(line,300,fpi);  sscanf(line,"%lf",&ratio1);     // axis ratio for body b/a
        fgets(line,300,fpi);  sscanf(line,"%lf",&ratio2);     // second axis ratio   c/a
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegax);     // initial body spin
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegay);     // initial body spin
        fgets(line,300,fpi);  sscanf(line,"%lf",&omegaz);     // initial body spin
        fgets(line,300,fpi);  sscanf(line,"%lf",&obliquity_deg); // obliquity degrees
        fgets(line,300,fpi);  sscanf(line,"%d",&npointmass); // number of point masses
        for (int ip=0;ip<npointmass;ip++){
           fgets(line,300,fpi);  sscanf(line,"%lf %lf %lf %lf %lf",
             mp+ip,rad+ip,itaua+ip,itaue+ip,itmig+ip); 
           fgets(line,300,fpi);  sscanf(line,"%lf %lf %lf %lf %lf %lf",
             aa+ip,ee+ip,ii+ip,longnode+ip,argperi+ip,meananom+ip);
        }

     }
     // double obliquity = obliquity_deg*M_PI/180.0;   // in radians
     npert = 0;

/// end parameteres things to set /////////////////////////

        r->dt=dt; // set integration timestep
	const double boxsize = 3.8*rball;    // display
	reb_configure_box(r,boxsize,1,1,1);
	// r->softening      = b_distance/100.0;	// Gravitational softening length
	r->softening      = 0.0;


   // properties of springs
   spring_mush.gamma     = gamma_fac*gamma_all; // initial damping coefficient
   spring_mush.ks        = ks; // spring constant
   // spring_mush.smax      = 1e6; // not used currently
   double mush_distance=b_distance*mush_fac; 
       // distance for connecting and reconnecting springs

   
   FILE *fpr;
   char fname[200];
   sprintf(fname,"%s_run.txt",froot);  // store parameters for simulation
   fpr = fopen(fname,"w");

   NS=0; // start with no springs 

// for volume to be the same, adjusting here!!!!
   double volume_ratio = pow(rball,3.0)*ratio1*ratio2;  // neglecting 4pi/3 factor
   double vol_radius = pow(volume_ratio,1.0/3.0);

   rball /= vol_radius; // volume radius used to compute body semi-major axis
           // assuming that semi-major axis is rball
   // fprintf(fpr,"vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3 
   fprintf(fpr,"a %.3f\n",rball); 
   fprintf(fpr,"b %.3f\n",rball*ratio1); 
   fprintf(fpr,"c %.3f\n",rball*ratio2); 
   volume_ratio = pow(rball,3.0)*ratio1*ratio2;  // neglecting 4pi/3 factor
   fprintf(fpr,"vol_ratio %.6f\n",volume_ratio); // with respect to 4pi/3 
   // so I can check that it is set to 1
 
   // if ratio1 ==  1.0 then oblate
   // both ratios should be <1
   // if ratio1 == ratio2 then prolate

   // create resolved body particle distribution
// the body is orientated with principal axes along xyz axes
   if (lattice_type==0){
      // rand_football_from_sphere(r,b_distance,rball,rball*ratio1, rball*ratio2,mball );
      rand_football(r,b_distance,rball,rball*ratio1, rball*ratio2,mball );
   }
   if (lattice_type ==1){
      fill_hcp(r, b_distance, rball , rball*ratio1, rball*ratio2, mball);
   }
   if (lattice_type ==2){
      fill_cubic(r, b_distance, rball , rball*ratio1, rball*ratio2, mball);
   }
// for oblate z is up is narrow symmetry axis, xy  in plane is circular
// for prolate x is long symmetry axis in the plane, z is still up

   int il=0;  // index range for resolved body
   int ih=r->N;

   subtractcom(r,il,ih);  // move reference frame to resolved body 
   subtractcov(r,il,ih); // center of velocity subtracted 
   // spin it
   // this routine spins the body with Omega three components
   spin(r,il, ih, omegax, omegay, omegaz); 
       // can spin about non principal axis
   subtractcov(r,il,ih); // center of velocity subtracted 

   // rotate_body(r, il, ih, 0.0, -atan2(llz,lly),0); 
   // rotate by Euler angles, including velocities

   // make springs, all pairs connected within interparticle distance mush_distance
   connect_springs_dist(r,mush_distance, 0, r->N, spring_mush); 

   // assume minor semi is rball*ratio2
   double ddr = rball*ratio2 - 0.5*mush_distance;
   ddr = 0.4; // mini radius  for computing Young modulus
   double Emush = Young_mush(r,il,ih, 0.0, ddr); // compute from springs out to ddr
   double Emush_big = Young_mush_big(r,il,ih);
   printf("ddr = %.3f mush_distance =%.3f \n",ddr,mush_distance);
   printf("Young's modulus %.6f\n",Emush);
   printf("Young's modulus big %.6f\n",Emush_big);
   fprintf(fpr,"Young's_modulus %.6f\n",Emush);
   fprintf(fpr,"Young's_modulus big %.6f\n",Emush_big);
   fprintf(fpr,"mush_distance %.4f\n",mush_distance);
   double LL = mean_L(r);  // mean spring length
   printf("mean spring length = %.4f\n",LL);
   fprintf(fpr,"mean_L (mean spring length)  %.4f\n",LL);
   // relaxation timescale
   // note no 2.5 here!
   double tau_relax = 1.0*gamma_all*0.5*(mball/(r->N -npert))/spring_mush.ks; // Kelvin Voigt relaxation time
// factor of 0.5 is consistent with damping force law use of reduced mass 
// this is equation 31 by Frouard+16
   printf("relaxation time %.3e\n",tau_relax);
   fprintf(fpr,"relaxation_time  %.3e\n",tau_relax);
   double viscosity = tau_relax*Emush_big/2.5;
   printf("viscosity   %.3e\n",viscosity);
   fprintf(fpr,"viscosity   %.3e\n",viscosity);
  
   //compute moment of inertia matrix
   double Ixx,Iyy,Izz,Ixy,Iyz,Ixz;
   mom_inertia(r,0,r->N, &Ixx, &Iyy, &Izz,&Ixy, &Iyz, &Ixz);
   // get its eigenvalues eig1>eig2>eig3
   double eig1,eig2,eig3;
   eigenvalues(Ixx,Iyy,Izz,Ixy,Iyz,Ixz, &eig1, &eig2, &eig3);

   // compute angular momentum vector body with respect to its center
   // of mass position and velocity
   double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz);
   double Lang = sqrt(llx*llx + lly*lly + llz*llz);

   // compute angle between angular momentum and principal axis and precession rate 
   double theta = 0.0;
   double omega_3 = 0.0;
   double omega_prec = 0.0; // precession rate 
   if (ratio1 ==1){ // oblate
      double hsquare = ratio2*ratio2; // is (c/a)^2
      theta = acos(llz/Lang); // axis of symmetry is z axis
      omega_3 = Lang/eig1; // eig1 = I3 is largest moment of inertia, corresponding to smallest body axis
      omega_prec = (1.0 - hsquare)/(1.0 + hsquare);
      omega_prec *= omega_3*cos(theta);
      printf("oblate\n");
   }
   if (ratio1 ==ratio2){ // prolate
      double hsquare = 1.0/(ratio2*ratio2); // is a/c
      theta = acos(llx/Lang); // axis of symmetry is x axis
      omega_3 = Lang/eig3; // eig3 is smallest moment,  is I parallel, largest body axis
      omega_prec = (1.0 - hsquare)/(1.0 + hsquare);
      omega_prec *= omega_3*cos(theta);
      printf("prolate\n");
   }
   printf("theta (deg)= %.6f (radians)= %.6f\n",theta*180.0/M_PI,theta);
   fprintf(fpr,"theta = %.6f radians, %.6f (degrees)\n",theta,theta*180.0/M_PI);
   printf("omega_prec= %.3e\n",omega_prec);
   fprintf(fpr,"omega_prec= %.3e\n",omega_prec);
   printf("omega_3= %.6f\n",omega_3);
   fprintf(fpr,"omega_3= %.6f\n",omega_3);
 
// do rotations afterwards!
// -------------------------------
  // I would like to orient the spinning body so that angular momentum is up
  // Euler angles: alpha,beta,gamma
  //   first angle alpha: about z in xy plane
  //   second angle beta: about x' axis in yz plane
  //   lastly gamma: rotate about z'' axis in xy plane
  // in our viewer up is y direction


   // compute angular momentum vector body with respect to its center
   // of mass position and velocity
   // double llx,lly,llz;
   measure_L(r,il, ih, &llx, &lly, &llz);
   // printf("llx lly llz = %.2f %.2f %.2f\n",llx,lly,llz);
   rotate_body(r, il, ih, 0.0, -atan2(llz,lly),0); // rotate by Euler angles, including velocities
   measure_L(r,il, ih, &llx, &lly, &llz);
   // printf("llx lly llz = %.3f %.3f %.3f\n",llx,lly,llz);
   rotate_body(r, il, ih,-atan2(lly,llx),0,0); 
   measure_L(r,il, ih, &llx, &lly, &llz);
   // printf("llx lly llz = %.3f %.3f %.3f\n",llx,lly,llz);
   rotate_body(r, il, ih,0,0,M_PI/2.0); 
   measure_L(r,il, ih, &llx, &lly, &llz);
   printf("llx lly llz = %.3f %.3f %.3f\n",llx,lly,llz);
   // ---------------- body should now be oriented with L up (in y direction)


   double barchi = 1.0*fabs(omega_prec)*tau_relax;  // initial value of barchi
   fprintf(fpr,"barchi  %.6e\n",barchi);
   printf("barchi %.6e\n",barchi);

   // ratio of numbers of particles to numbers of springs for resolved body
   double Nratio = (double)NS/(ih-il);
   printf("N=%d  NS=%d NS/N=%.1f\n", r->N, NS, Nratio);
   fprintf(fpr,"N=%d  NS=%d NS/N=%.1f\n", r->N, NS,Nratio);
   fclose(fpr);
   
   // r->particles[0].x += 1.1; // kick one particle
   reb_springs(r); // pass spring index list to display
   r->heartbeat = heartbeat;
   subtractcov(r,il,ih); // center of velocity subtracted 
   centerbody(r,il,ih);  // move reference frame to resolved body 


   // max integration time
   if (tmax ==0.0)
      reb_integrate(r, INFINITY);
   else
      reb_integrate(r, tmax);
}


#define NSPACE 50
void heartbeat(struct reb_simulation* const r){
        static int first=0;
        static char extendedfile[50];
        // static char pointmassfile[NPMAX*NSPACE];
        if (first==0){
           first=1;
           sprintf(extendedfile,"%s_ext.txt",froot);
        }
	if (reb_output_check(r,10.0*r->dt)){
		reb_output_timing(r,0);
	}
        if (fabs(r->t - t_damp) < 0.9*r->dt) set_gamma(gamma_all); 
            // damp initial bounce only 
            // reset gamma only at t near t_damp
	
         // stuff to do every timestep
         centerbody(r,0,r->N);  // move reference frame to resolved body for display


	 if (reb_output_check(r,t_print)) {
            print_extended(r,0,r->N,extendedfile); // orbital info and stuff
         }


}

// make a spring index list for display
void reb_springs(struct reb_simulation* const r){
   r->NS = NS;
   r->springs_ii = malloc(NS*sizeof(int));
   r->springs_jj = malloc(NS*sizeof(int));
   for(int i=0;i<NS;i++){
     r->springs_ii[i] = springs[i].i;
     r->springs_jj[i] = springs[i].j;
   }
}


