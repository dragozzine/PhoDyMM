#define demcmc_compile 0
#define celerite_compile 1

// This code written in C by Sean M Mills

// To compile lcout:
// make sure demcmc_compile is defined as 0
// Compile with:
// g++ -w -O3 -o lcout -I/home/smills/celerite/celerite/cpp/include -I/home/smills/celerite/celerite/cpp/lib/eigen_3.3.3 -lm -lgsl -lgslcblas -fpermissive phodymm.cpp 
// or
// gcc -w -O3 -o lcout.c -lm -lgsl -lgslcblas -fpermissive phodymm.cpp
// Run with:
// ./lcout demcmc.in kep35.pldin [[-rv0=rvs0.file] [-rv1=rvs1.file] ... ]

// To compile demcmc:
// make sure demcmc_compile is defined as 1
// mpic++ -w -Ofast -o demcmc -I/home/smills/celerite/celerite/cpp/include -I/home/smills/celerite/celerite/cpp/lib/eigen_3.3.3 -lm -lgsl -lgslcblas -lmpi -fpermissive phodymm.cpp
// or 
// mpic++ -w -Ofast -o demcmc -lm -lgsl -lgslcblas -lmpi -fpermissive phodymm.cpp
// Run with:
// mpirun ./demcmc demcmc.in kep.pldin

// To compile longterm stability
// make sure demcmc_compile is defined as 3
// Compile with:
// $ gcc -Ofast -o stability -lgsl -lgslcblas -fopenmp phodymm.c
// g++ -w -O3 -o lcout -lm -lgsl -lgslcblas -fpermissive phodymm.cpp

#if (demcmc_compile==1) 
#include <mpi.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include <gsl/gsl_multimin.h>

#if (celerite_compile == 1)
// For CELERITE
#include <cmath>
#include <iostream>
#include <Eigen/Core> 
#include "celerite/celerite.h"
using Eigen::VectorXd;
#endif

int CELERITE;
int NCELERITE=4;
int RVCELERITE;
int NRVCELERITE=4;

// These variables are now all defined in the input file instead of here
// turn on N temperature annealing 
int NTEMPS;
// turn on one sided inclination distribution
int IGT90;
// turn on positive dilution
int INCPRIOR;
// Limb darkening (Pal+2011, or Maxsted+2018)
int LDLAW;
// dilution capped?
int DIGT0;
// sqrt(e) as parameter?
int SQRTE;
// restrict masses to > 0
int MGT0;
double MASSHIGH; 
int MASSPRIOR;
double *MASSPRIORCENTERS;
double *MASSPRIORSIGMAS;
// density cuts
int DENSITYCUTON;
double *MAXDENSITY;// g/cm^3
double MSUNGRAMS = 1.98855e33; // g
double RSUNCM = 6.95700e10; // cm
// eccentricity cuts. 
int ECUTON;
double *EMAX;
// turn on eccentricity prior.
int EPRIOR = 0;
int *EPRIORV;
// this sigma is the sigma parameter in a Rayleigh distribution
double ESIGMA;
// define rayleigh distribution
double rayleighpdf (double x) {
  double sigma = ESIGMA;
  if (x < 0.0) {
    printf(" rayleigh pdf requires positive x value\n");
    exit(0);
  }
  double f =  x / (sigma*sigma) * exp( -x*x / (2*sigma*sigma) );
  return f;
}
// define normal distribution
double normalpdf (double x) {
  double sigma = ESIGMA;
  double f = (1.0 / (sigma * sqrt(2.*M_PI))) *
        exp(-x*x / (2.0 * sigma * sigma));
  return f;
}


// Spectroscopy Constraints - Currently only works for single star. Don't use it for multistar. 
int SPECTROSCOPY;
// Assumes assymetric Gaussian
double SPECRADIUS=0.;
double SPECERRPOS=0.;
double SPECERRNEG=0.;

// Spectroscopy Constraints - Currently only works for single star. Don't use it for multistar. 
int MASSSPECTROSCOPY;
// Assumes assymetric Gaussian
double SPECMASS=0.;
double MASSSPECERRPOS=0.;
double MASSSPECERRNEG=0.;

int RANK;
int SIZE;

double PRINTEPOCH = 800.0;
//// Global Variables
// Initialized in main
int RESTART;
int RVS;
// Initialized in getinput function call
int MULTISTAR;
int NBODIES;
int PPERPLAN;
int PSTAR = 5;  // Number of parameters for central star
int PPERWALKER;
char *OUTSTR;
int NPL;
double EPOCH;
long NWALKERS;
long NSTEPS;
int CADENCESWITCH;
char TFILE[1000];
int *PARFIX;
double T0;
double T1;
unsigned long OUTPUTINTERVAL;
unsigned long SEED;
double DISPERSE;
double OPTIMAL;
double RELAX;
int NPERBIN;
double BINWIDTH;
double *PSTEP;
double *SSTEP;
double *STEP;
int STEPFLAG;
double *CSTEP;
int BIMODF;
int *BIMODLIST;
double OFFSETMULT;
double OFFSETMIN;
double OFFSETMINOUT;
double DIST2DIVISOR;
int LTE;
int SPLITINCO;
int XYZFLAG;
int *XYZLIST;
int *RVARR;
char **RVFARR;

int RVJITTERFLAG;
int NSTARSRV;
int NTELESCOPES;
int RVJITTERTOT = 0;

int TTVJITTERFLAG;
int TTVJITTERTOT = 0;
int NSTARSTTV;
int NTELESCOPESTTV;

int OOO = 2;
int CONVERT=1;

#if (demcmc_compile==3)
double TINTERVAL = 1000.0;
double AFRACSTABLE = 0.10;
#endif


// Variables for Fitting additional TTVs 
// as of now the ttv input files MUST be SORTED numerically and at same epoch
int TTVCHISQ;
long **NTTV;  //number
double **TTTV; // time
double **ETTV; //error
double **MTTV; //modeled


// Often used system constants
const int SOFD = sizeof(double);
const int SOFDS = sizeof(double*);
const int SOFI = sizeof(int);
const int SOFIS = sizeof(int*);


///* Integration parameters */
#define DY 1e-14                   ///* Error allowed in parameter values per timestep. */
#define HStart 1e-5                ///* Timestep to start.  If you get NaNs, try reducing this. */

///* Some physical constants */
#define G 2.9591220363e-4         ///*  Newton's constant, AU^3*days^-2 */
#define RSUNAU  0.0046491            ///* solar radius in au */
#define REARTHORSUN 0.009171    ///* earth radius divided by solar radius */
#define MPSTOAUPD 5.77548327363993e-7    ///* meters per second to au per day conversion factor */
#define MSOMJ 1.04737701464237e3  ///* solar mass in terms of jupiter masses */ 
#define CAUPD 173.144632674240  ///* speed of light in AU per day */

//Note in dpintegrator and related code:
///* y vector is [x,y,z,v_x,v_y,v_z] */
///* f is d/dt y */


// Check if two doubles are equal (at the DBL_EPSILON level)
int dbleq (double a, double b) {

  return fabs(a-b) < DBL_EPSILON;

}


// Find the greater of two values
int compare (const void* a, const void* b) {

  double dbl_a = * ( (double*) a);
  double dbl_b = * ( (double*) b);

  if (dbl_a < dbl_b) return -1;
  else if (dbl_b < dbl_a) return 1;
  else return 0;

}


// set a 6*npl element vector equal to another 6*npl element vector
int seteq (int npl, double ynew[], const double y[]) {

  int i;
  for(i=0;i<6*npl;i++) ynew[i]=y[i];

}


// N-body interaction between planets with position vector y and masses in params
int func (double t, const double y[], double f[], void *params) {

  int i,j;
  int i1,i2;
  double * npl_masses = (double *)params;
  double * masses = &npl_masses[1];
  double dnpl = npl_masses[0];
  const int npl = dnpl;
  double gmc1, gmc2, gm1, gm2;
  double rc1m3, rc2m3, r12m3;

  for(i=0; i<npl; i++) {
    gmc1 = G*(masses[0]+masses[i+1]);
    rc1m3 = pow(pow(y[i*6+0],2)+pow(y[i*6+1],2)+pow(y[i*6+2],2),-3.0/2);
    for(j=0; j<3; j++) {
      f[i*6+j] = y[i*6+3+j];  /* x dot = v */
      f[i*6+3+j] = -gmc1*y[i*6+j]*rc1m3;   /* Pull of the star. */
    }   
  }
 
  /* Interaction between each pair of planets. */
  /* Eqn 6.8,6.9 in Murray and Dermott */
  /* Astrocentric coordinates are used (i.e. center of star is origin) */
  for(i1=0; i1<npl-1; i1++) {
    gm1 = G*masses[i1+1];
    rc1m3 = pow(pow(y[i1*6+0],2)+pow(y[i1*6+1],2)+pow(y[i1*6+2],2),-3.0/2);
    for(i2=i1+1; i2<npl; i2++) {
      gm2 = G*masses[i2+1];
      rc2m3 = pow(pow(y[i2*6+0],2)+pow(y[i2*6+1],2)+pow(y[i2*6+2],2),-3.0/2);
      r12m3 = pow(pow(y[i1*6+0]-y[i2*6+0],2)+pow(y[i1*6+1]-y[i2*6+1],2)+pow(y[i1*6+2]-y[i2*6+2],2),-3.0/2);
  
      for(j=0; j<3; j++) f[i1*6+3+j] += -gm2*( (y[i1*6+j]-y[i2*6+j])*r12m3 + y[i2*6+j]*rc2m3 );
      for(j=0; j<3; j++) f[i2*6+3+j] += -gm1*( (y[i2*6+j]-y[i1*6+j])*r12m3 + y[i1*6+j]*rc1m3 );
    }
  }

  return GSL_SUCCESS;
}


// We don't use the Jacobian, but still need to define a pointer for gsl
void *jac;




// Helper for checking if certain parameters are out of 
int check_boundaries( double *p0local, long nw) {
  
  int notallowed = 0;
  const int npl = NPL;
  const int pperwalker = PPERWALKER;
  const int pperplan = PPERPLAN;
  const int sofd = SOFD;
  int ip;


 
  for (ip=0; ip<npl; ip++) { 
    // make sure i and Omega angles are not cycling through:
    if ( p0local[nw*pperwalker+ip*pperplan+4] < 0.0 || p0local[nw*pperwalker+ip*pperplan+4] > 180.0 ) notallowed=1;
    if ( p0local[nw*pperwalker+ip*pperplan+5] < -180.0 || p0local[nw*pperwalker+ip*pperplan+5] > 180.0 ) notallowed=1;
    // make sure i>=90 or i<= 90
    if ( (IGT90==1) && (p0local[nw*pperwalker+ip*pperplan+4] < 90.0) ) notallowed=1;
    if ( (IGT90==2) && (p0local[nw*pperwalker+ip*pperplan+4] > 90.0) ) notallowed=1;
    // make sure m>=0
    if ( MGT0 && (p0local[nw*pperwalker+ip*pperplan+6] < 0.0) ) notallowed=1;
    //make sure density is allowed
    if (DENSITYCUTON) {
      double massg = p0local[nw*pperwalker+ip*pperplan+6] / MSOMJ * MSUNGRAMS; 
      double radcm = p0local[nw*pperwalker+npl*pperplan+1] * p0local[nw*pperwalker+ip*pperplan+7] * RSUNCM;
      double rhogcc = massg / (4./3.*M_PI*radcm*radcm*radcm);
      if (rhogcc > MAXDENSITY[ip]) notallowed=1;
    }
    // make sure radius ratios are positive
    if (p0local[nw*pperwalker+ip*pperplan+7] < 0.0) notallowed=1;
  }

  if ( DIGT0 && (p0local[nw*pperwalker+npl*pperplan+4] < 0.0) ) notallowed=1;

  // check that RVJITTER >= 0
  if (RVS) {
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<RVJITTERTOT; ki++) {
        if (p0local[nw*pperwalker+npl*pperplan+(5+ki)] < 0.0) {
          notallowed=1;
        }
      }
    }
  }
  
  // check that celerite terms >= 0
  if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        if (p0local[nw*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+ki] < 0.0) {
          notallowed=1;
        }
      }
  }
  if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        if (p0local[nw*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+CELERITE*4+ki] < 0.0) {
          notallowed=1;
        }
      }
  }
  
  double* evector = malloc(npl*sofd);
  if (ECUTON) {
    if (SQRTE) {
      int i0;
      for (i0=0; i0<npl; i0++) {
        evector[i0] = pow(sqrt( pow(p0local[nw*pperwalker+i0*pperplan+2], 2) + pow(p0local[nw*pperwalker+i0*pperplan+3], 2) ), 2);
      }
    } else {
      int i0;
      for (i0=0; i0<npl; i0++) {
        evector[i0] = sqrt( pow(p0local[nw*pperwalker+i0*pperplan+2], 2) + pow(p0local[nw*pperwalker+i0*pperplan+3], 2) );
      }
    }
    // make sure e_d cos  w_d is within its range
    double* emax = EMAX;
    int i0;
    for (i0=0; i0<npl; i0++) {
      if (evector[i0] > emax[i0] ) notallowed=1;
    }
  }
  
  // make sure ld coefficients in range
  // this prior only applies to the Pal+2011 code
  if ( (LDLAW == 0) && (p0local[nw*pperwalker+npl*pperplan+2] < 0.0 || p0local[nw*pperwalker+npl*pperplan+2] > 1.0 || p0local[nw*pperwalker+npl*pperplan+3] < 0.0 || p0local[nw*pperwalker+npl*pperplan+3] > 1.0) ) notallowed=1;
 
  free(evector);
 
  return notallowed;
}



// Helper for computing the priors
double compute_priors(double *p0local, int i) {

    double neg2logliketemp = 0.;
    const int pperwalker = PPERWALKER;
    const int pperplan = PPERPLAN;
    const int sofd = SOFD;
    const int npl = NPL;

    double photoradius;
    photoradius = p0local[i*pperwalker+npl*pperplan+1]; 
    double photomass;
    photomass = p0local[i*pperwalker+npl*pperplan+0]; 
    double* evector; 
    evector = malloc(npl*sofd);

    if (EPRIOR) {
      if (SQRTE) {
        int i0;
        for (i0=0; i0<npl; i0++) {
          evector[i0] = pow(sqrt( pow(p0local[i*pperwalker+i0*pperplan+2], 2) + pow(p0local[i*pperwalker+i0*pperplan+3], 2) ), 2);
        }
      } else {
        int i0;
        for (i0=0; i0<npl; i0++) {
          evector[i0] = sqrt( pow(p0local[i*pperwalker+i0*pperplan+2], 2) + pow(p0local[i*pperwalker+i0*pperplan+3], 2) );
        }
      }
      int i0;
      for (i0=0; i0<NPL; i0++) {
        if (EPRIORV[i0]) {
          double priorprob;
          if (EPRIOR==1) {
            priorprob = rayleighpdf(evector[i0]);
          } else if (EPRIOR==2) {
            priorprob = normalpdf(evector[i0]);
          }
          neg2logliketemp += -2.0*log( priorprob );
        }
      }
    }

    if (MASSPRIOR) {
      int i0;
      for (i0=0; i0<NPL; i0++) {
        double massi = p0local[i*pperwalker+i0*pperplan+6] / MSOMJ;
        double massratioi = massi / photomass;
        neg2logliketemp += -2.*log(1./(sqrt(2.*M_PI)*MASSPRIORSIGMAS[i0]));
        neg2logliketemp += pow((massratioi - MASSPRIORCENTERS[i0])/MASSPRIORSIGMAS[i0], 2);
      } 
    }

    if (SPECTROSCOPY) {
      double normalization = 1./(sqrt(2.*M_PI));
      normalization *= 2./(SPECERRPOS+SPECERRNEG); // normalized asymmetric Gaussian
      neg2logliketemp += -2.*log(normalization);
      if (photoradius > SPECRADIUS) {
        neg2logliketemp += pow( (photoradius - SPECRADIUS) / SPECERRPOS, 2 );
      }
      else {
        neg2logliketemp += pow( (photoradius - SPECRADIUS) / SPECERRNEG, 2 );
      }
    }

    if (MASSSPECTROSCOPY) {
      double normalization = 1./(sqrt(2.*M_PI));
      normalization *= 2./(MASSSPECERRPOS+MASSSPECERRNEG); // normalized asymmetric Gaussian
      neg2logliketemp += -2.*log(normalization);
      if (photomass > SPECMASS) neg2logliketemp += pow( (photomass - SPECMASS) / MASSSPECERRPOS, 2 );
      else neg2logliketemp += pow( (photomass - SPECMASS) / MASSSPECERRNEG, 2 );
    }

    if (INCPRIOR) {
      int i0;
      for (i0=0; i0<npl; i0++) {
        neg2logliketemp += -2.0*log( sin(p0local[i*pperwalker+i0*pperplan+4] *M_PI/180.) ); 
      }
    }

  free(evector);

  return neg2logliketemp;

}



#if (celerite_compile == 1)
// Wrapper for computing a celerite fit and returning the effective chi^2
double celerite_fit(double*** flux_rvs, double* p0local, int i, int rvflag, int verbose)  { 

  const int npl = NPL;
  const int pperplan = PPERPLAN;
  const int pperwalker = PPERWALKER;
  const int sofd = SOFD;

  double neg2loglike;  

  double *xs = flux_rvs[rvflag][0];
  long maxil = (long) xs[0];
  double *trueys = flux_rvs[rvflag][1];
  double *modelys = flux_rvs[rvflag][2];
  double *es = flux_rvs[rvflag][3];
  double *diffys = malloc(sofd*maxil);
  int il;
  for (il=0; il<maxil; il++) { 
    diffys[il] = trueys[il+1]-modelys[il+1];
  }
  double *yvarp = malloc(sofd*maxil);
  for (il=0; il<maxil; il++) { 
    yvarp[il] = es[il+1]*es[il+1]; 
  }
  double *xp = &xs[1]; 
  
  int j_real = 0;
  int j_complex;
  double jitter, k1, k2, k3, S0, w0, Q;
  jitter = p0local[pperwalker*i+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE*rvflag+0];
  S0 = p0local[pperwalker*i+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE*rvflag+1];
  w0 = p0local[pperwalker*i+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE*rvflag+2];
  Q = p0local[pperwalker*i+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE*rvflag+3];
  if (Q >= 0.5) {
    j_complex = 1;
  } else {
    j_complex = 2;
  }
  VectorXd a_real(j_real),
          c_real(j_real),
          a_comp(j_complex),
          b_comp(j_complex),
          c_comp(j_complex),
          d_comp(j_complex);
  if (Q >= 0.5) {
    k1 = S0*w0*Q;
    k2 = sqrt(4.*Q*Q - 1.);
    k3 = w0/(2.*Q);
    a_comp << k1;
    b_comp << k1/k2;
    c_comp << k3;
    d_comp << k3*k2;
  } else {
    j_complex = 2;
    k1 = 0.5*S0*w0*Q;
    k2 = sqrt(1. - 4.*Q*Q);
    k3 = w0/(2.*Q);
    a_comp << k1*(1. + 1./k2), k1*(1. - 1./k2);
    b_comp << 0., 0.;
    c_comp << k3*(1. - k2), k3*(1. + k2);
    d_comp << 0., 0.;
  }
  
  if (verbose) {  
    printf("%lf %lf %lf %lf\n", jitter, S0, w0, Q);
  }

  VectorXd x = VectorXd::Map(xp, maxil);
  VectorXd yvar = VectorXd::Map(yvarp, maxil);
  VectorXd dy = VectorXd::Map(diffys, maxil);
  
  celerite::solver::CholeskySolver<double> solver;
  solver.compute(
        jitter,
        a_real, c_real,
        a_comp, b_comp, c_comp, d_comp,
        x, yvar  // Note: this is the measurement _variance_
    );
  
  // see l.186-192 in celerite.py
  double logdet, diffs, llike;
  logdet = solver.log_determinant();
  diffs = solver.dot_solve(dy); 
  llike = -0.5 * (diffs + logdet); 
  
  neg2loglike = diffs+logdet;
  
  if (verbose) { 
    printf("Celerite params:\n");
    printf("a+comp=%25.17lf\n", a_comp[0]);
    printf("b+comp=%25.17lf\n", b_comp[0]);
    printf("c+comp=%25.17lf\n", c_comp[0]);
    printf("d+comp=%25.17lf\n", d_comp[0]);
    printf("diffs=%lf\n", diffs);
    printf("logdet=%lf\n", logdet);
    printf("llike=%lf\n", llike);
  
    VectorXd prediction = solver.predict(dy, x);

    printf("predicted\n");
 
    char tmtefstr[1000];
    //memset(tmtefstr, '\0', sizeof(tmtefstr));
    char str1[100];
    if (rvflag) {
      strcpy(str1, "rv");
    } else {
      strcpy(str1, "lc");
    }
    strcpy(tmtefstr, str1);
    strcat(tmtefstr, "_");
    strcat(tmtefstr, OUTSTR);
    strcat(tmtefstr, ".gp_");
    strcat(tmtefstr, str1);
    strcat(tmtefstr, "out");
    printf("Saving GP fit to: %s\n", tmtefstr); 
    FILE *tmtef;
    tmtef = fopen(tmtefstr, "a");
    long ijk;
    for (ijk=0; ijk<maxil; ijk++) {
      fprintf(tmtef, "%.12lf \n", prediction[ijk]);
    }
    fclose(tmtef);// = openf(tmtstr,"w");
    
    printf("saved\n");
  }

  free(yvarp);
  free(diffys);

  return neg2loglike;

}
#endif 
#if (celerite_compile == 0)
// dummy function if celerite is not used to compile the code 
double celerite_fit(double*** flux_rvs, double* p0local, int i, int rvflag, int verbose) {} 
#endif
 



// Compute lightcurve if only 1 star / luminous object
double ***dpintegrator_single (double ***int_in, double **tfe, double **tve, double **nte, int *cadencelist) {

  const int cadenceswitch = CADENCESWITCH;
  double t0 = T0;
  double t1 = T1;
  int nperbin = NPERBIN;
  double binwidth = BINWIDTH;
  const int lte = LTE;

  const int sofi = SOFI;
  const int sofd = SOFD;
  const int sofds = SOFDS;

  long kk = (long) tfe[0][0];
  int nplplus1 = int_in[0][0][0];
  const int npl = nplplus1-1;
  double tstart = int_in[1][0][0]; //epoch
  double rstar = int_in[3][0][0];
  double c1 = int_in[3][0][1];
  double c2 = int_in[3][0][2];
  double dilute = int_in[3][0][3];


  const gsl_odeiv_step_type * T 
  /*   = gsl_odeiv_step_bsimp;   14 s */
  = gsl_odeiv_step_rk8pd;   /* 3 s */
  /*  = gsl_odeiv_step_rkf45;    14 s */
  /*  = gsl_odeiv_step_rk4;  26 s */

  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 6*npl);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new (DY, 0.0);
  gsl_odeiv_evolve * e 
    = gsl_odeiv_evolve_alloc (6*npl);
  
  long i;
  double mu[npl+1];
  double npl_mu[npl+2];
  gsl_odeiv_system sys = {func, jac, 6*npl, npl_mu};

  double t; /* this is the given epoch. */

  double hhere;  
  double y[6*npl], yin[6*npl];
  //note:rad is in units of rp/rstar
  double *rad = malloc(npl*sofd);

  mu[0] = int_in[2][0][0]; 
 
  for (i=0; i<npl; i++) {
    mu[i+1] = int_in[2][i+1][0];
    int ih;
    for (ih=0; ih<6; ih++) {
      yin[i*6+ih] = int_in[2][i+1][ih+1];
    }
    rad[i] = int_in[2][i+1][7];  
  }

  memcpy(&npl_mu[1], mu, (npl+1)*sofd);
  npl_mu[0] = npl;


  double mtot=mu[0];
  for(i=0;i<npl;i++) mtot += mu[i+1];

  double yhone[6*npl], dydt[6*npl], yerrhone[6*npl];
  double thone,tstep,dist,vtrans;
  double vstar,astar;


  int transitcount[npl];
  int pl;
  for (pl=0;pl<npl;pl++) transitcount[pl]=-1;

#if ( demcmc_compile==3)
  double msys[npl];
  msys[0]=mu[0]+mu[1];
  for (i=1; i<npl; i++) msys[i] = msys[i-1] + mu[i+1];

  double *statetokep(double x, double y, double z, double vx, double vy, double vz, double m); //prototype
  double *keptoorb(double x, double y, double z, double vx, double vy, double vz, double m); //prototype

  double *amax = malloc(npl*sofd);
  double *amin = malloc(npl*sofd);
  for (i=0; i<npl; i++) {
    double *kepelementsin = statetokep(yin[0+i*6], yin[1+i*6], yin[2+i*6], yin[3+i*6], yin[4+i*6], yin[5+i*6], msys[i]);
    double aorig = kepelementsin[0];
    free(kepelementsin);
    amax[i] = aorig*(1.0+AFRACSTABLE);
    amin[i] = aorig*(1.0-AFRACSTABLE);
  }
#endif
 
 /* Setup forward integration */
  seteq(npl, y, yin);
  double dp[npl], ddpdt, dpold[npl] ;  /* dp = dot-product = x.v */ 
  for(pl=0; pl<npl; pl++) dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
  double ps2, dist2;  /* ps2 = projected separation squared. */ 
  t = tstart;
  double h = HStart;


  int k;

  long maxtransits = 10000;
  int toomany=0;
  int ntarrelem = 6;
  double *transitarr = malloc((maxtransits+1)*ntarrelem*sofd);
  if (transitarr == NULL) {
    printf("Allocation Error\n");
    exit(0);
  } 
  long ntransits = 0;

  long vv = 0;

  // This messes up reduced chi squared if your rv values are outside of the integration time!!!  
  // But we are only comparing relative chi squareds
  // Plus just don't include rv times that aren't within your integration bounds, because that's silly
  double *rvarr;
  double *rvtarr;
  int onrv = 0;
  int rvcount, startrvcount; 
  if (RVS) {
    vv = (long) tve[0][0];
    rvarr = calloc(vv+1, sofd);
    rvarr[0] = (double) vv;
    rvtarr = malloc((vv+2)*sofd);
    rvtarr[0] = -HUGE_VAL;
    rvtarr[vv+1] = HUGE_VAL;
    memcpy(&rvtarr[1], &tve[0][1], vv*sofd);
    startrvcount = 1;
    while (rvtarr[startrvcount] < tstart) startrvcount+=1;
    rvcount = startrvcount;
  }
  long vvt;
  double *ttvarr;
  double **ttvarr2;
  if (TTVCHISQ) {
    ttvarr2 = malloc(npl*sofds);
    for (i=0; i<npl; i++) {
      ttvarr2[i] = malloc(NTTV[i][0]*sofd);
    }
    vvt = (long) nte[0][0];
    ttvarr = malloc((vvt+1)*sofd);
    ttvarr[0] = (double) vvt;

  }


#if ( demcmc_compile==3)

  int stabilitystatus=0;

  double **orbarray = malloc(npl*sofds);
  int onorbout = 0;
  double orbt = tstart+TINTERVAL/2.0;
  int densestart = 1;
  double denseinterval=0.1;
  double densemax=5000.0;
  if (densestart) {
    orbt = tstart;// + denseinterval;// /2.0;
  }
  long nposorbouts=0;
  long nnegorbouts=0;

  char stabstr[1000];
  strcpy(stabstr, "stability_");
  strcat(stabstr, OUTSTR);
  strcat(stabstr, ".out");
  FILE *stabout = fopen(stabstr, "w");
  fprintf(stabout, "t\t\t\tpnum\t\tP\t\ta\t\te\t\ti\t\t\tlo\t\t\tbo\t\t\tf\n");
  char stabstr2[1000];
  strcpy(stabstr2, "stability_xyz_");
  strcat(stabstr2, OUTSTR);
  strcat(stabstr2, ".out");
  FILE *stabout2 = fopen(stabstr2, "w");
  fprintf(stabout2, "t\t\t\tpnum\t\tx\t\ty\t\tz\t\tvx\t\tvy\t\tvz\t\t\n");

  double tfirst, tlast;

#endif

#if ( demcmc_compile == 0 || demcmc_compile == 3)
  FILE **directory = malloc(npl*sizeof(FILE*));
  for (i=0; i<npl; i++) {
    char str[30];
    sprintf(str, "tbv00_%02i.out", i+1);
    directory[i]=fopen(str,"w");
  }
  double tbvmax = 36500.0;
#endif

  double **tmte;
  double **rvtmte;

  double dist2divisor=DIST2DIVISOR;
  double baryz;
  double eps=1e-10;

  double printepoch=PRINTEPOCH;

#if (demcmc_compile != 3)
  for (i=0; i<npl; i++) {
    int ih;
    for (ih=0; ih<6; ih++) {
      if (isnan(yin[i*6+ih])) {
        goto exitgracefully;
      }
    }
  }
#endif
 
  int *ttvi;
  int *ttviinit;
  if (TTVCHISQ) {
    ttvi = calloc(npl, sofi);
    ttviinit = malloc(npl*sofi);
    for (i=0; i<npl; i++) {
      while(NTTV[i][ttvi[i]+1] < 0) {
        ttvi[i]+=1;
      }
      ttviinit[i] = ttvi[i];
    }
  }

#if ( demcmc_compile==0 )
  printf("starting integration.\n");
  printf("t=%lf tepoch=%lf t1=%lf t0=%lf h=%lf\n", t, tstart, t1, t0, h); 
#endif


  // Main forward integration loop
  while (t < t1 ) {      

#if ( demcmc_compile==0 )
    int printtime=0;
    double htimehere;
    if (h+t > printepoch) {
      printtime=1;
      htimehere=h;
      h = printepoch-t;
      if (h<0) {
        printf("Error, bad printepoch\n");
      }
      printepoch=HUGE_VAL;
    }
#endif


    if (RVS) {
      onrv=0;
      if (h + t > rvtarr[rvcount]) {  
        onrv = 1;
        hhere = h;
        h = rvtarr[rvcount] - t;
      }
    }

#if ( demcmc_compile==3)
    onorbout = 0;
    if (h + t > orbt) {
      onorbout = 1;
      hhere = h;
      h = orbt - t;
    }
#endif

    int status = gsl_odeiv_evolve_apply (e, c, s,
                                         &sys, 
                                         &t, t1,
                                         &h, y);

    if (status != GSL_SUCCESS)
        break;

#if ( demcmc_compile==0 )
    if (printtime) {
      char outfile2str[80];
      strcpy(outfile2str, "xyz_adjusted_");
      strcat(outfile2str, OUTSTR);
      strcat(outfile2str, ".pldin");
      FILE *outfile2 = fopen(outfile2str, "a");
      fprintf(outfile2, "planet                 x                        y                      z                      v_x                   v_y                   v_z                       m                      rpors             ");
      fprintf(outfile2, "\n");
      double pnum = 0.1;
      for (i=0; i<NPL; i++) {
        fprintf(outfile2, "%1.1lf", pnum);
        int j;
        for (j=0; j<6; j++) {
          fprintf(outfile2, "\t%.15lf", y[6*i+j]);
        }
        fprintf(outfile2, "\t%.15lf", mu[i+1]*MSOMJ);
        fprintf(outfile2, "\t%.15lf", rad[i]);
  
        fprintf(outfile2, "\n");
        pnum+=0.1;
      }
      fprintf(outfile2, "%.15lf ; mstar\n", mu[0]);
      fprintf(outfile2, "%.15lf ; rstar\n", rstar);
      fprintf(outfile2, "%.15lf ; c1\n", c1);
      fprintf(outfile2, "%.15lf ; c2\n", c2);
      fprintf(outfile2, "%.15lf ; dilute\n", dilute);
      fprintf(outfile2, " ; These coordinates are stellar centric\n");
      fprintf(outfile2, " ; Tepoch = %0.15lf\n", PRINTEPOCH);

      fclose(outfile2);

      printtime=0;
      h=htimehere;
    }
#endif

#if ( demcmc_compile==3)
    if (onorbout == 1) {
      h = hhere;
      for (i=0; i<npl; i++) {
        orbarray[i] = statetokep(y[0+i*6], y[1+i*6], y[2+i*6], y[3+i*6], y[4+i*6], y[5+i*6], msys[i]);
        double *aeielem = keptoorb(orbarray[i][0], orbarray[i][1], orbarray[i][2], orbarray[i][3], orbarray[i][4], orbarray[i][5], msys[i]);
        fprintf(stabout, "%.13lf\t%i\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", orbt, i, aeielem[0], orbarray[i][0], orbarray[i][1], orbarray[i][2]*180.0/M_PI, orbarray[i][3]*180.0/M_PI, orbarray[i][4]*180.0/M_PI, orbarray[i][5]*180.0/M_PI);
        fprintf(stabout2, "%.13lf\t%i\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", orbt, i, y[0+i*6], y[1+i*6], y[2+i*6], y[3+i*6], y[4+i*6], y[5+i*6]);
        if (orbarray[i][0] > amax[i] || orbarray[i][0] < amin[i]) {
          stabilitystatus = 1; 
          tlast = orbt;
        }
        free(orbarray[i]);
        free(aeielem);
      }
      if (densestart) {
        if (t<densemax) {
          orbt += denseinterval;
        } else {
          if (t<365.*100000.) orbt += TINTERVAL;
          else orbt += TINTERVAL*50.;
        }
      } else {
        if (t<365.*100000.) orbt += TINTERVAL;
        else orbt += TINTERVAL*50.;
      }
      nposorbouts++;
      if (nposorbouts % 1 == 0) fflush(stabout);
      if (nposorbouts % 1 == 0) fflush(stabout2);
    }
    if (stabilitystatus == 1) { 
      printf("Went Unstable -> Break\n");
      break;
    }
#endif


    if (RVS) {
      if (onrv ==1) {
        h = hhere;
        func(t, y, dydt, mu);
        double vstar = 0;
        for (i=0; i<npl; i++) {
          vstar += -y[i*6+5]*mu[i+1]/mtot;
        }
        rvarr[rvcount] = vstar;    
        rvcount += 1;
      }
    }


    dist2 = pow(y[0],2)+pow(y[1],2)+pow(y[2],2);
    for(pl=0; pl<npl; pl++) {  /* Cycle through the planets, searching for a transit. */
      dpold[pl]=dp[pl];
      dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
      ps2=y[0+pl*6]*y[0+pl*6]+y[1+pl*6]*y[1+pl*6];

      if( /* 0 == 1 */ dp[pl]*dpold[pl] <= 0 && ps2 < dist2/dist2divisor && y[2+pl*6] < 0 ) { 
                     /* A minimum projected-separation occurred: "Tc".  
                   ps2 constraint makes sure it's when y^2+z^2 is small-- its on the face of the star.  
                   y[2] constraint means a primary eclipse, presuming that the observer is on the positive x=y[2] axis. */

        transitcount[pl] += 1;  /* Update the transit number. */

        seteq(npl, yhone, y);
        thone = t;
        i=0;
        do { 
          func (thone, yhone, dydt, npl_mu);
          ddpdt = dydt[0+pl*6]*dydt[0+pl*6] + dydt[1+pl*6]*dydt[1+pl*6] + yhone[0+pl*6]*dydt[3+pl*6] + yhone[1+pl*6]*dydt[4+pl*6];
          tstep = - dp[pl] / ddpdt;
          gsl_odeiv_step_apply (s, thone, tstep, yhone, yerrhone, dydt, NULL, &sys);
          thone += tstep;
          dp[pl]=yhone[0+pl*6]*yhone[3+pl*6]+yhone[1+pl*6]*yhone[4+pl*6];
          i++;
        } while (i<5 && fabs(tstep)>eps);
        
        dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];

        double zp=0;
        double zb=0;
        if (LTE) {
          baryz = 0;
          for(i=0; i<npl; i++) baryz += yhone[2+i*6]*mu[i+1];
          baryz /= mtot; 
          zp = yhone[2+pl*6]/CAUPD;
          zb = baryz/CAUPD;
          zp -= zb;
        }
        double t2 = thone - zb + zp;


#if ( (demcmc_compile == 0) )//|| (demcmc_compile == 3) )
        dist = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
        vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2)); 
        fprintf(directory[pl], "%6d  %.18e  %.10e  %.10e\n", transitcount[pl], t2, dist, vtrans);
#endif

        if (TTVCHISQ) {
          if (transitcount[pl] == NTTV[pl][ttvi[pl]+1]) {
#if (demcmc_compile==0)
            MTTV[pl][ttvi[pl]+1]=t2;
#endif
            ttvarr2[pl][ttvi[pl]]=t2;
            ttvi[pl]+=1;
          }
        }

#if ( (demcmc_compile == 3) )
        if (t2 < tbvmax) {
          double bsign=1.;
          if (yhone[1+pl*6] < 0.) bsign=-1;
          dist = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
          vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2)); 
          fprintf(directory[pl], "%6d  %.18e  %.10e  %.10e\n", transitcount[pl], t2, bsign*dist, vtrans);
        }
#endif

#if (demcmc_compile != 3)
        transitarr[ntransits*ntarrelem+0] = t2;
        transitarr[ntransits*ntarrelem+1] = yhone[0+pl*6]; 
        transitarr[ntransits*ntarrelem+2] = yhone[1+pl*6];
        transitarr[ntransits*ntarrelem+3] = yhone[3+pl*6];
        transitarr[ntransits*ntarrelem+4] = yhone[4+pl*6];
        transitarr[ntransits*ntarrelem+5] = rad[pl];
        ntransits++;
        if (ntransits > maxtransits) {
          printf("Too many transits - increase maxtransits\n");
          fflush(stdout);
          toomany=1;
          goto exitgracefully; 
        }
#endif
      }
    }
  }


  long postransits;
  postransits = ntransits;


  /* Setup backward integration */
  seteq(npl, y, yin);
  for(pl=0; pl<npl; pl++) dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6]; /* dp = dot-product = x.v */
  t = tstart;
  h = -HStart;
  for(pl=0;pl<npl;pl++) transitcount[pl] = 0;

  if (RVS) {
    rvcount = startrvcount-1;
  }

#if ( demcmc_compile==3)
  int posstability = 0;
  if (stabilitystatus==1) posstability=1;
  stabilitystatus=0;
  if (densestart) {
    orbt = t-denseinterval;// /2.0;
  } else {
    orbt = t-TINTERVAL;// /2.0;
  }
#endif
  
  if (TTVCHISQ) {
    for (i=0; i<npl; i++) {
      ttvi[i] = ttviinit[i]-1;
    }
  }

#if ( demcmc_compile==0 )
  printf("RE starting integration.\n");
  printf("t=%lf tepoch=%lf t1=%lf t0=%lf h=%lf\n", t, tstart, t1, t0, h); 
#endif

  while (t > t0+1) {      
    if (RVS) {
      onrv=0;
      if ( (h + t) < rvtarr[rvcount]) {  
        onrv = 1;
        hhere = h;
        h = rvtarr[rvcount] - t;
      }
    }

#if ( demcmc_compile==3)
    onorbout = 0;
    if (h + t < orbt) {
      onorbout = 1;
      hhere = h;
      h = orbt - t;
    }
#endif

    int status = gsl_odeiv_evolve_apply (e, c, s,
                                         &sys, 
                                         &t, t0,
                                         &h, y);

    if (status != GSL_SUCCESS)
        break;

#if ( demcmc_compile==3)
    if (onorbout == 1) {
      h = hhere;
      for (i=0; i<npl; i++) {
        orbarray[i] = statetokep(y[0+i*6], y[1+i*6], y[2+i*6], y[3+i*6], y[4+i*6], y[5+i*6], msys[i]);
        double *aeielem = keptoorb(orbarray[i][0], orbarray[i][1], orbarray[i][2], orbarray[i][3], orbarray[i][4], orbarray[i][5], msys[i]);
        fprintf(stabout, "%.13lf\t%i\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", orbt, i, aeielem[0], orbarray[i][0], orbarray[i][1], orbarray[i][2]*180.0/M_PI, orbarray[i][3]*180.0/M_PI, orbarray[i][4]*180.0/M_PI, orbarray[i][5]*180.0/M_PI);
        fprintf(stabout2, "%.13lf\t%i\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", orbt, i, y[0+i*6], y[1+i*6], y[2+i*6], y[3+i*6], y[4+i*6], y[5+i*6]);
        if (orbarray[i][0] > amax[i] || orbarray[i][0] < amin[i]) {
          stabilitystatus = 1; 
          tlast = orbt;
        }
        free(orbarray[i]);
        free(aeielem);
      }
      if (densestart) {
        orbt -= denseinterval;
      } else {
        orbt -= TINTERVAL;
      }
    }
    if (stabilitystatus == 1) { 
      printf("Went Unstable -> Break\n");
      break;
    }
#endif

    if (RVS) {
      if (onrv ==1) {
        h = hhere;
        func(t, y, dydt, mu);
        double vstar = 0;
        for (i=0; i<npl; i++) {
          vstar += -y[i*6+5]*mu[i+1]/mtot;
        }
        rvarr[rvcount] = vstar;    
        rvcount -= 1;
      }
    }


    dist2 = pow(y[0],2)+pow(y[1],2)+pow(y[2],2);

    for(pl=0; pl<npl; pl++) {  /* Cycle through the planets, searching for a transit. */
      dpold[pl]=dp[pl];
      dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];
      ps2=y[0+pl*6]*y[0+pl*6]+y[1+pl*6]*y[1+pl*6];

      if( /* 0 == 1 */ dp[pl]*dpold[pl] <= 0 && ps2 < dist2/dist2divisor && y[2+pl*6] < 0 ) { 
                /* A minimum projected-separation occurred: "Tc".  
                   ps2 constraint makes sure it's when y^2+z^2 is small-- its on the face of the star.  
                   y[0] constraint means a primary eclipse, presuming that the observer is on the positive x=y[0] axis. */

        transitcount[pl] -= 1;  /* Update the transit number. */

        seteq(npl, yhone, y);
        thone = t;
        i=0;
        do {
          func (thone, yhone, dydt, npl_mu);
          ddpdt = dydt[0+pl*6]*dydt[0+pl*6] + dydt[1+pl*6]*dydt[1+pl*6] + yhone[0+pl*6]*dydt[3+pl*6] + yhone[1+pl*6]*dydt[4+pl*6];
          tstep = - dp[pl] / ddpdt;
          gsl_odeiv_step_apply (s, thone, tstep, yhone, yerrhone, dydt, NULL, &sys);
          thone += tstep;
          dp[pl]=yhone[0+pl*6]*yhone[3+pl*6]+yhone[1+pl*6]*yhone[4+pl*6];
          i++;
        } while (i<5 && fabs(tstep)>eps); 
          dp[pl]=y[0+pl*6]*y[3+pl*6]+y[1+pl*6]*y[4+pl*6];

        double zp=0;
        double zb=0;
        if (LTE) {
          baryz = 0;
          for(i=0; i<npl; i++) baryz += yhone[2+i*6]*mu[i+1];
          baryz /= mtot; 
          zp = yhone[2+pl*6]/CAUPD;
          zb = baryz/CAUPD;
          zp -= zb;
        }
        double t2 = thone - zb + zp;

#if ( demcmc_compile == 0 )//|| demcmc_compile == 3)
        dist = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
        vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2)); 
        fprintf(directory[pl], "%6d  %.18e  %.10e  %.10e\n", transitcount[pl], t2, dist, vtrans);
#endif

        if (TTVCHISQ) {
          if (ttvi[pl] >= 0 && transitcount[pl] == NTTV[pl][ttvi[pl]+1]) { 
#if (demcmc_compile==0)
            MTTV[pl][ttvi[pl]+1]=t2;
#endif
            ttvarr2[pl][ttvi[pl]]=t2;
            ttvi[pl]-=1;
          }
        }
#if ( (demcmc_compile == 3) )
        if (t2 < tbvmax) {
          double bsign=1.;
          if (yhone[1+pl*6] < 0.) bsign=-1;
          dist = sqrt(pow(yhone[0+pl*6],2)+pow(yhone[1+pl*6],2));
          vtrans = sqrt(pow(yhone[3+pl*6],2)+pow(yhone[4+pl*6],2)); 
          fprintf(directory[pl], "%6d  %.18e  %.10e  %.10e\n", transitcount[pl], t2, bsign*dist, vtrans);
        }
#endif

#if (demcmc_compile != 3)
        transitarr[ntransits*ntarrelem+0] = t2;
        transitarr[ntransits*ntarrelem+1] = yhone[0+pl*6]; 
        transitarr[ntransits*ntarrelem+2] = yhone[1+pl*6];
        transitarr[ntransits*ntarrelem+3] = yhone[3+pl*6];
        transitarr[ntransits*ntarrelem+4] = yhone[4+pl*6];
        transitarr[ntransits*ntarrelem+5] = rad[pl];
        ntransits++;
        if (ntransits > maxtransits) {
          printf("Too many transits - increase maxtransits\n");
          fflush(stdout);
          toomany=1;
          goto exitgracefully; 
        }
#endif

      }
    }
  }

  if (RVS) {
    free(rvtarr);
  }

  long negtransits;
  negtransits = ntransits-postransits;
  
#if ( demcmc_compile==3)
  int negstability = 0;
  if (stabilitystatus==1) negstability=1;
  if (negstability) fprintf(stabout, "System went unstable at t=%lf as 1 or more planets had a deviation from its original semi-major axis by %lf or more\n", tfirst, AFRACSTABLE); 
  if (posstability) fprintf(stabout, "%lf : System went unstable at t=%lf as 1 or more planets had a deviation from its original semi-major axis by %lf or more\n", tlast, tlast, AFRACSTABLE);
  fclose(stabout);
  fclose(stabout2);
  free(amin); free(amax);
  free(orbarray);
#endif


#if ( demcmc_compile == 0 || demcmc_compile == 3)

  for (i=0; i<npl; i++){
    fclose(directory[i]);
  }
  free(directory);

#endif


#if (demcmc_compile != 3)
  long timelistlen; 
  timelistlen = kk;

  double *timelist;
  double *temptimelist;
  double *fluxlist;
  double *errlist;
 
  double *rvtimelist;
  double *rvlist;
  double *rverrlist;

  //array of times, measured, and theory data, and error
  tmte = malloc(4*sofds);
  tmte[2] = malloc((kk+1)*sofd);

  temptimelist= &tfe[0][1];
  fluxlist= &tfe[1][1];
  errlist= &tfe[2][1];
  tmte[0] = &tfe[0][0];
  tmte[1] = &tfe[1][0];
  tmte[3] = &tfe[2][0];

  int tele;
  int maxteles;
  maxteles=10;
  double *rvoffset;
  rvoffset = malloc(maxteles*sofd);
  if (RVS) {
    // Compute RV offset
    double rvdiff;
    long vvw;
    double weight;

    long vvt=1;
    tele=0;
    while (vvt<(vv+1)) {
      double numerator = 0;
      double denominator = 0;
      for (vvw=1; vvw<(vv+1); vvw++) {
        if ( (int) tve[4][vvw] == tele ) {
          rvdiff = rvarr[vvw] - tve[1][vvw];
          weight = 1.0 / pow(tve[2][vvw], 2);
          numerator += rvdiff*weight;
          denominator += weight;
        }
      }
      rvoffset[tele] = numerator/denominator;
      for (vvw=1; vvw<(vv+1); vvw++) {
        if ( (int) tve[4][vvw] == tele ) {
          rvarr[vvw] -= rvoffset[tele];
          vvt+=1;
        }
      }
      tele+=1;
    }
  }

  rvtmte = malloc(4*sofds);

  rvtimelist=&tve[0][1];
  rvlist=&tve[1][1];
  rverrlist=&tve[2][1];
  rvtmte[0] = &tve[0][0];
  rvtmte[1] = &tve[1][0];
  rvtmte[3] = &tve[2][0];
  rvtmte[2] = rvarr;
  
  
  double **ttvtmte;
  if (TTVCHISQ) {
    int ki;
    int ksofar=0;
    for (i=0; i<npl; i++) {
      for (ki=0; ki<NTTV[i][0]; ki++) {
        ttvarr[ksofar+ki+1] = ttvarr2[i][ki];
      }
      ksofar += NTTV[i][0];
    }
    for (i=0; i<npl; i++) {
      free(ttvarr2[i]);
    }
    free(ttvarr2);
    ttvtmte = malloc(4*sofds);
    ttvtmte[0] = &nte[0][0];
    ttvtmte[1] = &nte[1][0];
    ttvtmte[3] = &nte[2][0];
    ttvtmte[2] = ttvarr;
  }

  double *temptlist; 
  double *timedlc ( double *times, int *cadences, long ntimes, double **transitarr, int nplanets, double rstar, double c1, double c2); //prototype
  double *binnedlc ( double *times, int *cadences, long ntimes, double binwidth, int nperbin,  double **transitarr, int nplanets, double rstar, double c1, double c2);

  long marker, l;
  int ntran;
  l=0;  
  marker = 0;
  ntran=1;
  double vel;
  double lhs;
  double rhs;


  qsort(transitarr, ntransits, (ntarrelem*sofd), compare);

  double *tclist;
  double *order;

  long q;
  if (cadenceswitch==2) {
    if (OOO == 2) {
      OOO = 0;
      for (q=0; q<kk-1; q++) {
        if (temptimelist[q+1] < temptimelist[q]) {
          OOO = 1;
          break;
        }
      }
    }
    if (OOO) {
      
      timelist = malloc(kk*sofd);
      long ti;
      for (ti=0; ti<kk; ti++) {
        timelist[ti] = tfe[0][ti+1];
      }

      int nsort=3;
      tclist = malloc((nsort*sofd)*kk); 
      order = malloc(sofd*kk);
      for (q=0; q<kk; q++) {
        tclist[nsort*q] = timelist[q];
        tclist[nsort*q+1] = (double) q;
        tclist[nsort*q+2] = (double) cadencelist[q];
      }
      qsort(tclist, kk, sofd*nsort, compare);
      for (q=0; q<kk; q++) {
        timelist[q] = tclist[nsort*q];
        cadencelist[q] = (int) tclist[nsort*q+2];
        order[q] = (long) tclist[nsort*q+1];
      }
      free(tclist);
    } else {
      timelist = temptimelist;
    }
  } else {
    timelist = temptimelist;
  }

  double omdilute;
  omdilute = 1.0-dilute;

  // cycle through each transit w/ index l
  while (l<ntransits-1) {
    // velocity of planet relative to star
    vel = sqrt( pow(transitarr[l*ntarrelem+3],2) + pow(transitarr[l*ntarrelem+4],2) );
    if (vel <= 0) {
      printf("Error, invalid planet velocity\n");
      exit(0);
    }
    // check if transits overlap or are close
    lhs = transitarr[l*ntarrelem+0] + 5.*rstar*RSUNAU/vel;
    rhs = transitarr[(l+1)*ntarrelem+0];
    // if they don't, compute the light curve (either w/ 1 planet or however many previous overlaps there were)
    if ( lhs<= rhs) {
      double **temparr = malloc(sofds*ntran);
      double minvel=vel;
      int ii;
      for (ii=0; ii<ntran; ii++) {
        temparr[ii] = malloc(ntarrelem*sofd);
        memcpy(temparr[ii], &transitarr[((l-ntran)+1+ii)*ntarrelem], ntarrelem*sofd);
        minvel = fmin(minvel, sqrt(pow(temparr[ii][3],2) + pow(temparr[ii][4],2)));
      }

      double startt, stopt;
      double *trantlist;
      int *tranclist;
      startt = temparr[0][0] - 3.*rstar*RSUNAU/minvel;
      stopt = temparr[ntran-1][0] + 3.*rstar*RSUNAU/minvel;
      while (timelist[marker]<startt && marker<kk-1) {
        tmte[2][marker+1]=1.0;       
        marker++;
      }
      long ntimes=0; 
      while (timelist[marker+ntimes]<stopt && marker+ntimes<kk-1) {
        ntimes++; 
      } 
      if (ntimes!=0) {
        trantlist = &timelist[marker];
        if (cadenceswitch==2)  tranclist = &cadencelist[marker];
        if (cadenceswitch == 1 || cadenceswitch == 2) {
          temptlist = binnedlc( trantlist, tranclist, ntimes, binwidth, nperbin, temparr, ntran, rstar, c1, c2); 
        }
        if (cadenceswitch == 0) {
          temptlist = timedlc( trantlist, tranclist, ntimes, temparr, ntran, rstar, c1, c2);
        }
        if (temptlist!=NULL) {
          memcpy(&tmte[2][marker+1],&temptlist[0], ntimes*sofd);
          free(temptlist);
        }
      }
      marker+=ntimes;
      for (ii=0; ii<ntran; ii++) {
        free(temparr[ii]);
      } 
      free(temparr);
      ntran=1;
    } else {
      ntran+=1;
    }
    l++;
  }


  //last case - compute last transit (w/ any overlaps)
  if (ntransits>0) {
    l=ntransits-1;
    double **temparr = malloc(sofds*ntran);
    vel = sqrt( pow(transitarr[l*ntarrelem+3],2) + pow(transitarr[l*ntarrelem+4],2) );
    double minvel=vel;
    int ii;
    for (ii=0; ii<ntran; ii++) {
      temparr[ii] = malloc(ntarrelem*sofd);
      memcpy(temparr[ii], &transitarr[((l-ntran)+1+ii)*ntarrelem], ntarrelem*sofd);
      minvel = fmin(minvel, sqrt(pow(temparr[ii][3],2) + pow(temparr[ii][4],2)));
    }
    double startt, stopt;
    double *trantlist;
    int *tranclist;
    startt = temparr[0][0] - 3.*rstar*RSUNAU/minvel;
    stopt = temparr[ntran-1][0] + 3.*rstar*RSUNAU/minvel;
    while (timelist[marker]<startt && marker<kk-1) {
      tmte[2][marker+1]=1.0;       
      marker++;
    }
    long ntimes=0;
    while (timelist[marker+ntimes]<stopt && marker+ntimes<kk-1) {
      ntimes++; 
    } 
    if (ntimes!=0) {
      trantlist = &timelist[marker];
      if (cadenceswitch==2)  tranclist = &cadencelist[marker];
      if (cadenceswitch == 1 || cadenceswitch == 2) {
        temptlist = binnedlc( trantlist, tranclist, ntimes, binwidth, nperbin, temparr, ntran, rstar, c1, c2); 
      }
      if (cadenceswitch == 0) {
        temptlist = timedlc( trantlist, tranclist, ntimes, temparr, ntran, rstar, c1, c2);
      }
      if (temptlist!=NULL) {
        memcpy(&tmte[2][marker+1],&temptlist[0], ntimes*sofd);
        free(temptlist);
      }
    }
    marker+=ntimes;
  
    for (ii=0; ii<ntran; ii++) {
      free(temparr[ii]);
    } 
    free(temparr);
  }

  while (marker<kk) {
    tmte[2][marker+1] = 1.0;
    marker++;
  }
 
  if ( !(dbleq(dilute, 0.0)) ) {
    for (l=0; l<kk; l++) {
      tmte[2][l] = tmte[2][l]*omdilute+dilute;
    }
  }

  if (cadenceswitch==2) {
    if (OOO) { 
      int nsort=4;
      tclist = malloc((nsort*sofd)*kk); 
      for (q=0; q<kk; q++) {
        tclist[nsort*q] = (double) order[q]; //timelist[q];
        tclist[nsort*q+1] = tmte[2][q+1]; //  (double) q;
        tclist[nsort*q+3] = (double) cadencelist[q]; //  (double) q;
      }
      qsort(tclist, kk, sofd*nsort, compare);
      for (q=0; q<kk; q++) {
        tmte[2][q+1] = tclist[nsort*q+1];
        cadencelist[q] = (int) tclist[nsort*q+3];
      }
      free(tclist);
      free(order);
      free(timelist);
    }
  }




  tmte[0][0] = kk;
  tmte[1][0] = kk;
  tmte[2][0] = kk;
  tmte[3][0] = kk;



#if (demcmc_compile==0)

  char tmtefstr[1000];
  strcpy(tmtefstr, "lc_");
  strcat(tmtefstr, OUTSTR);
  strcat(tmtefstr, ".lcout");
  FILE *tmtef;
  tmtef = fopen(tmtefstr, "a");
  long ijk;
  if (CADENCESWITCH==2) {
    for (ijk=1; ijk<=kk; ijk++) {
      fprintf(tmtef, "%lf %lf %lf %lf %i\n", tmte[0][ijk], tmte[1][ijk], tmte[2][ijk], tmte[3][ijk], cadencelist[ijk-1]);
    }
  } else {
    for (ijk=1; ijk<=kk; ijk++) {
      fprintf(tmtef, "%lf %lf %lf %lf\n", tmte[0][ijk], tmte[1][ijk], tmte[2][ijk], tmte[3][ijk]);
    }
  }
  fclose(tmtef);

  int nbodies;
  nbodies=npl+1;
  if (RVS) {
    char rvtmtefstr[1000];
    char num[10];
    strcpy(rvtmtefstr, "rv");
    strcat(rvtmtefstr, "00");
    strcat(rvtmtefstr, "_");
    strcat(rvtmtefstr, OUTSTR);
    strcat(rvtmtefstr, ".rvout");
    FILE *rvtmtef = fopen(rvtmtefstr, "a");
    for (ijk=1; ijk<vv+1; ijk++) {
      fprintf(rvtmtef, "%lf\t%e\t%e\t%e\t%i\n", rvtmte[0][ijk], rvtmte[1][ijk]/MPSTOAUPD, rvtmte[2][ijk]/MPSTOAUPD, rvtmte[3][ijk]/MPSTOAUPD, (int) tve[4][ijk]);
    }
    int l;
    for (l=0; l<tele; l++) fprintf(rvtmtef, "RV offset %i = %lf\n", l, rvoffset[l]/MPSTOAUPD);
    fclose(rvtmtef);
  }
#endif


  free(rvoffset);

  exitgracefully:;
  if (toomany) {
    tmte = malloc(4*sofds);
    tmte[2] = calloc((kk+1),sofd);
    tmte[2][1] = HUGE_VAL;
    tmte[2][0] = (double) kk;
    tmte[0] = &tfe[0][0];
    tmte[1] = &tfe[1][0];
    tmte[3] = &tfe[2][0];
    tmte[0][0] = (double) kk;
    tmte[1][0] = (double) kk;
    tmte[3][0] = (double) kk;
    rvtmte = malloc(4*sofds);
    rvtmte[0] = &tve[0][0];
    rvtmte[1] = &tve[1][0];
    rvtmte[3] = &tve[2][0];
    if (RVS) {
      rvtmte[2] = calloc(vv+1, sofd);
      rvtmte[2][1] = HUGE_VAL;
      rvtmte[2][0] = (double) vv;
      rvtmte[0][0] = (double) vv;
      rvtmte[1][0] = (double) vv;
      rvtmte[3][0] = (double) vv;
    }
  }

  if (TTVCHISQ) {
    free(ttvi);
    free(ttviinit);
  }

  free(rad); 
  free(transitarr);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  double ***fl_rv = malloc(3*sizeof(double**));
  fl_rv[0] = tmte;
  fl_rv[1] = rvtmte;
  fl_rv[2] = ttvtmte;

  
  return fl_rv;

#else

  free(rad); 
  free(transitarr);
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  double ***empty;
  return empty;

#endif


}



// Read the formatted input file
// Sets global variables
int getinput(char fname[]) {

  const int sofd = SOFD;
  const int sofi = SOFI;

  FILE *inputf = fopen(fname, "r"); 
  if (inputf == NULL) {
    printf("Error: Bad Input File Name\n");
    printf("       Either no file was passed (see README for usage) or\n");
    printf("       the specified file was not found in the run directory.\n");
    exit(0);
  }

  OUTSTR = malloc(1000*sizeof(char));

  char buffer[1000];
  char type[100];
  char varname[100];
  int i;
  int j;
  for (i=0; i<8; i++) fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %s", type, varname, &OUTSTR[0]); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &NBODIES); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);

  NPL = NBODIES-1;
  if (NPL < 0) {
    printf("NPL = %i\n");
    printf("Error: At least one planet must be present.\n");
    exit(0);
  }

  fscanf(inputf, "%s %s %lf", type, varname, &EPOCH); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %li", type, varname, &NWALKERS); fgets(buffer, 1000, inputf);
#if (demcmc_compile == 1)
  if (NWALKERS <= 2) {
    printf("Error: At least 3 walkers are required for the DEMCMC algorithm to function\n");
    exit(0);
  }
  int ncores;
  MPI_Comm_size(MPI_COMM_WORLD, &ncores);
  if ((NWALKERS % ncores) > 0) {
    printf("ERROR: Nwalkers modulo Ncores must be 0.\n");
    printf("       Nwalkers=%i, Ncores=%i\n", NWALKERS, ncores);
    printf("       At least one of these must be changed.\n");
    exit(0);
  }
#endif
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %li", type, varname, &NSTEPS); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &CADENCESWITCH); fgets(buffer, 1000, inputf); 
  printf("cadenceswitch = %i\n", CADENCESWITCH);
  if (!(CADENCESWITCH == 0 || CADENCESWITCH == 1 || CADENCESWITCH == 2)) {
    printf("Error: Cadenceswitch takes only values of 0, 1, or 2.\n");
    exit(0);
  }
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %s", type, varname, TFILE); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &MULTISTAR); fgets(buffer, 1000, inputf); 

  if (MULTISTAR) {
    printf("WARNING!\n");
    printf("Multiple-stars are disabled in this version of the code\n");
  }
  if (MULTISTAR) PPERPLAN = 11;
  else PPERPLAN = 8; 
  printf("pperplan=%i\n", PPERPLAN);
  
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i %i %i", type, varname, &RVJITTERFLAG, &NSTARSRV, &NTELESCOPES); fgets(buffer, 1000, inputf);
  printf("rvjitterflag, NstarsRV, Ntelescopes =  %i, %i, %i\n", RVJITTERFLAG, NSTARSRV, NTELESCOPES);
  if (RVJITTERFLAG && NSTARSRV > 1) {
    printf("Error: Multiple-stars are disabled in this version of the code, but you set NSTARSRV=%i\n", NSTARSRV);
    exit(0);
  }
  if (RVJITTERFLAG && (NSTARSRV < 0 || NTELESCOPES < 1)) {
    printf("Error, invalid RV jitter parameter.\n");
    exit(0);
  }
  if (RVJITTERFLAG) {
    RVJITTERTOT = NSTARSRV*NTELESCOPES;
    PSTAR += RVJITTERTOT;//*2
  }
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &TTVCHISQ); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i %i %i", type, varname, &TTVJITTERFLAG, &NSTARSTTV, &NTELESCOPESTTV); fgets(buffer, 1000, inputf);
  if (TTVJITTERFLAG) {
    TTVJITTERTOT = NSTARSTTV*NTELESCOPESTTV;
    PSTAR += TTVJITTERTOT;//*2
  }
  if (TTVJITTERFLAG && NSTARSTTV > 1) {
    printf("TTVjitterflag, NstarsTTV, Ntelescopes =  %i, %i, %i\n", TTVJITTERFLAG, NSTARSTTV, NTELESCOPESTTV);
    printf("Error: Multiple-stars are disabled in this version of the code, but you set NSTARSTTV=%i\n", NSTARSTTV);
    exit(0);
  }
  if (TTVJITTERFLAG && (NSTARSTTV < 0 || NTELESCOPESTTV < 1)) {
    printf("TTVjitterflag, NstarsTTV, Ntelescopes =  %i, %i, %i\n", TTVJITTERFLAG, NSTARSTTV, NTELESCOPESTTV);
    printf("Error, invalid value for TTV jitter parameters.\n");
    exit(0);
  }

  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &CELERITE); fgets(buffer, 1000, inputf);
#if (celerite_compile == 0)
  if (CELERITE) {
    printf("Error: Celerite is set to True, but the code was compiled without celerite functionality\n");
    printf("Either recompile with celerite, or turn off celerite in the .in file.\n");
    exit(0);
  }
#endif
  if (CELERITE) {
    PSTAR += NCELERITE;//*2
  }
  printf("celerite=%i\n", CELERITE);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &RVCELERITE); fgets(buffer, 1000, inputf);
#if (celerite_compile == 0)
  if (RVCELERITE) {
    printf("Error: RVCelerite is set to True, but the code was compiled without celerite functionality\n");
    printf("Either recompile with celerite, or turn off RVcelerite in the .in file.\n");
    exit(0);
  }
#endif
  if (RVCELERITE) {
    PSTAR += NRVCELERITE;//*2
  }
  printf("RVcelerite=%i\n", RVCELERITE);

  printf("pstar=%i\n", PSTAR);

  PARFIX = malloc((PPERPLAN*NPL+PSTAR)*sofi);
  PSTEP = malloc(PPERPLAN*sofd);
  SSTEP = malloc(PSTAR*sofd);
  CSTEP = malloc((PPERPLAN*NPL+PSTAR)*sofd);
  BIMODLIST = malloc((PPERPLAN*NPL+PSTAR)*sofd);
  
  MAXDENSITY = malloc(NPL*sofd);
  EMAX = malloc(NPL*sofd); 
  EPRIORV = malloc(NPL*sofd);
  MASSPRIORCENTERS = malloc(NPL*sofd);
  MASSPRIORSIGMAS = malloc(NPL*sofd);

  const int npl = NPL;

  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &SPECTROSCOPY); fgets(buffer, 1000, inputf); 
  if (SPECTROSCOPY) {
    fscanf(inputf, "%lf", &SPECRADIUS); fscanf(inputf, "%lf", &SPECERRPOS); fscanf(inputf, "%lf", &SPECERRNEG); fgets(buffer, 1000, inputf);
    printf("spectroscopy = %i, %lf, %lf, %lf\n", SPECTROSCOPY, SPECRADIUS, SPECERRPOS, SPECERRNEG);
    if (SPECRADIUS < 0. || SPECERRPOS < 0. || SPECERRNEG < 0.) {
       printf("Error: all spectroscopy values must be >= 0\n");
       exit(0);
    }
  } else {
    fgets(buffer, 1000, inputf);
    printf("radius spectroscopy off\n");
  }
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &MASSSPECTROSCOPY); fgets(buffer, 1000, inputf); 
  if (MASSSPECTROSCOPY) {
    fscanf(inputf, "%lf", &SPECMASS); fscanf(inputf, "%lf", &MASSSPECERRPOS); fscanf(inputf, "%lf", &MASSSPECERRNEG); fgets(buffer, 1000, inputf);
    printf("MASS spectroscopy = %i, %lf, %lf, %lf\n", MASSSPECTROSCOPY, SPECMASS, MASSSPECERRPOS, MASSSPECERRNEG);
    if (SPECMASS < 0. || MASSSPECERRPOS < 0. || MASSSPECERRNEG < 0.) {
       printf("Error: all spectroscopy values must be >= 0\n");
       exit(0);
    }
  } else {
    fgets(buffer, 1000, inputf);
    printf("mass spectroscopy off\n");
  }
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &LDLAW); fgets(buffer, 1000, inputf); 
  printf("LDLAW = %i\n", LDLAW);
  if (LDLAW == 1) {
    printf("Note: this implementation of the Maxsted+2018 power2 law does not account for mutual transits\n");
  }
  if (!(LDLAW == 0 || LDLAW == 1)) {
    printf("LDLAW = %i\n", LDLAW);
    printf("Error: ldlaw must be 0 or 1\n");
    exit(0);
  }
  fgets(buffer, 1000, inputf); 
  fscanf(inputf, "%s %s %i", type, varname, &DIGT0); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &IGT90); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &INCPRIOR); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &MGT0); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &DENSITYCUTON); fgets(buffer, 1000, inputf);
  if (DENSITYCUTON) { 
    for (i=0; i<npl; i++) fscanf(inputf, "%lf", &MAXDENSITY[i]); fgets(buffer, 1000, inputf); //l.38
    for (i=0; i<npl; i++) {
      if (MAXDENSITY[i] <= 0) {
        printf("Warning: MAXDENSITY[%i] is <=0, this may cause a crash or hang.\n", i);
      }
    }
  } else {
    fgets(buffer, 1000, inputf);
  }
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &MASSPRIOR); fgets(buffer, 1000, inputf);
  if (!((MASSPRIOR == 0) || (MASSPRIOR == 1))) {
    printf("Error: MASSPRIOR = %i\n", MASSPRIOR);
    printf("       Massprior must be 0 or 1\n");
    exit(0); 
  }
  if (MASSPRIOR) { 
    for (i=0; i<npl; i++) fscanf(inputf, "%lf", &MASSPRIORCENTERS[i]); fgets(buffer, 1000, inputf); 
    for (i=0; i<npl; i++) fscanf(inputf, "%lf", &MASSPRIORSIGMAS[i]); fgets(buffer, 1000, inputf); 
    for (i=0; i<npl; i++) {
      printf("Mass ratio prior for planet %i: %lf +/- %lf\n", i+1, MASSPRIORCENTERS[i], MASSPRIORSIGMAS[i]);
      if (MASSPRIORCENTERS[i] < 0) {
        printf("Warning: MASSPRIORCENTERS[%i] is <0. This will work as expected but may be unphysical.\n", i);
      }
      if (MASSPRIORSIGMAS[i] <= 0) {
        printf("Error: MASSPRIORSIGMAS[%i] = %lf\n", i, MASSPRIORSIGMAS[i]);
        printf("       Mass uncertainties must be >= 0.\n");
        if (MASSPRIORSIGMAS[i] == 0) {
           printf("To fix a mass, use parfix in the .in file rather than setting MassPriorSigma = 0.\n");
        }
        exit(0);
      }
    }
  } else {
    fgets(buffer, 1000, inputf);
    fgets(buffer, 1000, inputf);
  }
  
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &SQRTE); fgets(buffer, 1000, inputf); 
  //printf("%s\n", buffer);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &ECUTON); fgets(buffer, 1000, inputf);
  if (ECUTON) {  
    for (i=0; i<npl; i++) fscanf(inputf, "%lf", &EMAX[i]); fgets(buffer, 1000, inputf);
    for (i=0; i<npl; i++) {
      if (EMAX[i] <= 0) {
        printf("Warning: EMAX[%i] is <=0, this may cause a crash or hang.\n", i);
      }
    }
  } else {
    fgets(buffer, 1000, inputf); 
  }
  fgets(buffer, 1000, inputf); //l.45
  fgets(buffer, 1000, inputf); //l.45
  fscanf(inputf, "%s %s %i", type, varname, &EPRIOR); fgets(buffer, 1000, inputf);
  if (EPRIOR) { 
    for (i=0; i<npl; i++) fscanf(inputf, "%i", &EPRIORV[i]); fgets(buffer, 1000, inputf);
  } else {
    fgets(buffer, 1000, inputf);
  }
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &ESIGMA); fgets(buffer, 1000, inputf);
  if (ESIGMA <= 0) {
     printf("Warning ESIGMA <= 0, this may cause a crash or hang.\n");
  } 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &NTEMPS); fgets(buffer, 1000, inputf); 
  printf("sqrte, ecuton, eprior[0], esigma = %i, %i, %i, %lf\n", SQRTE, ECUTON, EPRIORV[0], ESIGMA);
  if (NTEMPS != 0){
    printf("WARNING!\n");
    printf("This version of the code stets NTEMPS=0!\n");
  }

  for (i=0; i<9; i++) fgets(buffer, 1000, inputf);
  int nfixed = 0;
  for (i=0; i<npl; i++) {
    for (j=0; j<PPERPLAN; j++) {
      fscanf(inputf, "%i", &PARFIX[PPERPLAN*i+j]);
      if (PARFIX[PPERPLAN*i+j] == 1) {
        nfixed++;
      } else {
        if (PARFIX[PPERPLAN*i+j] != 0) {
          printf("Error: All parfix values must be 0 or 1\n");
          printf("A value of %i was read in.\n", PARFIX[PPERPLAN*i+j]);
          exit(0);
        }
      } 
    }
    fgets(buffer, 1000, inputf);
  }
  for (i=0; i<PSTAR; i++) {
    fscanf(inputf, "%i", &PARFIX[NPL*PPERPLAN+i]);
    if (PARFIX[NPL*PPERPLAN+i] == 1) {
      nfixed++;
    } else {
      if (PARFIX[NPL*PPERPLAN+i] != 0) {
        printf("Error: All parfix values must be 0 or 1\n");
        printf("A value of %i was read in.\n", PARFIX[NPL*PPERPLAN+i]);
        exit(0);
      }
    } 
  } 
  fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &SPLITINCO); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &BIMODF); fgets(buffer, 1000, inputf); 
  printf("bimodf = %i\n", BIMODF);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  for (i=0; i<npl; i++) {
    for (j=0; j<PPERPLAN; j++) fscanf(inputf, "%i", &BIMODLIST[PPERPLAN*i+j]); fgets(buffer, 1000, inputf);
  }
  for (i=0; i<PSTAR; i++) fscanf(inputf, "%i", &BIMODLIST[NPL*PPERPLAN+i]); fgets(buffer, 1000, inputf);
  // Infrequently Edited Parameters
  for (i=0; i<4; i++) fgets(buffer, 1000, inputf); 
  fscanf(inputf, "%s %s %lf", type, varname, &T0); fgets(buffer, 1000, inputf); 
  printf("t0 = %lf\n", T0);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &T1); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  printf("t1 = %lf\n", T1);
  if (T1 <= T0 || EPOCH < T0 || EPOCH > T1) {
    printf("t0 = %lf, epoch=%lf, t1=%lf\n", T0, EPOCH, T1);
    printf("Error: You must have t0 <= epoch <= t1 and t0 < t1.\n");
    exit(0);
  } 

  fscanf(inputf, "%s %s %lu", type, varname, &OUTPUTINTERVAL); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  if (OUTPUTINTERVAL <= 0) {
    printf("Outputinterval must be a positive integer\n");
    printf("The current value is: %li\n", OUTPUTINTERVAL);
    printf("If that is not the value you specified, make sure your .in file is formatted correctly.\n");
    exit(0);
  }
  fscanf(inputf, "%s %s %lu", type, varname, &SEED); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &DISPERSE); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &OPTIMAL); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &RELAX); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &NPERBIN); fgets(buffer, 1000, inputf);
  if (CADENCESWITCH > 0 && NPERBIN <= 0) {
    printf("NPERBIN = %i\n", NPERBIN);
    printf("Error: Nperbin must be >= 1\n");
    exit(0);
  } 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &BINWIDTH); fgets(buffer, 1000, inputf); 
  if (CADENCESWITCH > 0 && BINWIDTH <= 0) {
    printf("BINWIDTH = %lf\n", BINWIDTH);
    printf("Error: binwidth must be > 0\n");
    exit(0);
  } 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &LTE); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &OFFSETMULT); fgets(buffer, 1000, inputf); 
  fscanf(inputf, "%s %s %lf", type, varname, &OFFSETMIN); fgets(buffer, 1000, inputf); 
  fscanf(inputf, "%s %s %lf", type, varname, &OFFSETMINOUT); fgets(buffer, 1000, inputf); 
  fscanf(inputf, "%s %s %lf", type, varname, &DIST2DIVISOR); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %lf", type, varname, &PRINTEPOCH); fgets(buffer, 1000, inputf);
  if (PRINTEPOCH > T1 || PRINTEPOCH < EPOCH) { 
    printf("Warning: We must have: epoch <= printepoch <= t1\n");
    printf("         printepoch was outside this range so printepoch was set to epoch (this should not be a problem in most cases)\n"); 
    PRINTEPOCH = EPOCH;
  }
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s %i", type, varname, &XYZFLAG); fgets(buffer, 1000, inputf);
  if (XYZFLAG < 0 || XYZFLAG > 6) {
    printf("XYZFLAG = %i\n", XYZFLAG);
    printf("Error: XYZFLAG must be in {0, 1, 2, 3, 4, 5, 6}\n");
    exit(0);
  }
  fgets(buffer, 1000, inputf);
  if (XYZFLAG==1 || XYZFLAG==2 || XYZFLAG==3) {
    XYZLIST=malloc(npl*sofi);
    fscanf(inputf, "%s %s", type, varname);
    for (j=0; j<npl; j++) fscanf(inputf, "%i", &XYZLIST[j]);
    fgets(buffer, 1000, inputf);
  } else {
    XYZLIST = calloc(npl, sofi);
    fgets(buffer, 1000, inputf);
  }
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s", type, varname); for (j=0; j<PPERPLAN; j++) fscanf(inputf, "%lf", &PSTEP[j]); fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fgets(buffer, 1000, inputf);
  fscanf(inputf, "%s %s", type, varname); for (i=0; i<PSTAR; i++) fscanf(inputf, "%lf", &SSTEP[i]); fgets(buffer, 1000, inputf); 
  fgets(buffer, 1000, inputf);

  fscanf(inputf, "%s %s %i", type, varname, &STEPFLAG); fgets(buffer, 1000, inputf); 
  if (STEPFLAG) {
    for (i=0; i<9; i++) fgets(buffer, 1000, inputf);
    for (i=0; i<npl; i++) {
      for (j=0; j<PPERPLAN; j++) fscanf(inputf, "%lf", &CSTEP[PPERPLAN*i+j]); fgets(buffer, 1000, inputf);
    }
    for (i=0; i<PSTAR; i++) fscanf(inputf, "%lf", &CSTEP[NPL*PPERPLAN+i]);
    STEP = CSTEP;
  } else {
    STEP = malloc((PPERPLAN*NPL+PSTAR)*sofd);
    free(CSTEP);
    for (i=0; i<npl; i++) memcpy(&STEP[i*PPERPLAN], &PSTEP[0], PPERPLAN*sofd);
    memcpy(&STEP[NPL*PPERPLAN], &SSTEP[0], PSTAR*sofd);
  }
  
  for (j=0; j<PPERPLAN*NPL+PSTAR; j++) {
    if (STEP[j] <= 0. && PARFIX[j] != 1) {
      printf("Warning: STEP[%i] = %lf \n", j, STEP[j]);
      printf("Setting STEP <= 0 may cause the code to crash or hang if PARFIX != 1\n");
    }
  }

  PPERWALKER = NPL*PPERPLAN + PSTAR;
  int nfree = PPERWALKER-nfixed;
#if (demcmc_compile == 1)
  if (NWALKERS <= nfree) {
    printf("Warning: N_walkers = %i, N_free = %i\n", NWALKERS, nfree); 
    printf("         N_Walkers should be > N free parameters.\n"); 
    printf("         The DEMCMC algorithm may fail to perform correctly.\n");
  } 
#endif
  fclose(inputf);
  return 0;
}


// Compute mean longitude from f, e, pomega
double getlambda(double f, double e, double pomega) {
  double bige = 2*atan(sqrt((1-e)/(1+e)) * tan(f/2));
  double lambda = pomega + bige - e*sin(bige);
  return lambda;
}


// Compute true anomaly from lambda, e, pomega
double getf(double lambda, double e, double pomega) {
  double bigm = lambda-pomega;
  double bige = bigm;
  int i;
  for (i=0; i<=25; i++) {
    bige = bigm + e*sin(bige);
  }
  double sendback = 2.0*atan( sqrt(1.0+e)/sqrt(1.0-e) * tan(bige/2.0) );
  return sendback;
}


// Modulo pi function
double *pushpi(double x[], long int lenx) {

  long int cyclecount = 0;
  long int i;
  for (i=0; i<lenx;i++){
    x[i] += 2.0*M_PI*cyclecount;
    while (x[i] > M_PI) {
      x[i] -= 2.0*M_PI;
      cyclecount -= 1;
    }
    while (x[i] < -M_PI) {
      x[i] += 2.0*M_PI;
      cyclecount += 1;
    }
  }
 
  return x;
}


// Get planetary coordinates in the i, omega, Omega reference frame
double *rotatekep(double x, double y, double i, double lo, double bo) {

  const int sofd = SOFD;

  double z=0.;
  double x1 = cos(lo)*x-sin(lo)*y;
  double y1 = sin(lo)*x+cos(lo)*y;
  double z1 = z;
  double x2 = x1;
  double y2 = cos(i)*y1-sin(i)*z1;
  double z2 = sin(i)*y1+cos(i)*z1;
  double x3 = cos(bo)*x2-sin(bo)*y2;
  double y3 = sin(bo)*x2+cos(bo)*y2;
  double z3 = z2;

  double *arr3 = malloc(3*sofd);
  arr3[0]=x3;arr3[1]=y3;arr3[2]=z3;
  return arr3;
}


// compute x, y, z, vx, vy, vz array from a, e, i, omega, Omega, f, mass array
double *keptostate(double a, double e, double i, double lo, double bo, double f, double m) {

  const int sofd = SOFD;

  double mass = m;
  if (isnan(m)) {
    printf("mass error\n");
    exit(0);
  }
  double r = a*(1.0-pow(e,2))/(1.0+e*cos(f));
  double x0 = r*cos(f);
  double y0 = r*sin(f);
  double vx0 = -sqrt(mass/(a*(1.0-pow(e,2)))) * sin(f);
  double vy0 = sqrt(mass/(a*(1.0-pow(e,2)))) * (e+cos(f));

  double *statearr = malloc(6*sofd);  
  double *state1arr, *state2arr;
  state1arr = rotatekep(x0, y0, i, lo, bo);
  state2arr = rotatekep(vx0, vy0, i, lo, bo);
  memcpy(statearr,state1arr,3*sofd);
  memcpy(statearr+3,state2arr,3*sofd);
  free(state1arr);free(state2arr);

  return statearr;

}


// Compute  a, e, i, omega, Omega, f array from x, y, zy, vx, vy, vz, mass array
double *statetokep (double x, double y, double z, double vx, double vy, double vz, double m) {

  double mass = G*m;

  double rsq = x*x + y*y + z*z;
  double r = sqrt(rsq);
  double vsq = vx*vx + vy*vy + vz*vz;
  double hx = y*vz - z*vy;
  double hy = z*vx - x*vz;
  double hz = x*vy - y*vx;
  double hsq = hx*hx + hy*hy + hz*hz;
  double h = sqrt(hsq);

  double rdot;
  if (vsq <= hsq/rsq) rdot=0.0;
  else rdot = sqrt( vsq - hsq/rsq );
  double rrdot = x*vx + y*vy + z*vz;
  if (rrdot < 0.0) rdot *= -1.0;

  double a = 1.0/(2.0/r - vsq/mass);
  double e, e1;
  e1 = 1.0 - hsq/(a*mass); // This little check is due to rounding errors when e=0.0 exactly causing negative square roots. 
  if (e1 <= 0) e=0.0;
  else e=sqrt(e1);
  double i = acos(hz/h);

  double sini = sin(i);
  double bo, lof;
  if (dbleq(sini, 0.0)) {
    bo = 0.0;
    lof = atan2(y/r, x/r);
  } else {
    double sinbo = hx / (h*sini);
    double cosbo = -hy / (h*sini);
    bo = atan2(sinbo, cosbo);

    double sinlof = z/(r*sini);
    double coslof = (1.0/cos(bo)) * (x/r + sin(bo)*sinlof*hz/h);
    lof = atan2(sinlof, coslof);
  }

  double lo, f;
  if (dbleq(e, 0.0)) {
    lo = 0.0;
    f = lof;
  } else {
    double sinf = a*(1.0-e*e)/(h*e)*rdot;
    double cosf = (1.0/e)*(a*(1.0-e*e)/r - 1.0);
    f = atan2(sinf, cosf);
    double lofmf = lof-f;
    double *lotemp = pushpi(&lofmf, 1);
    lo = lotemp[0];

  }


  double sofd = SOFD;
  double *kepelements = malloc(6*sofd); 
  kepelements[0] = a;
  kepelements[1] = e;
  kepelements[2] = i;
  kepelements[3] = lo;
  kepelements[4] = bo;
  kepelements[5] = f;

  return kepelements;

}


// Compute P, T0, e, i, Omega, omega array from a, e, i, omega, Omega, f, mass array
double *keptoorb(double a, double e, double i, double lo, double bo, double f, double mass) {

  double period = sqrt( (pow(a, 3)*4.0*M_PI*M_PI) / (G * mass) );
  double Tepoch = EPOCH;
  double E0 = 2.0*atan(sqrt( (1.0-e)/(1.0+e) ) * tan((M_PI/2.0-lo)/2.0) ); 
  double Ef = 2.0*atan(sqrt( (1.0-e)/(1.0+e) ) * tan(f/2.0) );
  double t0 = period/(2.0*M_PI) * ( 2.0*M_PI/period*Tepoch - Ef + E0 + e*sin(Ef) - e*sin(E0) );
  if (t0 < Tepoch) {
    do t0 += period;
    while (t0 < Tepoch);
  } else {
    while (t0 >= (Tepoch+period)) t0 -= period;
  }

  int sofd = SOFD;
  double *orbelements = malloc(6*sofd);
  orbelements[0] = period;
  orbelements[1] = t0;
  orbelements[2] = e;
  orbelements[3] = i;
  orbelements[4] = bo;
  orbelements[5] = lo;

  return orbelements;

}


// Compute Eccentric Anomaly from M and e
double getE(double M, double e) {
  int i;
  double bige=M;
  for (i=0; i<=25; i++) {
    bige = M + e*sin(bige);
  }
  return bige;
}


// Compute P, T0, e, i, Omega, omega vector from a, e, i, omega, Omega, MeanAnomoly, mass vector
double *keptoorbmean(double a, double e, double i, double lo, double bo, double M, double mass) {

  double period = sqrt( (pow(a, 3)*4.0*M_PI*M_PI) / (G * mass) );
  double Tepoch = EPOCH;
  double E0 = 2.0*atan(sqrt( (1.0-e)/(1.0+e) ) * tan((M_PI/2.0-lo)/2.0) ); 
  double Ef = getE(M, e); 
  double t0 = period/(2.0*M_PI) * ( 2.0*M_PI/period*Tepoch - Ef + E0 + e*sin(Ef) - e*sin(E0) );
  if (t0 < Tepoch) {
    do t0 += period;
    while (t0 < Tepoch);
  } else {
    while (t0 >= (Tepoch+period)) t0 -= period;
  }

  int sofd = SOFD;
  double *orbelements = malloc(6*sofd);
  orbelements[0] = period;
  orbelements[1] = t0;
  orbelements[2] = e;
  orbelements[3] = i;
  orbelements[4] = bo;
  orbelements[5] = lo;

  return orbelements;

}


// Compute Orbital xyzvxvyvz array from inputted orbital elements vector
double ***dsetup2 (double *p, const int npl){
  const int sofd = SOFD;
  const int sofds = SOFDS; 
 
  int pperplan = PPERPLAN;
  int pstar = PSTAR;

  double epoch = EPOCH;

  double brightstar=0;
  double bsum=0;
  int i;

  double ms = p[npl*pperplan+0];
  double rstar = p[npl*pperplan+1];
  double c1 = p[npl*pperplan+2];
  double c2 = p[npl*pperplan+3];
  double dilute = p[npl*pperplan+4];
#if (demcmc_compile == 0)
  if (ms <= 0.) {
    printf("Warning: Mstar <= 0\n", i);
    printf("This may cause a crash or hang.\n");
  }
  if (rstar <= 0.) {
    printf("Warning: Rstar <= 0\n", i);
    printf("This may cause a crash or hang.\n");
  }
  if (dilute < 0. || dilute > 1.) {
    printf("Warning: dilute value is not 0<=dilute<=1\n", i);
    printf("This may cause a crash or hang.\n");
  }
#endif

  double bigg = 1.0e0; //Newton's constant
  double ghere = G; //2.9591220363e-4; 
  double jos = 1.0/MSOMJ;  //9.545e-4; //M_jup/M_sol

  double *mp = malloc(npl*sofd);
  double *mpjup = malloc(npl*sofd);
  double *msys = malloc((npl+1)*sofd);
  msys[0] = ms;  

  double *a = malloc(npl*sofd);
  double *e = malloc(npl*sofd);
  double *inc = malloc(npl*sofd);
  double *bo = malloc(npl*sofd); 
  double *lo = malloc(npl*sofd);
  double *lambda = malloc(npl*sofd);
  double *f = malloc(npl*sofd);   

  for (i=0;i<npl; i++) {
    if (SQRTE) {
      e[i] = pow( sqrt(pow(p[i*pperplan+2],2)+pow(p[i*pperplan+3],2)), 2);
    } else {
      e[i] = sqrt(pow(p[i*pperplan+2],2)+pow(p[i*pperplan+3],2));
    }
    inc[i] = p[i*pperplan+4]*M_PI/180;
    bo[i] = p[i*pperplan+5]*M_PI/180;
    lo[i] = atan2(p[i*pperplan+3],p[i*pperplan+2]);
    mp[i] = p[i*pperplan+6];

#if (demcmc_compile == 0)
    if (mp[i] < 0.) {
      printf("Warning: mass[%i] < 0\n", i);
    }
    if (e[i] > 1.) {
      printf("Warning: e[%i] > 1\n", i);
      printf("This may cause a crash or hang.\n");
    }
#endif
 
    mpjup[i] = mp[i]*jos;       //          ; M_Jup
    msys[i+1] = msys[i]+mpjup[i];
    a[i] = cbrt(ghere*(msys[i+1])) * pow(cbrt(p[i*pperplan+0]),2) * pow(cbrt(2*M_PI),-2);

    double pomega = bo[i]+lo[i];
    double lambda0 = getlambda( (M_PI/2-lo[i]), e[i], pomega);
    double m0 = lambda0-pomega;
    double me = m0 + 2*M_PI*(epoch - p[i*pperplan+1])/p[i*pperplan+0]; 
    double mepomega = me+pomega;
    double *lambdaepoint = pushpi(&mepomega,1);
    double lambdae = lambdaepoint[0];
    f[i] = getf(lambdae, e[i], pomega);

  }

  double **state = malloc(npl*sofds);
  for (i=0; i<npl;i++) {
    state[i] = keptostate(a[i],e[i],inc[i],lo[i],bo[i],f[i],(ghere*msys[i+1]));
    int j;
    for (j=0; j<6; j++) {
      state[i][j] = -state[i][j];
    }
  }


  // jacobian
  double *sum;
  if (XYZFLAG==1) {
    for (i=0; i<npl; i++) {
      if(XYZLIST[i]) {
        int j;
        for (j=0; j<6; j++) state[i][j] = p[i*pperplan+j];
      }
    }
  }

#if (demcmc_compile == 0) 
  if (CONVERT && (XYZFLAG != 2) && (XYZFLAG != 3)) {
    char outfile2str[80];
    strcpy(outfile2str, "xyz_out_");
    strcat(outfile2str, OUTSTR);
    strcat(outfile2str, ".pldin");
    FILE *outfile2 = fopen(outfile2str, "a");
    fprintf(outfile2, "planet                 x                        y                      z                      v_x                   v_y                   v_z                       m                      rpors             ");
    fprintf(outfile2, "\n");
    double pnum = 0.1;
    for (i=0; i<NPL; i++) {
      fprintf(outfile2, "%1.1lf", pnum);
      int j;
      for (j=0; j<6; j++) {
        fprintf(outfile2, "\t%.15lf", state[i][j]);
      }
      fprintf(outfile2, "\t%.15lf", mpjup[i]/jos);
      fprintf(outfile2, "\t%.15lf", p[i*pperplan+7]);
      fprintf(outfile2, "\n");
      pnum+=0.1;
    }
    fprintf(outfile2, "%.15lf ; mstar\n", ms);
    fprintf(outfile2, "%.15lf ; rstar\n", rstar);
    fprintf(outfile2, "%.15lf ; c1\n", c1);
    fprintf(outfile2, "%.15lf ; c2\n", c2);
    fprintf(outfile2, "%.15lf ; dilute\n", dilute);
    if (RVJITTERFLAG) {
      for (i=0; i<RVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+i]);
      }
    }
    if (TTVJITTERFLAG) {
      for (i=0; i<TTVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+RVJITTERTOT+i]);
      }
    }
    if (CELERITE) {
      for (i=0; i<NCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+i]);
      }
    }
    if (RVCELERITE) {
      for (i=0; i<NRVCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; RV celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+i]);
      }
    }
    fprintf(outfile2, " ; These coordinates are Jacobian \n");

    fclose(outfile2);
  }
#endif


  
  double **statenew = malloc(npl*sofds);
  for (i=0; i<npl; i++) {
    statenew[i] = malloc(6*sofd);
  }

  memcpy(statenew[0], state[0], 6*sofd);
  int j;
  sum = calloc(6,sofd);
  for (i=1; i<npl; i++){
    int j;
    for (j=0; j<6; j++) {
      sum[j] += state[i-1][j]*mpjup[i-1]/msys[i];
      statenew[i][j] = state[i][j]+sum[j];
    }
  }
  free(sum);


  //barycentric
  if (XYZFLAG==3) {
    printf("barycentric coordinate input is broken at the moment. Try using stellar centric (xyzflag=2) instead.\n");
    exit(0);
    double *starpos = calloc(6,sofd);
    int j;
    for (j=0; j<6; j++) {
      for (i=0; i<npl; i++) {
        if (XYZLIST[i]) {
          starpos[j] -= p[i*pperplan+j]*mpjup[i];
        } else {
          starpos[j] -= state[i][j]*mpjup[i];
        }
      }
      starpos[j] /= ms;
    }
    for (i=0; i<npl; i++) {
      if (XYZLIST[i]) {
        int j;
        for (j=0; j<6; j++) {
          statenew[i][j] = p[i*pperplan+j]-starpos[j];
        }
      }
    }
    free(starpos);
  }

  //stellarcentric
  if (XYZFLAG==2) {
    for (i=0; i<npl; i++) {
      if (XYZLIST[i]) {
        int j;
        for (j=0; j<6; j++) {
          statenew[i][j] = p[i*pperplan+j];
        }
      }
    }
  }

#if (demcmc_compile == 0) 
  if (CONVERT) {
    char outfile2str[80];
    strcpy(outfile2str, "xyz_out_");
    strcat(outfile2str, OUTSTR);
    strcat(outfile2str, ".pldin");
    FILE *outfile2 = fopen(outfile2str, "a");
    fprintf(outfile2, "planet                 x                        y                      z                      v_x                   v_y                   v_z                       m                      rpors             ");
    fprintf(outfile2, "\n");
    double pnum = 0.1;
    for (i=0; i<NPL; i++) {
      fprintf(outfile2, "%1.1lf", pnum);
      int j;
      for (j=0; j<6; j++) {
        fprintf(outfile2, "\t%.15lf", statenew[i][j]);
      }
      fprintf(outfile2, "\t%.15lf", mpjup[i]/jos);
      fprintf(outfile2, "\t%.15lf", p[i*pperplan+7]);
      fprintf(outfile2, "\n");
      pnum+=0.1;
    }
    fprintf(outfile2, "%.15lf ; mstar\n", ms);
    fprintf(outfile2, "%.15lf ; rstar\n", rstar);
    fprintf(outfile2, "%.15lf ; c1\n", c1);
    fprintf(outfile2, "%.15lf ; c2\n", c2);
    fprintf(outfile2, "%.15lf ; dilute\n", dilute);
    if (RVJITTERFLAG) {
      for (i=0; i<RVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+i]);
      }
    }
    if (TTVJITTERFLAG) {
      for (i=0; i<TTVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+RVJITTERTOT+i]);
      }
    }
    if (CELERITE) {
      for (i=0; i<NCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+i]);
      }
    }
    if (RVCELERITE) {
      for (i=0; i<NRVCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; RV celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+i]);
      }
    }
    fprintf(outfile2, " ; These coordinates are stellar centric\n");

    fclose(outfile2);
  }
  if (CONVERT) {

    double *bary = calloc(6,sofd);
    int j;
    for (j=0; j<6; j++) {
      double mtot=ms;
      for (i=0; i<npl; i++) {
        if (XYZLIST[i]) {
          bary[j] += p[i*pperplan+j]*mpjup[i];
        } else {
          bary[j] += statenew[i][j]*mpjup[i];
        }
        mtot+=mpjup[i];
      }
      bary[j] /= mtot;
    }

    char outfile2str[80];
    strcpy(outfile2str, "xyz_out_");
    strcat(outfile2str, OUTSTR);
    strcat(outfile2str, ".pldin");
    FILE *outfile2 = fopen(outfile2str, "a");
    fprintf(outfile2, "planet                 x                        y                      z                      v_x                   v_y                   v_z                       m                      rpors             ");
    fprintf(outfile2, "\n");
    double pnum = 0.0;

    fprintf(outfile2, "%1.1lf", pnum);
    for (j=0; j<6; j++) {
      fprintf(outfile2, "\t%.15lf", -bary[j]);
    }
    fprintf(outfile2, "\t%.15lf", ms/jos);
    fprintf(outfile2, "\t%.15lf", 1.0);
    fprintf(outfile2, "\t%.15lf", brightstar);
    fprintf(outfile2, "\t%.15lf", c1);
    fprintf(outfile2, "\t%.15lf\n", c2);

    pnum+=0.1;
    for (i=0; i<NPL; i++) {
      fprintf(outfile2, "%1.1lf", pnum);
      int j;
      for (j=0; j<6; j++) {
        fprintf(outfile2, "\t%.15lf", statenew[i][j]-bary[j]);
      }
      fprintf(outfile2, "\t%.15lf", mpjup[i]/jos);
      fprintf(outfile2, "\t%.15lf", p[i*pperplan+7]);
      fprintf(outfile2, "\n");
      pnum+=0.1;
    }
    fprintf(outfile2, "%.15lf ; mstar\n", ms);
    fprintf(outfile2, "%.15lf ; rstar\n", rstar);
    fprintf(outfile2, "%.15lf ; c1\n", c1);
    fprintf(outfile2, "%.15lf ; c2\n", c2);
    fprintf(outfile2, "%.15lf ; dilute\n", dilute);
    if (RVJITTERFLAG) {
      for (i=0; i<RVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+i]);
      }
    }
    if (TTVJITTERFLAG) {
      for (i=0; i<TTVJITTERTOT; i++) {
        fprintf(outfile2, "%.15lf ; rvjitter \n", p[npl*pperplan+5+RVJITTERTOT+i]);
      }
    }
    if (CELERITE) {
      for (i=0; i<NCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+i]);
      }
    }
    if (RVCELERITE) {
      for (i=0; i<NRVCELERITE; i++) {
        fprintf(outfile2, "%.15lf ; RV celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+i]);
      }
    }
    fprintf(outfile2, " ; These coordinates are barycentric\n");

    free(bary);
    fclose(outfile2);
  }
#endif

  double ***integration_in = malloc(4*sizeof(double**));
  integration_in[0] = malloc(sofds);
  integration_in[0][0] = malloc(sofd);
  integration_in[0][0][0] = npl+1;
  integration_in[1] = malloc(sofds);
  integration_in[1][0] = malloc(sofd);
  integration_in[1][0][0] = epoch;
  integration_in[2] = malloc((npl+1)*sofds);
  for (i=0; i<npl+1; i++) {
    integration_in[2][i] = malloc(pperplan*sofd);
  }
  integration_in[2][0][0] = bigg*ms;
  integration_in[2][0][1] = 0.;
  integration_in[2][0][2] = 0.;
  integration_in[2][0][3] = 0.;
  integration_in[2][0][4] = 0.;
  integration_in[2][0][5] = 0.;
  integration_in[2][0][6] = 0.;
  integration_in[2][0][7] = 1.;
  for (i=0; i<npl; i++) {
    integration_in[2][i+1][0] = bigg*mpjup[i];
    int j;
    for (j=0; j<6; j++) {
      integration_in[2][i+1][j+1] = statenew[i][j];
    }
    integration_in[2][i+1][7] = p[i*pperplan+7]; 
  }
  integration_in[3] = malloc(sofds);
  integration_in[3][0] = malloc(4*sofd);
  integration_in[3][0][0] = rstar;
  integration_in[3][0][1] = c1;
  integration_in[3][0][2] = c2;
  integration_in[3][0][3] = dilute;




  free(mp);free(mpjup);free(msys);free(a);free(e);free(inc);free(bo);free(lo);free(lambda);free(f);
  for (i=0; i<npl; i++) {
    free(state[i]);free(statenew[i]);
  }
  free(state);free(statenew);

  return integration_in;

}




// Computes (observation - theory)/error for a [time, obs, theory, err] vector
double *devoerr (double **tmte) {

  const int sofd = SOFD;

  long ntimes = (long) tmte[0][0];
  double *resid = malloc((ntimes+1)*sofd);

  resid[0] = (double) ntimes;
  long i;
  for (i=0; i<ntimes; i++) {
    resid[i+1] = (tmte[1][i+1]-tmte[2][i+1])/tmte[3][i+1];
  }

  return resid;

}



////// These routines are from Pal 2008
//// from icirc.c
////
typedef struct
{	double	x0,y0;
	double	r;
} circle;

typedef struct
{	int	cidx;
	double	phi0,dphi;
	int	noidx;
	int	*oidxs;
} arc;

int icirc_arclist_intersections(circle *circles,int ncircle,arc **routs,int *rnout)
{
    int	i,j;
    
    arc	*arcs,*aouts;
    int	*acnt;
    int	naout;
    
    arcs=(arc *)malloc(sizeof(arc)*ncircle*ncircle*2);
    acnt=(int *)malloc(sizeof(int)*ncircle);
    
    for ( i=0 ; i<ncircle ; i++ )
    {	acnt[i]=0;		}
    
    for ( i=0 ; i<ncircle ; i++ )
    { for ( j=0 ; j<ncircle ; j++ )
    {	double	xa,ya,xb,yb,ra,rb;
        double	dx,dy,d;
        double	w,phia,phi0;
        
        if ( i==j )
            continue;
        
        xa=circles[i].x0;
        ya=circles[i].y0;
        ra=circles[i].r;
        xb=circles[j].x0;
        yb=circles[j].y0;
        rb=circles[j].r;
        dx=xb-xa;
        dy=yb-ya;
        d=sqrt(dx*dx+dy*dy);
        if ( ra+rb <= d )
            continue;
        else if ( d+ra <= rb )
            continue;
        else if ( d+rb <= ra )
            continue;
        w=(ra*ra+d*d-rb*rb)/(2*ra*d);
        if ( ! ( -1.0 <= w && w <= 1.0 ) )
            continue;
        phia=acos(w);
        
        phi0=atan2(dy,dx);
        if ( phi0 < 0.0 )	phi0+=2*M_PI;
		
        if ( acnt[i] <= 0 )
        {	arc	*a;
            a=&arcs[2*ncircle*i];
            a[0].phi0=phi0-phia;
            a[0].dphi=2*phia;
            a[1].phi0=phi0+phia;
            a[1].dphi=2*(M_PI-phia);
            acnt[i]=2;
        }
        else
        {	arc	*a;
            double	wp[2],w,dw;
            int	k,n,l;
            wp[0]=phi0-phia;
            wp[1]=phi0+phia;
            a=&arcs[2*ncircle*i];
            n=acnt[i];
            for ( k=0 ; k<2 ; k++ )
            {	w=wp[k];
                for ( l=0 ; l<n ; l++ )
                {	dw=w-a[l].phi0;
                    while ( dw<0.0 )	dw+=2*M_PI;
                    while ( 2*M_PI<=dw )	dw-=2*M_PI;
                    if ( dw<a[l].dphi )
                        break;
                }
                if ( l<n )
                {	memmove(a+l+1,a+l,sizeof(arc)*(n-l));
                    a[l+1].phi0=a[l].phi0+dw;
                    a[l+1].dphi=a[l].dphi-dw;
                    a[l].dphi=dw;
                    n++;
                }
            }
            acnt[i]=n;
        }
        
    }
    }
    
    naout=0;
    for ( i=0 ; i<ncircle ; i++ )
    {	if ( acnt[i] <= 0 )
		naout++;
	else
		naout+=acnt[i];
    }
    aouts=(arc *)malloc(sizeof(arc)*naout);
    j=0;
    for ( i=0 ; i<ncircle ; i++ )
    {	int	k;
        if ( acnt[i] <= 0 )
        {	aouts[j].cidx=i;
            aouts[j].phi0=0.0;
            aouts[j].dphi=2*M_PI;
            j++;
        }
        else
        {	for ( k=0 ; k<acnt[i] ; k++ )
        {	aouts[j].cidx=i;
			aouts[j].phi0=arcs[2*ncircle*i+k].phi0;
			aouts[j].dphi=arcs[2*ncircle*i+k].dphi;
			j++;
        }
        }
    }
    for ( j=0 ; j<naout ; j++ )
    {	double	x,y,dx,dy;
        int	k;
        
        i=aouts[j].cidx;
        if ( acnt[i] <= 0 )
        {	x=circles[i].x0+circles[i].r;
            y=circles[i].y0;
        }
        else
        {	double	phi;
            phi=aouts[j].phi0+0.5*aouts[j].dphi;
            x=circles[i].x0+circles[i].r*cos(phi);
            y=circles[i].y0+circles[i].r*sin(phi);
        }
        aouts[j].noidx=0;
        aouts[j].oidxs=NULL;
        for ( k=0 ; k<ncircle ; k++ )
        {	if ( i==k )	continue;
            dx=x-circles[k].x0;
            dy=y-circles[k].y0;
            if ( dx*dx+dy*dy < circles[k].r*circles[k].r )
            {	aouts[j].oidxs=(int *)realloc(aouts[j].oidxs,sizeof(int)*(aouts[j].noidx+1));
                *(aouts[j].oidxs+aouts[j].noidx)=k;
                aouts[j].noidx++;
            }
        }
    }
    
    if ( routs != NULL )	*routs=aouts;
    if ( rnout != NULL )	*rnout=naout;
    
    free(acnt);
    free(arcs);
    
    return(0);
}

int icirc_arclist_free(arc *arcs,int narc)
{
    int	i;
    for ( i=0 ; i<narc ; i++ )
    {	if ( arcs[i].oidxs != NULL )
		free(arcs[i].oidxs);
    }
    free(arcs);
    return(0);
}


// From elliptic.c
//
//

#define		FMIN(a,b)	((a)<(b)?(a):(b))
#define		FMAX(a,b)	((a)<(b)?(a):(b))
#define		SQR(a)		((a)*(a))

#define C1 0.3
#define C2 (1.0/7.0)
#define C3 0.375
#define C4 (9.0/22.0)

double carlson_elliptic_rc(double x,double y)
{
    double alamb,ave,s,w,xt,yt,ans;
    
    if ( y > 0.0 )
    {	xt=x;
        yt=y;
        w=1.0;
    }
    else
    {	xt=x-y;
        yt = -y;
        w=sqrt(x)/sqrt(xt);
    }
    do
    {	alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        ave=(1.0/3.0)*(xt+yt+yt);
        s=(yt-ave)/ave;
    } while ( fabs(s) > 0.0012 );
    
    ans=w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
    
    return(ans);
}

#undef	C4
#undef	C3
#undef	C2
#undef	C1

double carlson_elliptic_rf(double x,double y,double z)
{
    double	alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;
    xt=x;
    yt=y;
    zt=z;
    do
    {	sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=(1.0/3.0)*(xt+yt+zt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
    } while ( fabs(delx) > 0.0025 || fabs(dely) > 0.0025 || fabs(delz) > 0.0025 );
    e2=delx*dely-delz*delz;
    e3=delx*dely*delz;
    return((1.0+((1.0/24.0)*e2-(0.1)-(3.0/44.0)*e3)*e2+(1.0/14.0)*e3)/sqrt(ave));
}

#define C1 (3.0/14.0)
#define C2 (1.0/6.0)
#define C3 (9.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.25*C3)
#define C6 (1.5*C4)

double carlson_elliptic_rd(double x,double y,double z)
{
    double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,
	sqrtx,sqrty,sqrtz,sum,xt,yt,zt,ans;
    
    xt=x;
    yt=y;
    zt=z;
    sum=0.0;
    fac=1.0;
    do
    {	sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        sum+=fac/(sqrtz*(zt+alamb));
        fac=0.25*fac;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=0.2*(xt+yt+3.0*zt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
    } while ( fabs(delx) > 0.0015 || fabs(dely) > 0.0015 || fabs(delz) > 0.0015 );
    ea=delx*dely;
    eb=delz*delz;
    ec=ea-eb;
    ed=ea-6.0*eb;
    ee=ed+ec+ec;
    ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
                     +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
    return(ans);
}

#undef	C6
#undef	C5
#undef	C4
#undef	C3
#undef	C2
#undef	C1

#define C1 (3.0/14.0)
#define C2 (1.0/3.0)
#define C3 (3.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.75*C3)
#define C6 (1.5*C4)
#define C7 (0.5*C2)
#define C8 (C3+C3)

double carlson_elliptic_rj(double x,double y,double z,double p)
{
    double	a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
	ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
    
    sum=0.0;
    fac=1.0;
    if ( p > 0.0 )
    {	xt=x;
        yt=y;
        zt=z;
        pt=p;
        a=b=rcx=0.0;
    }
    else
    {	xt=FMIN(FMIN(x,y),z);
        zt=FMAX(FMAX(x,y),z);
        yt=x+y+z-xt-zt;
        a=1.0/(yt-p);
        b=a*(zt-yt)*(yt-xt);
        pt=yt+b;
        rho=xt*zt/yt;
        tau=p*pt/yt;
        rcx=carlson_elliptic_rc(rho,tau);
    }
    do
    {	sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
        beta=pt*SQR(pt+alamb);
        sum += fac*carlson_elliptic_rc(alpha,beta);
        fac=0.25*fac;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        pt=0.25*(pt+alamb);
        ave=0.2*(xt+yt+zt+pt+pt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
        delp=(ave-pt)/ave;
    } while ( fabs(delx)>0.0015 || fabs(dely)>0.0015 || fabs(delz)>0.0015 || fabs(delp)>0.0015 );
    ea=delx*(dely+delz)+dely*delz;
    eb=delx*dely*delz;
    ec=delp*delp;
    ed=ea-3.0*ec;
    ee=eb+2.0*delp*(ea-ec);
    
    ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
                     +delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
    
    if ( p <= 0.0 ) ans=a*(b*ans+3.0*(rcx-carlson_elliptic_rf(xt,yt,zt)));
    
    return(ans);
}

#undef	C6
#undef	C5
#undef	C4
#undef	C3
#undef	C2
#undef	C1
#undef	C8
#undef	C7

#undef			SQR
#undef			FMAX
#undef			FMIN


// From mttr.c
//
//

// This function returns F[phi'] (Eq 34) from r=r, rho=c, x=phi'
double mttr_integral_primitive(double r,double c,double x)
{
    double	q2,s2,d2,sx,cx,w;
    double	rf,rd,rj;
    double	beta;
    double	iret;
    int         d2neq0, cneqr;
    
    // Eq 24-26 with r=r, rho=c, x=phi'
    q2=r*r+c*c+2*r*c*cos(x);
    d2=r*r+c*c-2*r*c;
    s2=r*r+c*c+2*r*c;


    double epsilon=3.0*DBL_EPSILON;
    d2neq0 = !((d2 < epsilon) && ( -d2 < epsilon));
    cneqr = !(((c-r) < epsilon) && (-(c-r) < epsilon));

    sx=sin(x/2);
    cx=cos(x/2);
    
    if ( 1.0<q2 )	q2=1.0;
    
    w=(1-q2)/(1-d2);
    if ( w<0.0 )	w=0.0;
    // Eq 31-33
    rf=carlson_elliptic_rf(w,sx*sx,1);
    rd=carlson_elliptic_rd(w,sx*sx,1);
    if ( d2neq0 && cneqr)	rj=carlson_elliptic_rj(w,sx*sx,1,q2/d2);
    else		        rj=0.0;
    
    // Eq 34 line 1
    beta=atan2((c-r)*sx,(c+r)*cx);
    iret=-beta/3;
    // Eq 34 line 2 first term
    iret+=x/6;
    
    w=cx/sqrt(1-d2);
    
    // Eq (34) lines 2-6
    iret+=
    +2.0/ 9.0*c*r*sin(x)*sqrt(1-q2)
    +1.0/ 3.0*(1+2*r*r*r*r-4*r*r)*w*rf
    +2.0/ 9.0*r*c*(4-7*r*r-c*c+5*r*c)*w*rf
    -4.0/27.0*r*c*(4-7*r*r-c*c)*w*cx*cx*rd;
    if ( d2neq0 && cneqr )
    //if ( d2neq0 )
        iret += 1.0/3.0*w*(r+c)/(r-c)*(rf-(q2-d2)/(3*d2)*rj);
    else
        iret -= 1.0/3.0*w*(r+c)*(q2-d2)*M_PI/(2*q2*sqrt(q2));

    return(iret);
}

double mttr_integral_definite(double r,double c,double x0,double dx)
{
    double	dc,nx;
    double	ret;
    
    
    // Eq (22) or (21) ?
    if ( c<=0.0 )
    {	if ( r<1.0 )
		ret=(1-(1-r*r)*sqrt(1-r*r))*dx/3.0;
	else /* this case implies r==1: */
		ret=dx/3.0;
        
        return(ret);
    }
    
    if ( dx<0.0 )
    {	x0+=dx;
        dx=-dx;
    }
    while ( x0<0.0 )	x0+=2*M_PI;
    while ( 2*M_PI<=x0 )	x0-=2*M_PI;
    
    ret=0.0;
    while ( 0.0<dx )
    {	dc=2*M_PI-x0;
        if ( dx<dc )	dc=dx,nx=x0+dx;
        else		nx=0.0;
        
        ret+=mttr_integral_primitive(r,c,x0+dc)-mttr_integral_primitive(r,c,x0);
        
        x0=nx;
        dx-=dc;
    }
   
    return(ret);
}

/*****************************************************************************/

// This function takes n circles, computes their overlap, and returns
// The c's are coefficients of the stellar flux terms
//c1 is constant term 
//c2 is polynomial term
double mttr_flux_general(circle *circles,int ncircle,double c0,double c1,double c2)
{
    arc	*arcs,*a;
    int	i,narc;
    double	fc,f0;
    
    // Get circle intersections into *arcs
    icirc_arclist_intersections(circles,ncircle,&arcs,&narc);
    
    fc=0.0;
    
    for ( i=0 ; i<narc ; i++ )
    {	double	sign,x0,y0,r,p0,dp,p1,df0,df1,df2;
        double	x2,y2,r2;
        
        a=&arcs[i];
        if ( a->cidx==0 && a->noidx<=0 )
            sign=+1;
        else if ( a->cidx != 0 && a->noidx==1 && a->oidxs[0]==0 )
            sign=-1;
        else
            continue;
        
        x0=circles[a->cidx].x0;
        y0=circles[a->cidx].y0;
        r =circles[a->cidx].r;
        p0=a->phi0;
        dp=a->dphi;
        p1=p0+dp;
        
        x2=x0*x0;
        y2=y0*y0;
        r2=r*r;
        
        // Eq 8 last line
        df0=0.5*r*(x0*(sin(p1)-sin(p0))+y0*(-cos(p1)+cos(p0)))+0.5*r*r*dp;
        
        // If c1, then constant term
        if ( c1 != 0.0 )
        {	double	delta,rho;
            delta=atan2(y0,x0);
            rho=sqrt(x2+y2);
            // this avoids edge cases where rho ~= r and the integral doesn't work properly
            if ( fabs(rho) > 1e-7 && fabs(r) > 1e-7  && fabs(rho-r) < 1e-7)  rho += 2e-7;
            df1=mttr_integral_definite(r,rho,p0-delta,dp);
        }
        else
            df1=0.0;
        
        // If c2 then polynomial term
        if ( c2 != 0.0 )
            // Eq (18)
        {	df2=(r/48.0)*(	+(24*(x2+y2)+12*r2)*r*dp
                          -4*y0*(6*x2+2*y2+9*r2)*(cos(p1)-cos(p0))
                          -24*x0*y0*r*(cos(2*p1)-cos(2*p0))
                          -4*y0*r2*(cos(3*p1)-cos(3*p0))
                          +4*x0*(2*x2+6*y2+9*r2)*(sin(p1)-sin(p0))
                          -4*x0*r2*(sin(3*p1)-sin(3*p0))
                          -r2*r*(sin(4*p1)-sin(4*p0)) );
        }
        else
            df2=0.0;
        
        fc += sign*(c0*df0+c1*df1+c2*df2);
        
        
    }
    
    // normalize
    f0=2.0*M_PI*(c0/2.0+c1/3.0+c2/4.0);
    
    icirc_arclist_free(arcs,narc);
    
    if ( 0.0<f0 )
        return(fc/f0);
    else
        return(0.0);
}


//// These routines are from Maxsted 2018

double clip(double a, double b, double c)
{
        if (a < b) 
                return b;
        else if (a > c) 
                return c;
        else 
                return a;
}

double q1(double z, double p, double c, double a, double g, double I_0)
{
        double zt = clip(z, 0,1-p);
        double s = 1-zt*zt;
        double c0 = (1-c+c*pow(s,g));
        double c2 = 0.5*a*c*pow(s,(g-2))*((a-1)*zt*zt-1);
        return 1-I_0*M_PI*p*p*(c0 + 0.25*p*p*c2 - 0.125*a*c*p*p*pow(s,(g-1)));
}

double q2(double z, double p, double c, double a, double g, double I_0, double eps)
{
        double zt = clip(z, 1-p,1+p);
        double d = clip((zt*zt - p*p + 1)/(2*zt),0,1);
        double ra = 0.5*(zt-p+d);
        double rb = 0.5*(1+d);
        double sa = clip(1-ra*ra,eps,1);
        double sb = clip(1-rb*rb,eps,1);
        double q = clip((zt-d)/p,-1,1);
        double w2 = p*p-(d-zt)*(d-zt);
        double w = sqrt(clip(w2,eps,1));
        double c0 = 1 - c + c*pow(sa,g);
        double c1 = -a*c*ra*pow(sa,(g-1));
        double c2 = 0.5*a*c*pow(sa,(g-2))*((a-1)*ra*ra-1);
        double a0 = c0 + c1*(zt-ra) + c2*(zt-ra)*(zt-ra);
        double a1 = c1+2*c2*(zt-ra);
        double aq = acos(q);
        double J1 =  (a0*(d-zt)-(2./3.)*a1*w2 + 0.25*c2*(d-zt)*(2.0*(d-zt)*(d-zt)-p*p))*w + (a0*p*p + 0.25*c2*pow(p,4))*aq ;
        double J2 = a*c*pow(sa,(g-1))*pow(p,4)*(0.125*aq + (1./12.)*q*(q*q-2.5)*sqrt(clip(1-q*q,0.0,1.0)) );
        double d0 = 1 - c + c*pow(sb,g);
        double d1 = -a*c*rb*pow(sb,(g-1));
        double K1 = (d0-rb*d1)*acos(d) + ((rb*d+(2./3.)*(1-d*d))*d1 - d*d0)*sqrt(clip(1-d*d,0.0,1.0));
        double K2 = (1/3)*c*a*pow(sb,(g+0.5))*(1-d);
        return 1 - I_0*(J1 - J2 + K1 - K2);
}

double Flux_drop_analytical_power_2(double d_radius, double k, double c, double a, double f, double eps)
{
        /*
        Calculate the analytical flux drop por the power-2 law.
        
        Parameters
        d_radius : double
                Projected seperation of centers in units of stellar radii.
        k : double
                Ratio of the radii.
        c : double
                The first power-2 coefficient.
        a : double
                The second power-2 coefficient.
        f : double
                The flux from which to drop light from.
        eps : double
                Factor (1e-9)
        */
        double I_0 = (a+2)/(M_PI*(a-c*a+2));
        double g = 0.5*a;

        //else if (abs(d_radius-1) < k) return q2(d_radius, k, c, a, g, I_0, 1e-9);
        if (d_radius < 1-k) return q1(d_radius, k, c, a, g, I_0);
        if (fabs(d_radius-1) < k) return q2(d_radius, k, c, a, g, I_0, 1e-9);
        else return 1.0;
}




// Lightcurve output functions
//
// This function computes the lightcurve for a single time given the positions of the planets and properties of the star
double onetlc (int nplanets, circle* system, double rstar, double c1, double c2) {
        
    if (nplanets==0) return 1.0;
 
    double flux;

    // Pal+2011 mutual event law
    if (LDLAW == 0) {
    
        double c0 = 1.0;
        double g0, g1, g2;
        g0 = c0-c1-2.0*c2;
        g1 = c1+2.0*c2;
        g2 = c2;
    
        flux = mttr_flux_general(system, nplanets+1, g0, g1, g2);// /nflux; 

    // Maxsted+2018 power2 law
    } else if (LDLAW == 1) {

        double runningflux = 1.0;
        int i;
        for (i=0; i<nplanets; i++) {
            //double d_radius = sqrt(pow(rxy[i][1]/rstar, 2) + pow(rxy[i][2]/rstar, 2));
            double d_radius = sqrt(pow(system[i+1].x0, 2) + pow(system[i+1].y0, 2));
            double k = system[i+1].r;
            double c = c1;
            double a = c2;
            double eps = 1e-9;
            double f = 1.0; // This doesn't seem to do anything in their code

            double thisplanetflux = Flux_drop_analytical_power_2(d_radius, k, c, a, f, eps);
            double thisdrop = 1.0 - thisplanetflux;
            runningflux -= thisdrop;
        } 

        flux = runningflux;

    } else {
        printf("Error: Currently the only LDLAW options are 1 or 2\n");
        printf("       Please select one of these\n");
        exit(0);
    } 
     

    return flux;

}


// Computes the lightcurve for a list of times (and cadences) for an array of planet positions as a function of time and the stellar properties
double *timedlc ( double *times, int *cadences, long ntimes, double **transitarr, int nplanets, double rstar, double c1, double c2) {

    const int sofd = SOFD;
    const int cadenceswitch = CADENCESWITCH;

    if (ntimes==0) {
      return NULL;
    }

    double *fluxlist = malloc(ntimes*sofd);
    const double rsinau = RSUNAU; //0.0046491; // solar radius in au
    const double rstarau = rstar*rsinau;


    double xstart[nplanets];
    double ystart[nplanets];
    int i;
    for (i=0; i<nplanets; i++) {
       xstart[i] = (transitarr[i][1] - transitarr[i][3]*(transitarr[i][0]-times[0]))/(rstarau);
       ystart[i] = (transitarr[i][2] - transitarr[i][4]*(transitarr[i][0]-times[0]))/(rstarau);
    }

    circle sun = {0.,0., 1.};
    circle system[nplanets+1];
    system[0] = sun;
    for (i=0; i< nplanets; i++) {
        system[i+1].x0 = xstart[i];
        system[i+1].y0 = ystart[i];
        system[i+1].r = transitarr[i][5];
    }

 
    double flux;   
    double t_cur;
    double t_next;
    long n=0;
    if (cadenceswitch==2) {
      long j=0;
      long jj=0;
      for (n=0; n<ntimes-1; n++) {
        while (cadences[j]==0) j++;
        t_cur = times[j];
        jj=j+1;
        while (cadences[j]==0) jj++;
        t_next = times[n+1];
        flux = onetlc (nplanets, system, rstarau, c1, c2);
        fluxlist[n] = flux;
        for (i=0; i<nplanets; i++) {
            system[i+1].x0 += transitarr[i][3]*(t_next-t_cur)/(rstarau);
            system[i+1].y0 += transitarr[i][4]*(t_next-t_cur)/(rstarau);
        }
        j=jj;
      }
      fluxlist[ntimes-1] = onetlc (nplanets, system, rstarau, c1, c2);


    } else {
      for (n=0; n<ntimes-1; n++) {
        t_cur = times[n];
        t_next = times[n+1];
        flux = onetlc (nplanets, system, rstarau, c1, c2);
        fluxlist[n] = flux;
        for (i=0; i<nplanets; i++) {
            system[i+1].x0 += transitarr[i][3]*(t_next-t_cur)/(rstarau);
            system[i+1].y0 += transitarr[i][4]*(t_next-t_cur)/(rstarau);
        }
      }
      fluxlist[ntimes-1] = onetlc(nplanets, system, rstarau, c1, c2);
    }

    return fluxlist;

}


// Computes and bins the lightcurve for a list of times (and cadences) for an array of planet positions as a function of time and the stellar properties
double *binnedlc ( double *times, int *cadences, long ntimes, double binwidth, int nperbin,  double **transitarr, int nplanets, double rstar, double c1, double c2) {

    //need to take into account extra 1's on either side of transit if bin is wider than buffer.  
    const int sofd = SOFD;
    const int cadenceswitch = CADENCESWITCH;

    if (ntimes==0) {
      return NULL;
    }

    double *fluxlist = malloc(ntimes*sofd);
    double rsinau = RSUNAU; //0.0046491; // solar radius in au
    long maxcalls = ntimes*nperbin;// + 1;
    double rstarau = rsinau*rstar; 
    double *fulltimelist=malloc(maxcalls*sofd);
    long j = 0;
    long jj = 0;
    if (cadenceswitch==2) {
      for (j=0; j<ntimes; j++) {
        if (cadences[j]==1) {
          int k;
          for (k=0; k<nperbin; k++) {
            fulltimelist[jj]=times[j]-binwidth/2+binwidth/nperbin*k+binwidth/(2*nperbin);
            jj++;
          }
        } else {
          fulltimelist[jj] = times[j];
          jj++;
        }
      }
    } else {
      for(j=0; j<ntimes; j++) {
        int k;
        for (k=0; k<nperbin; k++) {
          fulltimelist[nperbin*j+k]=times[j]-binwidth/2+binwidth/nperbin*k+binwidth/(2*nperbin);
        }
      }
    }

    double xstart[nplanets];
    double ystart[nplanets];
    int i;
    for (i=0; i<nplanets; i++) {
       xstart[i] = (transitarr[i][1] - transitarr[i][3]*(transitarr[i][0]-fulltimelist[0]))/rstarau;
       ystart[i] = (transitarr[i][2] - transitarr[i][4]*(transitarr[i][0]-fulltimelist[0]))/rstarau;
    }

    circle sun = {0.,0., 1.};
    circle system[nplanets+1];
    system[0] = sun;
    for (i=0; i< nplanets; i++) {
        system[i+1].x0 = xstart[i];
        system[i+1].y0 = ystart[i];
        system[i+1].r = transitarr[i][5];
    }

 
    double flux;   
    long n=0;
    if (cadenceswitch==2) {
      int nlong=0;
      for (n=0; n<ntimes; n++) {
        int nn;
        double binnedflux=0;
        if (cadences[n] == 1) {
          for (nn=0; nn<nperbin; nn++) {
            binnedflux += onetlc (nplanets, system, rstarau, c1, c2);
            //if (nn < (nperbin-1)) {
              double t_cur = fulltimelist[(n-nlong) + nlong*nperbin + nn];
              double t_next = fulltimelist[(n-nlong) + nlong*nperbin + nn + 1];
              for (i=0; i<nplanets; i++) {
                system[i+1].x0 += transitarr[i][3]*(t_next-t_cur)/rstarau;
                system[i+1].y0 += transitarr[i][4]*(t_next-t_cur)/rstarau;
              }
            //}
          }
          binnedflux = binnedflux/nperbin;
          nlong++;
        } else {
          binnedflux += onetlc (nplanets, system, rstarau, c1, c2);
          //if (nn < (nperbin-1)) {
            double t_cur = fulltimelist[(n-nlong) + nlong*nperbin];
            double t_next = fulltimelist[(n-nlong) + nlong*nperbin + 1];
            for (i=0; i<nplanets; i++) {
              system[i+1].x0 += transitarr[i][3]*(t_next-t_cur)/rstarau;
              system[i+1].y0 += transitarr[i][4]*(t_next-t_cur)/rstarau;
            }
          //}
        }
        fluxlist[n] = binnedflux;
      }
    } else {
      for (n=0; n<ntimes; n++) {
        int nn;
        double binnedflux=0;
        for (nn=0; nn<nperbin; nn++) {
          binnedflux += onetlc (nplanets, system, rstarau, c1, c2);
          //if (nn < (nperbin-1)) {
            double t_cur = fulltimelist[n*nperbin+nn];
            double t_next = fulltimelist[n*nperbin+nn+1];
            for (i=0; i<nplanets; i++) {
              system[i+1].x0 += transitarr[i][3]*(t_next-t_cur)/rstarau;
              system[i+1].y0 += transitarr[i][4]*(t_next-t_cur)/rstarau;
            }
          //}
        }
        binnedflux = binnedflux/nperbin;
        fluxlist[n] = binnedflux;
      }
    }

    free(fulltimelist);
//exit(0); 
    return fluxlist;

}




// This runs the DEMCMC 
int demcmc(char aei[], char chainres[], char bsqres[], char gres[]) {

  // Load in global vars (This isn't really necessary but is convenient to ensure not changing global vars)
  const long nwalkers = NWALKERS;
  const int nsteps = NSTEPS;
  const int cadenceswitch = CADENCESWITCH;
  const char *tfilename = TFILE;
  const int *parfix = PARFIX;
  const unsigned long seed = SEED;
  const double disperse = DISPERSE;
  const double optimal = OPTIMAL;
  const double relax = RELAX;
  const double *step = STEP;
  const int bimodf = BIMODF;
  const int *bimodlist = BIMODLIST;
  const int pperplan = PPERPLAN;
  const int pstar = PSTAR;
  const int npl = NPL;
  const int nbodies = NBODIES;
  const int splitincO = SPLITINCO;

  const int sofd = SOFD;
  const int sofds = SOFDS;
  const int sofi = SOFI;
  const int sofis = SOFIS;

  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng *rnw;
  gsl_rng_env_setup();
  //T=gsl_rng_ranlxd2;
  T=gsl_rng_taus2;

  r=gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  rnw=gsl_rng_alloc(T);
  
  if (TTVCHISQ) {

    int i;
    FILE **ttvfiles=malloc(npl*sizeof(FILE*));
    for (i=0; i<npl; i++) {
      char tempchar[1000];
      sprintf(tempchar, "ttv_%02i.in", i+1);
      printf("ttvfilei = %s\n", tempchar);
      ttvfiles[i]=fopen(tempchar, "r");
    }

    long maxttvs = 1000;
    // NTTV[i] is the list of the transit index for each transit with a measured time for planet i
    // Like TTTV, ETTV, and MTTV, the first entry is the total number of data points
    NTTV=malloc(npl*sizeof(long*));
    for (i=0; i<npl; i++) NTTV[i] = malloc(maxttvs*sizeof(long));
    TTTV=malloc(npl*sofds);
    for (i=0; i<npl; i++) TTTV[i] = malloc(maxttvs*sofd);
    ETTV=malloc(npl*sofds);
    for (i=0; i<npl; i++) ETTV[i] = malloc(maxttvs*sofd);
#if (demcmc_compile==0)
    MTTV=malloc(npl*sofds);
    for (i=0; i<npl; i++) MTTV[i] = malloc(maxttvs*sofd);
#endif

    for (i=0; i<npl; i++) {
      if (ttvfiles[i] == NULL) {
        printf("Bad ttv file name");
        exit(0);
      }
      long tt=0;
      while (fscanf(ttvfiles[i], "%li %lf %lf", &NTTV[i][tt+1], &TTTV[i][tt+1], &ETTV[i][tt+1]) == 3) { 
        if (tt>=maxttvs-1) {
          FILE *sout = fopen("demcmc.stdout", "a");
          fprintf(sout, "Too many TTVs, adjust maxttvs or use correct file\n");
          fclose(sout);
          exit(0);
        }
        printf("ttvread in planet %i num %li index %i time %lf err %lf\n", i, tt, NTTV[i][tt+1], TTTV[i][tt+1], ETTV[i][tt+1]);
        tt++;
      }
      NTTV[i][0]=tt;
      TTTV[i][0]=tt;
      ETTV[i][0]=tt;
      NTTV[i][tt+1]=LONG_MAX;
      TTTV[i][tt+1]=HUGE_VAL;
      ETTV[i][tt+1]=HUGE_VAL;
    }
 
    for (i=0; i<npl; i++) {
      fclose(ttvfiles[i]);
    }
    free(ttvfiles);
  }

  FILE *sout;
  double *p = malloc((npl*pperplan + pstar)*sofd);  

  if (!RESTART) {
 
    int i;
    double *planet1 = malloc(npl*sofd);
    double *period1 = malloc(npl*sofd);
    double *t01 = malloc(npl*sofd);
    double *e1 = malloc(npl*sofd);
    double *inc1 = malloc(npl*sofd);
    double *bo1 = malloc(npl*sofd);
    double *lo1 = malloc(npl*sofd);
    double *mp1 = malloc(npl*sofd);
    double *rpors1 = malloc(npl*sofd);
    double *brightness1;
    double *c1bin1;
    double *c2bin1;
    double *celeriteps = malloc(4*sofd);
    double *rvceleriteps = malloc(4*sofd);
    double ms;
    double c0;
    double c1;
    double rstar; 
    double dilute;
    double *jittersize;
    if (RVJITTERFLAG) {
      jittersize=malloc(RVJITTERTOT*sofd);
    }
    double *jittersizettv;
    if (TTVJITTERFLAG) {
      jittersizettv=malloc(TTVJITTERTOT*sofd);
    }
 
    char buffer[1000];
  
    FILE *aeifile = fopen(aei, "r");
    if (aeifile == NULL) {
      printf("Bad pldin Input File Name");
      exit(0);
    }
  
    fgets(buffer, 1000, aeifile);
    for (i=0; i<npl; i++) {
      fscanf(aeifile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &planet1[i], &period1[i], &t01[i], &e1[i], &inc1[i], &bo1[i], &lo1[i], &mp1[i], &rpors1[i]);
      printf("%lf\n", period1[i]);
      p[i*pperplan+0] = period1[i];
      p[i*pperplan+1] = t01[i];
      if (XYZFLAG==4 || XYZFLAG==5 || XYZFLAG==6) {
        p[i*pperplan+2] = e1[i];
        p[i*pperplan+3] = inc1[i];
        p[i*pperplan+4] = bo1[i];
        p[i*pperplan+5] = lo1[i];
      } else {
        if (XYZLIST[i]) {
          p[i*pperplan+2] = e1[i];
          p[i*pperplan+3] = inc1[i];
          p[i*pperplan+4] = bo1[i];
          p[i*pperplan+5] = lo1[i];
        } else {
          if (SQRTE) {
            p[i*pperplan+2] = sqrt(e1[i]) * cos(lo1[i]*M_PI/180);
            p[i*pperplan+3] = sqrt(e1[i]) * sin(lo1[i]*M_PI/180);
          } else {
            p[i*pperplan+2] = e1[i] * cos(lo1[i]*M_PI/180);
            p[i*pperplan+3] = e1[i] * sin(lo1[i]*M_PI/180);
          }
          p[i*pperplan+4] = inc1[i];
          p[i*pperplan+5] = bo1[i];
        }
      }
      p[i*pperplan+6] = mp1[i];
      p[i*pperplan+7] = rpors1[i];
      fgets(buffer, 1000, aeifile); // This line usually unnecessarily advances the file pointer to a new line. 
                                    // Unless  you are not multistar and but have too many inputs. It then  saves you from bad read-ins 
    }
    fscanf(aeifile, "%lf", &ms);
    fgets(buffer, 1000, aeifile);
    fscanf(aeifile, "%lf", &rstar);
    fgets(buffer, 1000, aeifile);
    fscanf(aeifile, "%lf", &c0);
    fgets(buffer, 1000, aeifile);
    fscanf(aeifile, "%lf", &c1);
    fgets(buffer, 1000, aeifile);
    fscanf(aeifile, "%lf", &dilute);
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<RVJITTERTOT; ki++) {
        fgets(buffer, 1000, aeifile);
        fscanf(aeifile, "%lf", &jittersize[ki]);
      }
    }
    if (TTVJITTERFLAG) {
      int ki;
      for (ki=0; ki<TTVJITTERTOT; ki++) {
        fgets(buffer, 1000, aeifile);
        fscanf(aeifile, "%lf", &jittersizettv[ki]);
      }
    }
    if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        fgets(buffer, 1000, aeifile);
        fscanf(aeifile, "%lf", &celeriteps[ki]);
      }
    }
    if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        fgets(buffer, 1000, aeifile);
        fscanf(aeifile, "%lf", &rvceleriteps[ki]);
      }
    }
 
    p[npl*pperplan+0] = ms;
    p[npl*pperplan+1] = rstar;
    p[npl*pperplan+2] = c0;
    p[npl*pperplan+3] = c1;
    p[npl*pperplan+4] = dilute; 
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<RVJITTERTOT; ki++) {
        p[npl*pperplan+5+ki] = jittersize[ki];
      }
    }
    if (TTVJITTERFLAG) {
      int ki;
      for (ki=0; ki<TTVJITTERTOT; ki++) {
        p[npl*pperplan+5+RVJITTERTOT+ki] = jittersizettv[ki];
      }
    }
    if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+ki] = celeriteps[ki];
      }
    }
    if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+ki] = rvceleriteps[ki];
      }
    }

    // if you gave me input in aei instead of PT0e, change it:
    if (XYZFLAG==5) {
      double masstot[npl+1]; 
      masstot[0] = ms;
      for (i=0; i<npl; i++) {
        masstot[i+1] = masstot[i];
        masstot[i+1] += p[i*pperplan+6]/MSOMJ;
      }
      for (i=0; i<npl; i++) {
        double *orbelements = keptoorb(p[i*pperplan+0], p[i*pperplan+1], p[i*pperplan+2]*M_PI/180., p[i*pperplan+3]*M_PI/180., p[i*pperplan+4]*M_PI/180., p[i*pperplan+5]*M_PI/180., masstot[i+1]);
        p[i*pperplan+0] = orbelements[0];
        p[i*pperplan+1] = orbelements[1];
        if (SQRTE) {
          p[i*pperplan+2] = sqrt(orbelements[2]) * cos(orbelements[5]);
          p[i*pperplan+3] = sqrt(orbelements[2]) * sin(orbelements[5]);
        } else {
          p[i*pperplan+2] = orbelements[2] * cos(orbelements[5]);
          p[i*pperplan+3] = orbelements[2] * sin(orbelements[5]);
        }
        p[i*pperplan+4] = orbelements[3]*180./M_PI;
        p[i*pperplan+5] = orbelements[4]*180./M_PI;

        free(orbelements);
      }
    }
    // aei but with mean anomaly not true 
    if (XYZFLAG==6) {
      double masstot[npl+1]; 
      masstot[0] = ms;
      for (i=0; i<npl; i++) {
        masstot[i+1] = masstot[i];
        masstot[i+1] += p[i*pperplan+6]/MSOMJ;
      }
      for (i=0; i<npl; i++) {
        double *orbelements = keptoorbmean(p[i*pperplan+0], p[i*pperplan+1], p[i*pperplan+2]*M_PI/180., p[i*pperplan+3]*M_PI/180., p[i*pperplan+4]*M_PI/180., p[i*pperplan+5]*M_PI/180., masstot[i+1]);
        p[i*pperplan+0] = orbelements[0];
        p[i*pperplan+1] = orbelements[1];
        if (SQRTE) {
          p[i*pperplan+2] = sqrt(orbelements[2]) * cos(orbelements[5]);
          p[i*pperplan+3] = sqrt(orbelements[2]) * sin(orbelements[5]);
        } else {
          p[i*pperplan+2] = orbelements[2] * cos(orbelements[5]);
          p[i*pperplan+3] = orbelements[2] * sin(orbelements[5]);
        }
        p[i*pperplan+4] = orbelements[3]*180./M_PI;
        p[i*pperplan+5] = orbelements[4]*180./M_PI;
        free(orbelements);
      }
    }

    fclose(aeifile);

    free(planet1); free(period1); free(t01); free (e1); free(inc1); free(bo1); free(lo1); free(mp1); free(rpors1); free(celeriteps); free(rvceleriteps); 
    if (RVJITTERFLAG) {
      free(jittersize);
    }
    if (TTVJITTERFLAG) {
      free(jittersizettv);
    }

  }

  int nparam = pperplan*npl+pstar;

  // read in list of times
  double *timelist;
  double *fluxlist;
  double *errlist;
  int *cadencelist;
  long timelistlen = 2000000;
  timelist = malloc(timelistlen*sofd);
  errlist = malloc(timelistlen*sofd);
  fluxlist = malloc(timelistlen*sofd);
  if (fluxlist==NULL) {
    sout = fopen("demcmc.stdout", "a");
    fprintf(sout, "malloc error\n");
    fclose(sout);
    exit(0);
  }

  FILE *tfile = fopen(tfilename, "r");

  printf ("%s\n", tfilename);

  if (tfile == NULL) {
    sout = fopen("demcmc.stdout", "a");
    printf("Error opening tfile\n");
    fprintf(sout,"Error opening tfile\n");
    fclose(sout);
    exit(0);
  }
  long num;
  double flux1, err1;
  long kk=0;

  if (CADENCESWITCH==0 || CADENCESWITCH==1) {
    while (fscanf(tfile, "%ld %lf %lf %lf %lf %lf", &num, &timelist[kk+1], &flux1, &err1, &fluxlist[kk+1], &errlist[kk+1]) == 6) { 
      if (kk>=timelistlen-1) {
        timelistlen+=1000000;
        timelist = realloc(timelist, timelistlen*sofd);
        fluxlist = realloc(fluxlist, timelistlen*sofd);
        errlist = realloc(errlist, timelistlen*sofd);
        if (timelist==NULL) {
          sout = fopen("demcmc.stdout", "a");
          fprintf(sout, "timelist allocation failure\n");
          fclose(sout);
          exit(0);
        }
      }
      kk++;
    }
  } else {
    cadencelist = malloc(timelistlen*sofi);
    while (fscanf(tfile, "%ld %lf %lf %lf %lf %lf %i", &num, &timelist[kk+1], &flux1, &err1, &fluxlist[kk+1], &errlist[kk+1], &cadencelist[kk]) == 7) { 
      if (kk>=timelistlen-1) {
        timelistlen+=1000000;
        timelist = realloc(timelist, timelistlen*sofd);
        fluxlist = realloc(fluxlist, timelistlen*sofd);
        cadencelist = realloc(cadencelist, timelistlen*sofi);
        errlist = realloc(errlist, timelistlen*sofd);
        
        if (timelist==NULL) {
          sout = fopen("demcmc.stdout", "a");
          fprintf(sout, "timelist allocation failure\n");
          fclose(sout);
          exit(0);
        }
      }
      kk++;
    }
    printf ("\n%s read\n", tfilename);
  }
  fclose(tfile);
  
  timelist[0]=kk;
  fluxlist[0]=kk;
  errlist[0]=kk;

  double **tfe = malloc(3*sofds);
  tfe[0]=timelist;
  tfe[1]=fluxlist;
  tfe[2]=errlist;

  double *rvtimelist;
  double *rvlist;
  double *rverrlist;
  double *rvbodylist;
  int *rvtelelist;
  double *rvtelelistd;
  long rvlistlen = 5000; // Maximum number of RVs by default
  long vv=0;
  double **tve = malloc(5*sofds);

  printf("pre RV\n");
  if (RVS) {
    rvtimelist = malloc(rvlistlen*sofd);
    rvlist = malloc(rvlistlen*sofd);
    rverrlist = malloc(rvlistlen*sofd);
    rvbodylist = malloc(rvlistlen*sofd);
    rvtelelist = malloc(rvlistlen*sofi);
    rvtelelistd = malloc(rvlistlen*sofd);
    if (rverrlist==NULL) {
      sout = fopen("demcmc.stdout", "a");
      fprintf(sout, "malloc error\n");
      fclose(sout);
      exit(0);
    }

    int i;
    for(i=0; i<nbodies; i++) {
      printf("%i, %i\n", i, RVARR[i]);
      if(RVARR[i]) { 
        FILE *rvfile = fopen(RVFARR[i], "r");
        if (rvfile == NULL) {
          sout = fopen("demcmc.stdout", "a");
          fprintf(sout,"Error opening rvfile\n");
          fclose(sout);
          exit(0);
        } 
        while (fscanf(rvfile, "%lf %lf %lf %i", &rvtimelist[vv+1], &rvlist[vv+1], &rverrlist[vv+1], &rvtelelist[vv+1]) == 4) { 
          rvbodylist[vv+1] = (double) i;
          if (vv>=rvlistlen-1) {
            rvlistlen+=1000;
            rvtimelist = realloc(rvtimelist, rvlistlen*sofd);
            rvlist = realloc(rvlist, rvlistlen*sofd);
            rverrlist = realloc(rverrlist, rvlistlen*sofd);
            rvbodylist = realloc(rvbodylist, rvlistlen*sofd);
            rvtelelist = realloc(rvtelelist, rvlistlen*sofi);
            rvtelelistd = realloc(rvtelelistd, rvlistlen*sofd); 
            if (rverrlist==NULL || rvtelelistd==NULL) {
              sout = fopen("demcmc.stdout", "a");
              fprintf(sout, "rvlist allocation failure\n");
              fclose(sout);
              exit(0);
            }
          }
          vv++;
        }
        fclose(rvfile);
      }
    }
  
    printf("vv=%li\n",vv);
  
    rvtimelist[0] = (double) vv;
    rvlist[0] = (double) vv;
    rverrlist[0] = (double) vv;
    rvbodylist[0] = (double) vv;
    rvtelelist[0] = (int) vv;
    rvtelelistd[0] = (double) vv;

    long w;
    for (w=1; w<(vv+1); w++) {
      rvlist[w] = rvlist[w] * MPSTOAUPD;
      rverrlist[w] = rverrlist[w] * MPSTOAUPD;
    }
 
    double bigrvlist[vv*5];
    long z;
    for (z=1; z<(vv+1); z++) {
      bigrvlist[(z-1)*5+0] = rvtimelist[z];
      bigrvlist[(z-1)*5+1] = rvlist[z];
      bigrvlist[(z-1)*5+2] = rverrlist[z];
      bigrvlist[(z-1)*5+3] = rvbodylist[z];
      bigrvlist[(z-1)*5+4] = (double) rvtelelist[z];
    }
    qsort(bigrvlist, vv, 5*sofd, compare);
    for (z=1; z<(vv+1); z++) {
      rvtimelist[z] = bigrvlist[(z-1)*5+0];
      rvlist[z] = bigrvlist[(z-1)*5+1];
      rverrlist[z] = bigrvlist[(z-1)*5+2];
      rvbodylist[z] = bigrvlist[(z-1)*5+3]; 
      rvtelelistd[z] = bigrvlist[(z-1)*5+4];
    }

    tve[0]=rvtimelist;
    tve[1]=rvlist;
    tve[2]=rverrlist;
    tve[3]=rvbodylist;
    tve[4]=rvtelelistd;
    free(rvtelelist);    
  }

  long ttvlistlen = 5000; // Max number TTVs by default
  long vvt=0;
  double **nte;
  if (TTVCHISQ) {
    nte = malloc(3*sofds);
    int kij;
    for (kij=0; kij<3; kij++) {
      nte[kij] = malloc(ttvlistlen*sofd);
    }
    int i;
    int ki;
    int ksofar=0;
    for (i=0; i<npl; i++) {
      for (ki=0; ki<NTTV[i][0]; ki++) {
        nte[0][1+ki+ksofar] = (double) NTTV[i][1+ki]; 
        nte[1][1+ki+ksofar] = TTTV[i][1+ki];
        nte[2][1+ki+ksofar] = ETTV[i][1+ki];
      }
      ksofar += NTTV[i][0];
    }
    nte[0][0] = (double) ksofar;
    nte[1][0] = (double) ksofar;
    nte[2][0] = (double) ksofar;
  }

#if (demcmc_compile==1)
  int i;
  double ***p0 = malloc(nwalkers*sizeof(double**));
  for (i=0; i<nwalkers; i++){
    p0[i] = malloc((npl+1)*sofds);
    int j;
    for (j=0; j<npl; j++) {
      p0[i][j] = malloc(pperplan*sofd); 
    }
    p0[i][npl] = malloc(pstar*sofd);
  }

  int w;

  double gamma;
  double neg2loglikemin; 
  double *gammaN;

  int fixedpars = 0;
  for (i=0; i<nparam; i++) fixedpars += parfix[i];
  int ndim = nparam - fixedpars;
  if (ndim == 0) {
    printf("Warning! No free parameters!\n");
    ndim=1;
  }

  long jj = 0;


  if (RESTART) {
  
      FILE *restartf = fopen(chainres, "r");
      if (restartf == NULL) {
        printf("Error: %s\n", chainres); 
        printf("nofile\n");
        exit(0);
      }
      double ignore; 
      char ignorec[1000];
      for (i=0; i<nwalkers; i++) {
        int j; 
        for (j=0; j<npl; j++) {
          fscanf(restartf, " %lf %lf %lf %lf %lf %lf %lf %lf %lf", &ignore, &p0[i][j][0], &p0[i][j][1], &p0[i][j][2], &p0[i][j][3], &p0[i][j][4], &p0[i][j][5], &p0[i][j][6], &p0[i][j][7]);
        }
        fscanf(restartf, "%lf", &p0[i][npl][0]);
        fgets(ignorec, 1000, restartf);
        fscanf(restartf, "%lf", &p0[i][npl][1]);
        fgets(ignorec, 1000, restartf);
        fscanf(restartf, "%lf", &p0[i][npl][2]);
        fgets(ignorec, 1000, restartf);
        fscanf(restartf, "%lf", &p0[i][npl][3]);
        fgets(ignorec, 1000, restartf);
        fscanf(restartf, "%lf", &p0[i][npl][4]);
        fgets(ignorec, 1000, restartf);
        if (RVJITTERFLAG) {
          int ki;
          for (ki=0; ki<(RVJITTERTOT); ki++) {
            fscanf(restartf, "%lf", &p0[i][npl][5+ki]);
            fgets(ignorec, 1000, restartf);
          }
        }
        if (TTVJITTERFLAG) {
          int ki;
          for (ki=0; ki<(TTVJITTERTOT); ki++) {
            fscanf(restartf, "%lf", &p0[i][npl][5+RVJITTERTOT+ki]);
            fgets(ignorec, 1000, restartf);
          }
        }
        if (CELERITE) {
          int ki;
          for (ki=0; ki<NCELERITE; ki++) {
            fscanf(restartf, "%lf", &p0[i][npl][5+RVJITTERTOT+TTVJITTERTOT+ki]);
            fgets(ignorec, 1000, restartf);
          }
        }
        if (RVCELERITE) {
          int ki;
          for (ki=0; ki<NRVCELERITE; ki++) {
            fscanf(restartf, "%lf", &p0[i][npl][5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+ki]);
            fgets(ignorec, 1000, restartf);
          }
        }
        fgets(ignorec, 1000, restartf);
      }
      fclose(restartf);
      printf("Read in Restart pldin\n");
     
      FILE *regammaf = fopen(gres,"r");
      long genres;
      double grac, gammares;
      fscanf(regammaf, "%li %lf %lf", &genres, &grac, &gammares); 
      fclose(regammaf);
     
      printf("Read in Restart gamma\n");
      // rounding generation
      genres = genres/10;
      genres = genres*10;
      jj = genres+1;
      // optimal multiplier
      gamma = gammares;
     
      FILE *rebsqf = fopen(bsqres, "r");
      char garbage[10000];
      for (i=0; i<(1+npl+pstar+1); i++) {
        fgets(garbage, 10000, rebsqf);
      }
      char cgarbage;
      for (i=0; i<11; i++) {
        cgarbage = fgetc(rebsqf);
      }
      double reneg2loglikemin;
      fscanf(rebsqf, "%lf", &reneg2loglikemin);
      fclose(rebsqf);
      neg2loglikemin = reneg2loglikemin;
     
      printf("Read in Restart best chi sq\n");
      printf("neg2loglikemin=%lf\n", neg2loglikemin);
      printf("gamma=%lf\n", gamma); 
    
  
  } else {  //if not RESTART

    int j;
    for (j=0; j<npl; j++) {
      if ( ((IGT90==1) && p[j*pperplan+4] < 90.) || ((IGT90==2) && p[j*pperplan+4] > 90.) ) {
        printf("Warning: at least one planet has initial inclination that is not allowed\n");
        printf("         Compare the igt90 flag in the .in file to the planetary i values in the .pldin file\n");
        printf("         This may cause initialization to hang.\n");
        break;
      } 
    }

    printf("Initializing Walkers\n");
    // take small random steps to initialize walkers
    for (i=0; i<nwalkers; i++) {
      int j;
      for (j=0; j<npl; j++) {
        int k;
        for (k=0; k<pperplan; k++) {
          // make sure all inclinations still > 90.0
          do {
            double epsilon = (1-parfix[j*pperplan+k])*step[j*pperplan+k]*gsl_ran_gaussian(r, 1.0);
            // don't split inclinations....
            //if (splitincO && k==4 && ((i+nwalkers/2)/nwalkers)) p0[i][j][k] = 180. - p[j][k] + epsilon; 
            //else if (splitincO && k==5 && ((i+nwalkers/2)/nwalkers)) p0[i][j][k] = -p[j][k] + epsilon;
            //else p0[i][j][k] = p[j][k] + epsilon;
            p0[i][j][k] = p[j*pperplan+k] + epsilon;
            if ( (int) ceil(disperse) ) p0[i][j][k] += disperse*(1-parfix[j*pperplan+k])*step[j*pperplan+k]*gsl_ran_gaussian(r, 1.0);
          } while ( ((IGT90==1) && (k==4 && p0[i][j][k] < 90.0)) || ((IGT90==2) && (k==4 && p0[i][j][k] > 90.0)) || (MGT0 && (k==6 && p0[i][j][k] < 0.0)) );
        }
      }
      for (j=0; j<pstar; j++) {
        do {
        p0[i][npl][j] = p[npl*pperplan+j] + (1-parfix[npl*pperplan+j])*step[npl*pperplan+j]*gsl_ran_gaussian(r, 1.0);
        if ( (int) ceil(disperse) ) p0[i][npl][j] += disperse*(1-parfix[npl*pperplan+j])*step[npl*pperplan+j]*gsl_ran_gaussian(r, 1.0);
        } while ( (TTVJITTERFLAG || RVJITTERFLAG) && (j > 4 && p0[i][npl][j] < 0.0) ); 
      }
    } 
 
    // optimal multiplier
    gamma = 2.38 / sqrt(2.*ndim);

  }

  printf("Initialized Walkers\n");

  int pperwalker = PPERWALKER; //npl*pperplan+pstar;
  int totalparams = nwalkers*pperwalker;
  double *p0local = malloc(totalparams*sofd);
  double **p0localN; 

  for (i=0; i<nwalkers; i++) {
    int j;
    for (j=0; j<npl; j++) {
      int k;
      for (k=0; k<pperplan; k++) {
        p0local[i*pperwalker + j*pperplan + k] = p0[i][j][k];
      }
    }
    for (j=0; j<pstar; j++) {
      p0local[i*pperwalker + npl*pperplan + j] = p0[i][npl][j];
    }
  }

#endif


  double ***dsetup2 (double *p, const int npl); //prototype 
  double ***dpintegrator_single (double ***int_in, double **tfe, double **tve, double **nte, int *cadencelist); //prototype 
  int dpintegrator_single_megno (double ***int_in); //prototype 
  double *devoerr (double **tmte); //prototype

#if (demcmc_compile==0)
  printf("Sanity check:\n");
  printf("p[][2]=%lf, orbelemets[2]=%lf\n", p[0], p[0]);

  double ***int_in = dsetup2(p, npl);
  printf("int_in %lf, %lf, %lf, %lf\n", int_in[0][0][0], int_in[1][0][0], int_in[2][0][0], int_in[2][1][0]);

  double ***flux_rvs; 
  flux_rvs = dpintegrator_single(int_in, tfe, tve, nte, cadencelist);
  printf("Integration Complete\n");
  
  double **ttvts = flux_rvs[2];
  double **flux = flux_rvs[0];
  double **radvs = flux_rvs[1];
  double *dev = devoerr(flux);
  double neg2loglike = 0;
  long il;
  long maxil = (long) dev[0];
  printf("kk=%li\n", maxil);
 
  if (! CELERITE) { 
    for (il=0; il<maxil; il++) neg2loglike += dev[il+1]*dev[il+1];
  } else { // if celerite
    neg2loglike = celerite_fit(flux_rvs, p, 0, 0, 1);
  }

  printf("pre rvjitter:\n");
  printf("neg2loglike=%lf\n", neg2loglike);

  if (RVS) {
    if (! RVCELERITE) { 
      double *newelist;
      if (RVJITTERFLAG) {
        long kj;
        long maxkj = (long) tve[0][0]; 
        newelist=malloc((maxkj+1)*sofd);
        newelist[0] = (double) maxkj;
        for (kj=0; kj<maxkj; kj++) {
          int jitterindex = (int) tve[3][1+kj]*NTELESCOPES + tve[4][1+kj];
          double sigmajitter = p[npl*pperplan+5+jitterindex]*MPSTOAUPD;
          double quadsum = sigmajitter*sigmajitter + radvs[3][1+kj]*radvs[3][1+kj];
          // double check this... factor of 1/2
          neg2loglike += log(quadsum / (radvs[3][1+kj]*radvs[3][1+kj]) );
          newelist[1+kj] = sqrt( quadsum );
        }
        radvs[3] = newelist;
      }
      double *rvdev = devoerr(radvs);
      long maxil = (long) rvdev[0];
      for (il=0; il<maxil; il++) neg2loglike += rvdev[il+1]*rvdev[il+1];
      free(rvdev);
    } else { // if rvcelerite
      neg2loglike += celerite_fit(flux_rvs, p, 0, 1, 1);
    }
  }

  if (TTVCHISQ) {
    double *newelistttv;
    if (TTVJITTERFLAG) {
      long kj;
      long maxkj = (long) ttvts[0][0];
      newelistttv = malloc((maxkj+1)*sofd);
      newelistttv[0] = (double) maxkj;
      for (kj=0; kj<maxkj; kj++) {
        int jitterindex = (kj < NTTV[0][0]) ? 0 : 1 ; 
        double sigmajitter = p[npl*pperplan+5+RVJITTERTOT+jitterindex];
        double quadsum = sigmajitter*sigmajitter + ttvts[3][1+kj]*ttvts[3][1+kj];
        // check that the index on ttvts should be 2
        neg2loglike += log(quadsum / (ttvts[3][1+kj]*ttvts[3][1+kj]) );
        newelistttv[1+kj] = sqrt(quadsum);
      }
      ttvts[3] = &newelistttv[0];
    }
    printf("%lf %lf %lf %lf\n", ttvts[0][0], ttvts[1][0], ttvts[2][0], ttvts[3][0]);
    double *ttvdev = devoerr(ttvts);
    printf("ttvts\n");
    long maxil = (long) ttvdev[0];
    for (il=0; il<maxil; il++) neg2loglike += ttvdev[il+1]*ttvdev[il+1];
    free (ttvdev);
    if (TTVJITTERFLAG) {
      free(newelistttv);
    }
  }
  printf("post ttvjitter:\n");
  printf("neg2loglike=%lf\n", neg2loglike);
  
//  double photoradius = int_in[3][0][0]; 
//  if (SPECTROSCOPY) {
//    if (photoradius > SPECRADIUS) neg2loglike += pow( (photoradius - SPECRADIUS) / SPECERRPOS, 2 );
//    else neg2loglike += pow( (photoradius - SPECRADIUS) / SPECERRNEG, 2 );
//  }
//  double photomass = int_in[2][0][0]; 
//  if (MASSSPECTROSCOPY) {
//    if (photomass > SPECMASS) neg2loglike += pow( (photomass - SPECMASS) / MASSSPECERRPOS, 2 );
//    else neg2loglike += pow( (photomass - SPECMASS) / MASSSPECERRNEG, 2 );
//  }
//  printf("neg2loglikenoinc=%lf\n", neg2loglike);
//  if (INCPRIOR) {
//    int i0;
//    for (i0=0; i0<npl; i0++) {
//      printf("%lf\n", neg2loglike);
//      neg2loglike += -2.0*log( sin(p[i0*pperplan+4] *M_PI/180.) ); 
//    }
//  }
//  printf("neg2loglikenoe=%lf\n", neg2loglike);
//  double* evector = malloc(npl*sofd);
//  if (ECUTON || EPRIOR) {
//    if (SQRTE) {
//      int i0;
//      for (i0=0; i0<npl; i0++) {
//        evector[i0] = pow(sqrt( pow(p[i0*pperplan+2], 2) + pow(p[i0*pperplan+3], 2) ), 2);
//      }
//    } else {
//      int i0;
//      for (i0=0; i0<npl; i0++) {
//        evector[i0] = sqrt( pow(p[i0*pperplan+2], 2) + pow(p[i0*pperplan+3], 2) );
//      }
//    }
//  }
//  if (EPRIOR) {
//    int i0;
//    for (i0=0; i0<NPL; i0++) {
//      if (EPRIORV[i0]) {
//        double priorprob;
//        if (EPRIOR==1) {
//          priorprob = rayleighpdf(evector[i0]);
//        } else if (EPRIOR==2) {
//          priorprob = normalpdf(evector[i0]);
//        }
//        neg2loglike += -2.0*log( priorprob );
//      }
//    }
//  }

  if (XYZFLAG != 0) {
    // generate p in normal format 
    // compute priors
    int i, j;
    double **aeiparam  = malloc(npl*sofds);
    double **orbparam  = malloc(npl*sofds);
    double masstot[npl+1]; 
    double stateorig[npl][6];
    masstot[0] = int_in[2][0][0];
    for (i=0; i<npl; i++) {
      masstot[i+1] = masstot[i];
      masstot[i+1] += int_in[2][i+1][0];
    }

    for (j=0; j<6; j++) {
      stateorig[0][j] = -int_in[2][1][j+1];
    }
    double *sum = calloc(6,sofd);
    for (i=1; i<npl; i++){
      for (j=0; j<6; j++) {
        sum[j] += int_in[2][i][j+1]*int_in[2][i][0]/masstot[i];
        stateorig[i][j] = -(int_in[2][i+1][j+1] - sum[j]);
      }
    }
    free(sum);
    
    for (i=0; i<npl; i++) {
      aeiparam[i] = statetokep(stateorig[i][0], stateorig[i][1], stateorig[i][2], stateorig[i][3], stateorig[i][4], stateorig[i][5], masstot[i+1]); 
      orbparam[i] = keptoorb(aeiparam[i][0], aeiparam[i][1], aeiparam[i][2], aeiparam[i][3], aeiparam[i][4], aeiparam[i][5], masstot[i+1]);
    }
    double *ptemp = malloc((npl*pperplan + pstar)*sofd); 
    for (i=0; i<npl; i++) {
      ptemp[i*pperplan+0] = orbparam[i][0];
      ptemp[i*pperplan+1] = orbparam[i][1];
      if (SQRTE) {
        ptemp[i*pperplan+2] = sqrt(orbparam[i][2]) * cos(orbparam[i][5]);
        ptemp[i*pperplan+3] = sqrt(orbparam[i][2]) * sin(orbparam[i][5]);
      } else {
        ptemp[i*pperplan+2] = orbparam[i][2] * cos(orbparam[i][5]);
        ptemp[i*pperplan+3] = orbparam[i][2] * sin(orbparam[i][5]);
      }
      ptemp[i*pperplan+4] = orbparam[i][3];
      ptemp[i*pperplan+5] = orbparam[i][4];
      ptemp[i*pperplan+6] = p[i*pperplan+6];
      ptemp[i*pperplan+7] = p[i*pperplan+7];
    } 
    for (i=0; i<pstar; i++) {
      ptemp[npl*pperplan+i] = p[npl*pperplan+i];
    }      
    
    for (i=0; i<npl; i++) free(aeiparam[i]);
    free(aeiparam);
    for (i=0; i<npl; i++) free(orbparam[i]);
    free(orbparam);
 
    neg2loglike += compute_priors(ptemp, 0);
    free(ptemp);

  } else {
    neg2loglike += compute_priors(p, 0);
  }

  printf("post priors:\n");
  printf("neg2loglike=%lf\n", neg2loglike);
  free(dev);

  if (CONVERT) {  
    int i, j;
    double **aeiparam  = malloc(npl*sofds);
    double **orbparam  = malloc(npl*sofds);
    double masstot[npl+1]; 
    double stateorig[npl][6];
    masstot[0] = int_in[2][0][0];
    for (i=0; i<npl; i++) {
      masstot[i+1] = masstot[i];
      masstot[i+1] += int_in[2][i+1][0];
    }

    for (j=0; j<6; j++) {
      stateorig[0][j] = -int_in[2][1][j+1];
    }
    double *sum = calloc(6,sofd);
    for (i=1; i<npl; i++){
      for (j=0; j<6; j++) {
        sum[j] += int_in[2][i][j+1]*int_in[2][i][0]/masstot[i];
        stateorig[i][j] = -(int_in[2][i+1][j+1] - sum[j]);
      }
    }
    free(sum);
    
    for (i=0; i<npl; i++) {
      aeiparam[i] = statetokep(stateorig[i][0], stateorig[i][1], stateorig[i][2], stateorig[i][3], stateorig[i][4], stateorig[i][5], masstot[i+1]); 
      orbparam[i] = keptoorb(aeiparam[i][0], aeiparam[i][1], aeiparam[i][2], aeiparam[i][3], aeiparam[i][4], aeiparam[i][5], masstot[i+1]);
    }

    char outfile2str[80];
    strcpy(outfile2str, "aei_out_");
    strcat(outfile2str, OUTSTR);
    strcat(outfile2str, ".pldin");
    FILE *outfile2 = fopen(outfile2str, "a");
    fprintf(outfile2, "planet         period (d)               T0 (d)                  e                   i (deg)                 Omega (deg)               omega(deg)               mp (mjup)              rpors           ");
    fprintf(outfile2, "\n");
    char ch = 'a';
    double pnum = 0.1;
    int ip;
    for (ip=0; ip<npl; ip++) {
      fprintf(outfile2, "%1.1lf", pnum);
      ch++; 
      pnum+=0.1;
      int ii;
      fprintf(outfile2, "\t%.15lf", orbparam[ip][0]); 
      fprintf(outfile2, "\t%.15lf", orbparam[ip][1]);
      fprintf(outfile2, "\t%.15lf", orbparam[ip][2]);
      fprintf(outfile2, "\t%.15lf", orbparam[ip][3]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", orbparam[ip][4]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", orbparam[ip][5]*180.0/M_PI);
      for (ii=6; ii<8; ii++) {
        fprintf(outfile2, "\t%.15lf", p[ip*pperplan+ii]); 
      }
      fprintf(outfile2, "\n");
    }
    fprintf(outfile2, "%.15lf ; Mstar (M_sol)\n", p[npl*pperplan+0]);
    fprintf(outfile2, "%.15lf ; Rstar (R_sol)\n", p[npl*pperplan+1]);
    fprintf(outfile2, "%.15lf ; c1 (linear limb darkening) \n", p[npl*pperplan+2]);
    fprintf(outfile2, "%.15lf ; c2 (quadratic limb darkening) \n", p[npl*pperplan+3]);
    fprintf(outfile2, "%.15lf ; dilution (frac light not from stars in system)\n", p[npl*pperplan+4]);
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(RVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; rv jitter %i\n", p[npl*pperplan+5+ki], ki);
      }
    }
    if (TTVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(TTVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; ttv jitter %i\n", p[npl*pperplan+5+RVJITTERTOT+ki], ki);
      }
    }
    if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+ki]);
      }
    }
    if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; rvcelerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+ki]);
      }
    }
    fprintf(outfile2, " ; Comments: These coordinates are jacobian (see Lee & Peale 2003).\n");
    fprintf(outfile2, " ; neg2loglike = %.15lf\n", neg2loglike);

    if (SQRTE) {
      fprintf(outfile2, "planet         period (d)               T0 (d)              sqrt[e] cos(omega)        sqrt[e] sin(omega)        i (deg)                 Omega (deg)            mp (mjup)              rpors           ");
    } else {
      fprintf(outfile2, "planet         period (d)               T0 (d)              e cos(omega)             e sin(omega)             i (deg)                 Omega (deg)            mp (mjup)              rpors           ");
    }
    fprintf(outfile2, "\n");
    ch = 'a';
    pnum = 0.1;
    for (ip=0; ip<npl; ip++) {
      fprintf(outfile2, "%1.1lf", pnum);
      ch++; 
      pnum+=0.1;
      int ii;
      fprintf(outfile2, "\t%.15lf", orbparam[ip][0]); 
      fprintf(outfile2, "\t%.15lf", orbparam[ip][1]);
      if (SQRTE) {
        fprintf(outfile2, "\t%.15lf", sqrt(orbparam[ip][2])*cos(orbparam[ip][5]));
        fprintf(outfile2, "\t%.15lf", sqrt(orbparam[ip][2])*sin(orbparam[ip][5]));
      } else {
        fprintf(outfile2, "\t%.15lf", orbparam[ip][2]*cos(orbparam[ip][5]));
        fprintf(outfile2, "\t%.15lf", orbparam[ip][2]*sin(orbparam[ip][5]));
      }
      fprintf(outfile2, "\t%.15lf", orbparam[ip][3]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", orbparam[ip][4]*180.0/M_PI);
      for (ii=6; ii<8; ii++) {
        fprintf(outfile2, "\t%.15lf", p[ip*pperplan+ii]); 
      }
      fprintf(outfile2, "\n");
    }
    fprintf(outfile2, "%.15lf ; Mstar (M_sol)\n", p[npl*pperplan+0]);
    fprintf(outfile2, "%.15lf ; Rstar (R_sol)\n", p[npl*pperplan+1]);
    fprintf(outfile2, "%.15lf ; c1 (linear limb darkening) \n", p[npl*pperplan+2]);
    fprintf(outfile2, "%.15lf ; c2 (quadratic limb darkening) \n", p[npl*pperplan+3]);
    fprintf(outfile2, "%.15lf ; dilution (frac light not from stars in system)\n", p[npl*pperplan+4]);
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(RVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; rv jitter %i\n", p[npl*pperplan+5+ki], ki);
      }
    }
    if (TTVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(TTVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; ttv jitter %i\n", p[npl*pperplan+5+RVJITTERTOT+ki], ki);
      }
    }
    if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+ki]);
      }
    }
    if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; rvcelerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+ki]);
      }
    }
    fprintf(outfile2, " ; Comments: These coordinates are jacobian (see Lee & Peale 2003).\n");
    fprintf(outfile2, " ; neg2loglike = %.15lf\n", neg2loglike);


    fprintf(outfile2, "planet         a (AU)                   e                       i (deg)                omega (deg)          Omega (deg)               f (deg)               mp (mjup)              rpors           ");
    fprintf(outfile2, "\n");
    ch = 'a';
    pnum = 0.1;
    for (ip=0; ip<npl; ip++) {
      fprintf(outfile2, "%1.1lf", pnum);
      ch++; 
      pnum+=0.1;
      int ii;
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][0]); 
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][1]);
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][2]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][3]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][4]*180.0/M_PI);
      fprintf(outfile2, "\t%.15lf", aeiparam[ip][5]*180.0/M_PI);
      for (ii=6; ii<8; ii++) {
        fprintf(outfile2, "\t%.15lf", p[ip*pperplan+ii]); 
      }
      fprintf(outfile2, "\n");
    }
    fprintf(outfile2, "%.15lf ; Mstar (M_sol)\n", p[npl*pperplan+0]);
    fprintf(outfile2, "%.15lf ; Rstar (R_sol)\n", p[npl*pperplan+1]);
    fprintf(outfile2, "%.15lf ; c1 (linear limb darkening) \n", p[npl*pperplan+2]);
    fprintf(outfile2, "%.15lf ; c2 (quadratic limb darkening) \n", p[npl*pperplan+3]);
    fprintf(outfile2, "%.15lf ; dilution (frac light not from stars in system)\n", p[npl*pperplan+4]);
    if (RVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(RVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; rv jitter %i\n", p[npl*pperplan+5+ki], ki);
      }
    }
    if (TTVJITTERFLAG) {
      int ki;
      for (ki=0; ki<(TTVJITTERTOT); ki++) {
        fprintf(outfile2, "%.15lf ; ttv jitter %i\n", p[npl*pperplan+5+RVJITTERTOT+ki], ki);
      }
    }
    if (CELERITE) {
      int ki;
      for (ki=0; ki<NCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; celerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+ki]);
      }
    }
    if (RVCELERITE) {
      int ki;
      for (ki=0; ki<NRVCELERITE; ki++) {
        fprintf(outfile2, "%.15lf ; rvcelerite \n", p[npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+ki]);
      }
    }
    fprintf(outfile2, " ; Comments: These coordinates are jacobian (see Lee & Peale 2003).\n");
    fprintf(outfile2, " ; neg2loglike = %.15lf\n", neg2loglike);

    fclose(outfile2);

    for (i=0; i<npl; i++) free(aeiparam[i]);
    free(aeiparam);
    for (i=0; i<npl; i++) free(orbparam[i]);
    free(orbparam);

  }

  int i;
  free(int_in[0][0]);
  free(int_in[0]);
  free(int_in[1][0]);
  free(int_in[1]);
  int i1;
  for(i1=0; i1<npl+1; i1++) free(int_in[2][i1]);
  free(int_in[2]);
  free(int_in[3][0]);
  free(int_in[3]);
  free(int_in);

  //Only 1 array is malloced!
  free(flux_rvs[0][2]);
  if (RVS) {
    free(flux_rvs[1][2]);
  }
  free(flux_rvs[0]);
  free(flux_rvs[1]);
  free(flux_rvs);

#elif (demcmc_compile==1)
  printf("Setup Complete\n");
  printf("Beginning DEMCMC (this may take a very long time)\n");
  // xi squared array (one value per walker)
  double *neg2loglike0 = malloc(nwalkers*sofd);
  double **neg2loglike0N;

  for (i=0; i<nwalkers; i++) {
    double ***int_in = dsetup2(&p0local[pperwalker*i], npl);
    //printf("Converted to XYZ\n");
    double ***flux_rvs; 
    flux_rvs = dpintegrator_single(int_in, tfe, tve, nte, cadencelist);
    //printf("Computed Flux\n");
    double **ttvts = flux_rvs[2];
    double **flux = flux_rvs[0];
    double **radvs = flux_rvs[1];
    double *dev = devoerr(flux);
    double neg2logliketemp = 0;
    long il;
    long maxil = (long) dev[0];
    if (! CELERITE) { 
      for (il=0; il<maxil; il++) neg2logliketemp += dev[il+1]*dev[il+1];
    } else { // if celerite
      neg2logliketemp = celerite_fit(flux_rvs, p0local, i, 0, 0); 
    }
    double *newelist;
 
    if (RVS) {
      if (! RVCELERITE) { 
        double *newelist;
        if (RVJITTERFLAG) {
          long kj;
          long maxkj = (long) tve[0][0]; 
          newelist=malloc((maxkj+1)*sofd);
          newelist[0] = (double) maxkj;
          for (kj=0; kj<maxkj; kj++) {
            int jitterindex = (int) tve[3][1+kj]*NTELESCOPES + tve[4][1+kj];
            //double sigmajitter = p[npl][5+jitterindex]*MPSTOAUPD;
            double sigmajitter = p0local[pperwalker*i+npl*pperplan+5+jitterindex]*MPSTOAUPD;
            double quadsum = sigmajitter*sigmajitter + radvs[3][1+kj]*radvs[3][1+kj];
            // double check this... factor of 1/2
            neg2logliketemp += log(quadsum / (radvs[3][1+kj]*radvs[3][1+kj]) );
            newelist[1+kj] = sqrt( quadsum );
          }
          radvs[3] = newelist;
        }

        double *rvdev = devoerr(radvs);
        long maxil = (long) rvdev[0];
        for (il=0; il<maxil; il++) neg2logliketemp += rvdev[il+1]*rvdev[il+1];
        free(rvdev);
  
  
      } else { // if rvcelerite
        neg2logliketemp += celerite_fit(flux_rvs, p0local, i, 1, 0);
      } 
    }

    if (TTVCHISQ) {
      double *newelistttv;
      if (TTVJITTERFLAG) {
        long kj;
        long maxkj = (long) ttvts[0][0];
        newelistttv = malloc((maxkj+1)*sofd);
        newelistttv[0] = (double) maxkj;
        for (kj=0; kj<maxkj; kj++) {
          int jitterindex = (kj < NTTV[0][0]) ? 0 : 1 ; 
          double sigmajitter = p0local[pperwalker*i+npl*pperplan+5+RVJITTERTOT+jitterindex];
          double quadsum = sigmajitter*sigmajitter + ttvts[3][1+kj]*ttvts[3][1+kj];
          // check that the index on ttvts should be 2
          neg2logliketemp += log(quadsum / (ttvts[3][1+kj]*ttvts[3][1+kj]) );
          newelistttv[1+kj] = sqrt(quadsum);
        }
        ttvts[3] = newelistttv;
      }
      double *ttvdev = devoerr(ttvts);
      long maxil = (long) ttvdev[0];
      for (il=0; il<maxil; il++) neg2logliketemp += ttvdev[il+1]*ttvdev[il+1];
      free (ttvdev);
      if (TTVJITTERFLAG) {
        free(newelistttv);
      }
    }

    neg2logliketemp += compute_priors(p0local, i); 

    celeritefail:

    neg2loglike0[i] = neg2logliketemp;
 
    //free(evector);

    free(int_in[0][0]);
    free(int_in[0]);
    free(int_in[1][0]);
    free(int_in[1]);
    int i1;
    for(i1=0; i1<npl+1; i1++) free(int_in[2][i1]);
    free(int_in[2]);
    free(int_in[3][0]);
    free(int_in[3]);
    free(int_in);

    //Only 1 array is malloced!
    free(flux_rvs[0][2]);
    if (RVS) {
      free(flux_rvs[1][2]);
      if (RVJITTERFLAG) {
        free(newelist);
      }
    }

    free(flux_rvs[0]);
    free(flux_rvs[1]);
    free(flux_rvs);

    free(dev);


  }

  

  if (! RESTART ) {
    // set lowest neg2loglike
    neg2loglikemin = HUGE_VAL;
  }


  int k = 0;
  int *acceptance = malloc(nwalkers*sofi);
  int **acceptanceN;
  FILE *outfile;

  double *p0localcopy = malloc(pperwalker*sofd);

  int *acceptanceglobal = malloc(nwalkers*sofi);
  double *p0global = malloc(totalparams*sofd);
  double *neg2loglike0global = malloc(nwalkers*sofd);
  int **acceptanceglobalN;
  double **p0globalN;
  double **neg2loglike0globalN;


  // loop over generations
  unsigned long nwcore = (unsigned long) RANK;
  unsigned long nwalkersul = (unsigned long) nwalkers;
  gsl_rng_set(rnw, seed*(nwcore+2));

  while (jj<nsteps) {
 
    if (RANK==0 && jj % 10 == 0) {
      sout = fopen("demcmc.stdout", "a");
      fprintf(sout, "begin gen %li\n", jj);
      fclose(sout);
    }

    //time testing
    struct timespec start, finish;
    double elapsed;
    if (RANK==0 && jj % 10 ==0) { 
      clock_gettime(CLOCK_MONOTONIC, &start);
    }

    int nwi;
    int ncores;
    MPI_Comm_size(MPI_COMM_WORLD, &ncores);
    int npercore = nwalkers / ncores;
    
    //if ((nwalkers % ncores) > 0) {
    //  //printf("WARNING: Nwalkers is not an integer multiple of ncores. This will reduce performance.\n");
    //  //npercore += 1;
    //  printf("ERROR: Nwalkers modulo Ncores must be 0.\n");
    //  printf("       Nwalkers=%i, Ncores=%i\n", nwalkers, ncores);
    //  printf("       At least one of these must be changed.\n");
    //  exit(0);
    //}
    //nwcore = (unsigned long) RANK;
    //printf("ncores %i nwalkers %i\n", ncores, nwalkers);
    
    unsigned long nw;
    long nwinit=nwcore*npercore;
    long nwfin=(nwcore+1)*npercore;
    // This loops allows you to have more walkers than cores
    for (nw=nwinit; nw < nwfin; nw++) {

      memcpy(p0localcopy, &p0local[nw*pperwalker], pperwalker*sofd);

      acceptance[nw] = 1;
  
      unsigned long nw1;
      do nw1 = gsl_rng_uniform_int(rnw, nwalkersul); while (nw1 == nw); 
      unsigned long nw2;
      do nw2 = gsl_rng_uniform_int(rnw, nwalkersul); while (nw2 == nw || nw2 == nw1); 
  
      int notallowed=0;
  
      int ip;
      if (bimodf && jj % bimodf == 0) {
        for (ip=0; ip<pstar; ip++) {
          p0local[nw*pperwalker+npl*pperplan+ip] += (gamma+(1-gamma)*bimodlist[npl*pperplan+ip])*(p0local[nw1*pperwalker+npl*pperplan+ip]-p0local[nw2*pperwalker+npl*pperplan+ip])*(1-parfix[npl*pperplan+ip])*(1+(1+(gamma-1)*bimodlist[npl*pperplan+ip])*gsl_ran_gaussian(rnw, 0.1));
        }
        for (ip=0; ip<npl; ip++) {
          int ii;
          for (ii=0; ii<pperplan; ii++) {
            p0local[nw*pperwalker+ip*pperplan+ii] += (gamma+(1-gamma)*bimodlist[ip*pperplan+ii])*(p0local[nw1*pperwalker+ip*pperplan+ii]-p0local[nw2*pperwalker+ip*pperplan+ii])*(1-parfix[ip*pperplan+ii])*(1+(1+(gamma-1)*bimodlist[ip*pperplan+ii])*gsl_ran_gaussian(rnw, 0.1));
          }
        } 
      } else {
        for (ip=0; ip<pstar; ip++) {
          p0local[nw*pperwalker+npl*pperplan+ip] += gamma*(p0local[nw1*pperwalker+npl*pperplan+ip]-p0local[nw2*pperwalker+npl*pperplan+ip])*(1-parfix[npl*pperplan+ip])*(1+gsl_ran_gaussian(rnw, 0.1));
        }
        for (ip=0; ip<npl; ip++) {
          int ii;
          for (ii=0; ii<pperplan; ii++) {
            p0local[nw*pperwalker+ip*pperplan+ii] += gamma*(p0local[nw1*pperwalker+ip*pperplan+ii]-p0local[nw2*pperwalker+ip*pperplan+ii])*(1-parfix[ip*pperplan+ii])*(1+gsl_ran_gaussian(rnw, 0.1));
          }
        } 
      }
     
      // Check for chains that have strayed past any hard boundaries. 
      // No need to integrate any of these.  
      notallowed += check_boundaries(p0local, nw); 
 
      if (notallowed) { 
          
        acceptance[nw] = 0;
  
      } else {
      
        double ***nint_in = dsetup2(&p0local[nw*pperwalker], npl);
  
        double ***nflux_rvs; 
        nflux_rvs = dpintegrator_single(nint_in, tfe, tve, nte, cadencelist);
        double **nttvts = nflux_rvs[2];
        double **nflux = nflux_rvs[0];
        double **nradvs = nflux_rvs[1];
        double *ndev = devoerr(nflux);
        double nneg2logliketemp = 0;
        long il;
        long maxil = (long) ndev[0];
        if (! CELERITE) { 
          for (il=0; il<maxil; il++) nneg2logliketemp += ndev[il+1]*ndev[il+1];
        } else { // if celerite
          nneg2logliketemp = celerite_fit(nflux_rvs, p0local, nw, 0, 0);
        }
        double *newelist;
        if (RVS) {
          if (! RVCELERITE) { 
            double *newelist;
            if (RVJITTERFLAG) {
              long kj;
              long maxkj = (long) tve[0][0]; 
              newelist=malloc((maxkj+1)*sofd);
              newelist[0] = (double) maxkj;
              for (kj=0; kj<maxkj; kj++) {
                int jitterindex = (int) tve[3][1+kj]*NTELESCOPES + tve[4][1+kj];
                double sigmajitter = p0local[pperwalker*nw+npl*pperplan+5+jitterindex]*MPSTOAUPD;
                double quadsum = sigmajitter*sigmajitter + nradvs[3][1+kj]*nradvs[3][1+kj];
                // double check this... factor of 1/2
                nneg2logliketemp += log(quadsum / (nradvs[3][1+kj]*nradvs[3][1+kj]) );
                newelist[1+kj] = sqrt( quadsum );
              }
              nradvs[3] = newelist;
            }
    
            double *rvdev = devoerr(nradvs);
            long maxil = (long) rvdev[0];
            for (il=0; il<maxil; il++) nneg2logliketemp += rvdev[il+1]*rvdev[il+1];
            free(rvdev);
  
          } else { // if rvcelerite
            nneg2logliketemp += celerite_fit(nflux_rvs, p0local, nw, 1, 0);
          }
        }
        if (TTVCHISQ) {
          double *newelistttv;
          if (TTVJITTERFLAG) {
            long kj;
            long maxkj = (long) nttvts[0][0];
            newelistttv = malloc((maxkj+1)*sofd);
            newelistttv[0] = (double) maxkj;
            for (kj=0; kj<maxkj; kj++) {
              int jitterindex = (kj < NTTV[0][0]) ? 0 : 1 ; 
              //double sigmajitter = p[npl][5+RVJITTERTOT+jitterindex];
              double sigmajitter = p0local[pperwalker*nw+npl*pperplan+5+RVJITTERTOT+jitterindex];
              double quadsum = sigmajitter*sigmajitter + nttvts[3][1+kj]*nttvts[3][1+kj];
              // check that the index on ttvts should be 2
              nneg2logliketemp += log(quadsum / (nttvts[3][1+kj]*nttvts[3][1+kj]) );
              newelistttv[1+kj] = sqrt(quadsum);
            }
            nttvts[3] = newelistttv;
          }
        
          double *ttvdev = devoerr(nttvts);
          long maxil = (long) ttvdev[0];
          for (il=0; il<maxil; il++) nneg2logliketemp += ttvdev[il+1]*ttvdev[il+1];
          free (ttvdev);
          if (TTVJITTERFLAG) {
            free(newelistttv);
          }
  
        }
       

        nneg2logliketemp += compute_priors(p0local, nw);
 
        double neg2loglike = nneg2logliketemp;

        free(nint_in[0][0]);
        free(nint_in[0]);
        free(nint_in[1][0]);
        free(nint_in[1]);
        int i1;
        for(i1=0; i1<npl+1; i1++) free(nint_in[2][i1]);
        free(nint_in[2]);
        free(nint_in[3][0]);
        free(nint_in[3]);
        free(nint_in);
    
        free(nflux_rvs[0][2]);
        if (RVS) {
          free(nflux_rvs[1][2]);
          if (RVJITTERFLAG) {
            free(newelist);
          }
        }
    
        free(nflux_rvs[0]);
        free(nflux_rvs[1]);
        free(nflux_rvs);
    
        free(ndev);
    
        // prob that you should take new state
        double prob;
        prob = exp((neg2loglike0[nw]-neg2loglike)/2.);
 
        double bar = gsl_rng_uniform(rnw);
    
        // accept new state?
        if (prob <= bar || isnan(prob)) {
          acceptance[nw] = 0;
        } else {
          neg2loglike0[nw] = neg2loglike;
        } 
  
      }
      //free(evector);
  
      // switch back to old ones if not accepted
      if (acceptance[nw] == 0) {
        memcpy(&p0local[nw*pperwalker], p0localcopy, pperwalker*sofd);
      }

    }

    if (RANK==0 && jj % 10 ==0) {
      //time testing
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      sout = fopen("demcmc.stdout", "a");
      fprintf(sout, "Gen %li  run time = %3.12lf secs\n", jj, elapsed);
      fclose(sout);
    }

    int nread;
    if (nwfin-nwalkers > 0) {
      nread = nwalkers-nwinit;
      if (nread < 0) nread = 0;
    } else {
      nread = npercore;
    }

    MPI_Allgather(&p0local[nwinit*pperwalker], pperwalker*nread, MPI_DOUBLE, p0global, pperwalker*nread, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Gather(&acceptance[nwinit], 1*nread, MPI_INT, acceptanceglobal, 1*nread, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Gather(&neg2loglike0[nwinit], 1*nread, MPI_DOUBLE, neg2loglike0global, 1*nread, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    memcpy(p0local, p0global, totalparams*sofd);

    int nwdex;
    if (RANK == 0) {
      memcpy(acceptance, acceptanceglobal, nwalkers*sofi);
      memcpy(neg2loglike0, neg2loglike0global, nwalkers*sofd);

      long naccept = 0;
      for (nwdex=0; nwdex<nwalkers; nwdex++) naccept += acceptance[nwdex];

      double fracaccept = 1.0*naccept/nwalkers;
      if (fracaccept >= optimal) {
        gamma *= (1+relax);
      } else {
        gamma *= (1-relax);
      }

      char gammafstr[80];
      strcpy(gammafstr, "gamma_");
      strcat(gammafstr, OUTSTR);
      strcat(gammafstr, ".txt");
      if (jj % 10 == 0) {
        FILE *gammaf = fopen(gammafstr, "a");
        fprintf(gammaf, "%li\t%lf\t%lf\n", jj, fracaccept, gamma);
        fclose(gammaf);
      }

      // print out occasionally
      if (jj % OUTPUTINTERVAL == 0) { 
        char outfilestr[80];
        strcpy(outfilestr, "demcmc_");
        strcat(outfilestr, OUTSTR);
        strcat(outfilestr, ".out");
        outfile = fopen(outfilestr, "a");
        
        int nwdex;
        for (nwdex=0; nwdex<nwalkers; nwdex++) {
          char ch = 'a';
          double pnum = 0.1;
          for (i=0; i<npl; i++) {
            fprintf(outfile, "%1.1lf", pnum);
            ch++; 
            pnum+=0.1;
            int ii;
            for (ii=0; ii<pperplan; ii++) {
              fprintf(outfile, "\t%.15lf", p0local[nwdex*pperwalker+i*pperplan+ii]);  
            }
            fprintf(outfile, "\n");
          }
          fprintf(outfile, "%.15lf ; Mstar (R_sol)\n", p0local[nwdex*pperwalker+npl*pperplan+0]);
          fprintf(outfile, "%.15lf ; Rstar (R_sol)\n", p0local[nwdex*pperwalker+npl*pperplan+1]);
          fprintf(outfile, "%.15lf ; c1 (linear limb darkening) \n", p0local[nwdex*pperwalker+npl*pperplan+2]);
          fprintf(outfile, "%.15lf ; c2 (quadratic limb darkening) \n", p0local[nwdex*pperwalker+npl*pperplan+3]);
          fprintf(outfile, "%.15lf ; dilution (frac light not from stars in system)\n", p0local[nwdex*pperwalker+npl*pperplan+4]);
          if (RVJITTERFLAG) {
            for (i=0; i<RVJITTERTOT; i++) {
              fprintf(outfile, "%.15lf ; rvjitter \n", p0local[nwdex*pperwalker+npl*pperplan+5+i]);
            }
          }
          if (TTVJITTERFLAG) {
            for (i=0; i<TTVJITTERTOT; i++) {
              fprintf(outfile, "%.15lf ; rvjitter \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+i]);
            }
          }
          if (CELERITE) {
            for (i=0; i<NCELERITE; i++) {
              fprintf(outfile, "%.15lf ; celerite \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+i]);
            }
          }
          if (RVCELERITE) {
            for (i=0; i<NRVCELERITE; i++) {
              fprintf(outfile, "%.15lf ; rvcelerite \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+i]);
            }
          }
          fprintf(outfile, "; neg2loglike = %18.11lf, %i, %li\n", neg2loglike0[nwdex], nwdex, jj);  
        }
        fclose(outfile);
      }

      FILE *outfile2;
      for (nwdex=0; nwdex<nwalkers; nwdex++) {
        if (neg2loglike0[nwdex] < neg2loglikemin) {
          neg2loglikemin = neg2loglike0[nwdex];
          sout = fopen("demcmc.stdout", "a");
          fprintf(sout, "Chain %lu has Best neg2loglike so far: %lf\n", nwdex, neg2loglike0[nwdex]);
          fclose(sout);
          char outfile2str[80];
          strcpy(outfile2str, "mcmc_bestchisq_");
          strcat(outfile2str, OUTSTR);
          strcat(outfile2str, ".aei");
          outfile2 = fopen(outfile2str, "a");
          fprintf(outfile2, "planet         period (d)               T0 (d)                  e                   i (deg)                 Omega (deg)               omega(deg)               mp (mjup)              rpors           ");
          fprintf(outfile2, "\n");
          char ch = 'a';
          double pnum = 0.1;
          int ip;
          for (ip=0; ip<npl; ip++) {
            //fprintf(outfile2, "%c", ch);
            fprintf(outfile2, "%1.1lf", pnum);
            ch++; 
            pnum+=0.1;
            int ii;
            for (ii=0; ii<2; ii++) {
              fprintf(outfile2, "\t%.15lf", p0local[nwdex*pperwalker+ip*pperplan+ii]); 
            }
            if (SQRTE) {
              fprintf(outfile2, "\t%18.15lf", pow( sqrt( pow(p0local[nwdex*pperwalker+ip*pperplan+2],2) + pow(p0local[nwdex*pperwalker+ip*pperplan+3],2) ), 2) );
            } else {
              fprintf(outfile2, "\t%18.15lf", sqrt( pow(p0local[nwdex*pperwalker+ip*pperplan+2],2) + pow(p0local[nwdex*pperwalker+ip*pperplan+3],2) ) );
            }
            for (ii=4; ii<6; ii++) {
              fprintf(outfile2, "\t%.15lf", p0local[nwdex*pperwalker+ip*pperplan+ii]); 
            }
            fprintf(outfile2, "\t%18.15lf", atan2( p0local[nwdex*pperwalker+ip*pperplan+3] , p0local[nwdex*pperwalker+ip*pperplan+2] ) * 180./M_PI);
            for (ii=6; ii<8; ii++) {
              fprintf(outfile2, "\t%.15lf", p0local[nwdex*pperwalker+ip*pperplan+ii]); 
            }
            fprintf(outfile2, "\n");
          }
          fprintf(outfile2, "%.15lf ; Mstar (M_sol)\n", p0local[nwdex*pperwalker+npl*pperplan+0]);
          fprintf(outfile2, "%.15lf ; Rstar (R_sol)\n", p0local[nwdex*pperwalker+npl*pperplan+1]);
          fprintf(outfile2, "%.15lf ; c1 (linear limb darkening) \n", p0local[nwdex*pperwalker+npl*pperplan+2]);
          fprintf(outfile2, "%.15lf ; c2 (quadratic limb darkening) \n", p0local[nwdex*pperwalker+npl*pperplan+3]);
          fprintf(outfile2, "%.15lf ; dilution (frac light not from stars in system)\n", p0local[nwdex*pperwalker+npl*pperplan+4]);
          if (RVJITTERFLAG) {
            for (i=0; i<RVJITTERTOT; i++) {
              fprintf(outfile2, "%.15lf ; rvjitter \n", p0local[nwdex*pperwalker+npl*pperplan+5+i]);
            }
          }
          if (TTVJITTERFLAG) {
            for (i=0; i<TTVJITTERTOT; i++) {
              fprintf(outfile2, "%.15lf ; rvjitter \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+i]);
            }
          }
          if (CELERITE) {
            for (i=0; i<NCELERITE; i++) {
              fprintf(outfile2, "%.15lf ; celerite \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+i]);
            }
          }
          if (RVCELERITE) {
            for (i=0; i<NRVCELERITE; i++) {
              fprintf(outfile2, "%.15lf ; rvcelerite \n", p0local[nwdex*pperwalker+npl*pperplan+5+RVJITTERTOT+TTVJITTERTOT+NCELERITE*CELERITE+i]);
            }
          }
          fprintf(outfile2, " ; Comments: These coordinates are jacobian (see Lee & Peale 2003).   gen = %li   ch=%li\n", jj, nwdex);
          fprintf(outfile2, " ; neg2loglike = %.15lf\n", neg2loglike0[nwdex]);
          fclose(outfile2);
        }
      }


    }

    MPI_Bcast(&gamma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if ((RANK==0) && (jj % 10 == 0 )) {
      //time testing
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      sout = fopen("demcmc.stdout", "a");
      fprintf(sout, "Gen %li  run time + message passing = %3.12lf secs\n", jj, elapsed);
      fclose(sout);
    }

    jj+=1;

  }


  free(acceptance);
  free(acceptanceglobal);
  free(p0local);
  free(p0global);
  free(p0localcopy);
  free(neg2loglike0global);


  for(i=0; i<nwalkers; i++) {
    int i1;
    for (i1=0; i1<npl; i1++) {
      free(p0[i][i1]);
    }
    free(p0[i][npl]);
    free(p0[i]);
  }
  free(p0);

  free(neg2loglike0);

#elif (demcmc_compile == 3)

  int i;
  double ***int_in = dsetup2(p, npl);
  double ***flux_rvs; 
  flux_rvs = dpintegrator_single(int_in, tfe, tve, nte, cadencelist);
  free(int_in[0][0]);
  free(int_in[0]);
  free(int_in[1][0]);
  free(int_in[1]);
  int i1;
  for(i1=0; i1<npl+1; i1++) free(int_in[2][i1]);
  free(int_in[2]);
  free(int_in[3][0]);
  free(int_in[3]);
  free(int_in);

#endif
  if (TTVCHISQ) {

#if (demcmc_compile==0) 

    for (i=0; i<npl; i++) {
      char tempchar[1000];
      sprintf(tempchar, "ttvmcmc00_%02i.out", i+1);
      FILE* ttvouti = fopen(tempchar, "w");
      int ii;
      for (ii=0; ii<NTTV[i][0]; ii++) {
        fprintf(ttvouti, "%8li \t %.15lf \t %.15lf \t %.15lf \n", NTTV[i][ii+1], TTTV[i][ii+1], MTTV[i][ii+1], ETTV[i][ii+1]);
      }
      fclose(ttvouti);
    }

    for (i=0; i<npl; i++) free(MTTV[i]);
    free(MTTV);
#endif

    for (i=0; i<npl; i++) free(NTTV[i]);
    free(NTTV);
    for (i=0; i<npl; i++) free(TTTV[i]);
    free(TTTV);
    for (i=0; i<npl; i++) free(ETTV[i]);
    free(ETTV);
  }

  free(p);

  for (i=0; i<3; i++) free(tfe[i]);
  free(tfe);
  if (RVS) {
    for (i=0; i<5; i++) free(tve[i]);
  }
  if (TTVCHISQ) {
    int kij;
    for (kij=0; kij<3; kij++) {
      free(nte[kij]);
    }
    free(nte);    
  }
  free(tve);
  if (CADENCESWITCH==2) {
    free(cadencelist); 
  }

  gsl_rng_free(r);
  gsl_rng_free(rnw);

  return 0;
}





int main (int argc, char *argv[]) {

#if (demcmc_compile==1)
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
  MPI_Comm_rank(MPI_COMM_WORLD, &SIZE);
#endif

  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);

  getinput(argv[1]);
  
  void * nullptr;

  int sofi = SOFI;
  RVARR = calloc(NBODIES, sofi);
  RVFARR = malloc(NBODIES*sizeof(char*));
  int i;
  for (i=0; i<NBODIES; i++) {
    RVFARR[i] = malloc(100*sizeof(char));
  }
  int rvscount=0;
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-' && argv[i][1]=='r' && argv[i][2]=='v') {
      int body = argv[i][3]-'0';
      RVARR[body] = 1;
      strcpy(RVFARR[body], &argv[i][5]); 
      RVS = 1;
      rvscount++;
    }
  }
   
  argc-=rvscount;

  if (argc == 3) {
    RESTART = 0;
    if (RVS) {
      demcmc(argv[2], nullptr, nullptr, nullptr); 
    } else {
      demcmc(argv[2], nullptr, nullptr, nullptr); 
    }
  } else if (argc == 6) {
    RESTART = 1;
    if (RVS) {
      demcmc(argv[2], argv[3+rvscount], argv[4+rvscount], argv[5+rvscount]);
    } else {
      demcmc(argv[2], argv[3], argv[4], argv[5]);
    }
  } else {
    printf("usage: $ ./demcmc demcmc.in kep11.aei [[-rv0=rvs0.file] [-rv1=rvs1.file] ... ] [demcmc.res mcmcbsq.res gamma.res]\n");
    exit(0);
  }

  
  for(i=0; i<NBODIES; i++) {
    free(RVFARR[i]);
  }
  free(RVFARR);
  free(PARFIX); free(PSTEP); free(SSTEP); free(STEP); free(BIMODLIST); free(OUTSTR); free(RVARR);
  free(MAXDENSITY); free(EMAX); free(EPRIORV); free(MASSPRIORCENTERS); free(MASSPRIORSIGMAS);

  free(XYZLIST);

  clock_gettime(CLOCK_MONOTONIC, &finish);

  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Total run time = %3.12lf secs\n", elapsed);

#if (demcmc_compile==1)
  MPI_Finalize();
#endif

  return 0;

}



