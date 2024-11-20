#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"

#define sech 1/cosh
/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status;
  REAL z=0, dz0=0, r = GetValue(DATAFILE,"rstretch",&status), dz2=0;

  if(dz!=NULL) {
    z = 0;
    for(int N = 0; N<Nkmax; N++) {
      dz0 = 1.5*(tanh((N-Nkmax/3.0)/(Nkmax/30.0))+1.0) + 1.5;
      dz[N] = dz0;
      z += dz0;
    }
    for(int N = 0; N<Nkmax; N++) {
      dz[N] /= z;
      dz[N] *= depth;
    }
  } else {
      printf("Error in GetDZ: dz in NULL\n");
      exit(1);
  }
}

double gauss(double x, double a, double b, double c) {
    return  a * exp(-pow(((x - b) / c), 2));
}

double d_gauss(double x, double a, double b, double c) {
    return  a * exp(-pow(((x - b) / c), 2)) * 2 * (b - x) / pow(c, 2);
}
  
/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  	double z;
	z = 0;

	double a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, a5, b5, c5, a6, b6, c6, a7, b7, c7, a8, b8, c8;

	a1 = 7.900666e+02;
	b1 = 2.717665e-02;
	c1 = -6.374170e-01;
	a2 = 3.754591e+02;
	b2 = 5.443289e-02;
	c2 = 1.232468e+00;
	a3 = 1.180416e+02;
	b3 = 9.870790e-02;
	c3 = 1.833616e+00;
	a4 = 6.888052e+01;
	b4 = 1.406913e-01;
	c4 = 2.538620e+00;
	a5 = 9.968110e+01;
	b5 = 2.027221e-01;
	c5 = 4.999398e+00;
	a6 = 1.130123e+02;
	b6 = 1.917498e-01;
	c6 = 2.483625e+00;
	a7 = 9.585347e+00;
	b7 = 2.925321e-01;
	c7 = 3.719423e+00;
	a8 = 5.146615e+00;
	b8 = 4.044217e-01;
	c8 = 8.558507e-01;

	z = a1 * sin(b1 *(x) / 1000 + c1) + a2 * sin(b2 * (x) / 1000 + c2) + a3 * sin(b3 * (x) / 1000 + c3) + a4 * sin(b4 * (x) / 1000 + c4) + a5 * sin(b5 * (x) / 1000 + c5) + a6 * sin(b6 * (x) / 1000 + c6) + a7 * sin(b7 * (x) / 1000 + c7) + a8 * sin(b8 * (x) / 1000 + c8);

  return 957 - z;
}

 /*
  * Function: ReturnFreeSurface
  * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
  * -------------------------------------------------------------
  * Helper function to create an initial free-surface. Used
  * in phys.c in the InitializePhysicalVariables function.
  *
  */
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
  return 0;
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) {
    double a, b, c, a1, b1, d;
	double dnsty = 0;
	a = 2.661278e+01;
	b = 3.123267e-05;
	c = -2.554030e+00;
	d = -4.231268e-02;

	dnsty = (a * exp(b * (-z)) + c * exp(d * (-z)) - 25.8);

    return dnsty/RHO0;
}

/*
 * Function: ReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth) {
  return 1;
}

/*
 * Function: ReturnHorizontalVelocity
 * Usage: grid->u[n]=ReturnHorizontalVelocity(grid->xv[n],grid->yv[n],
 *                                            grid->n1[n],grid->n2[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial velocity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
#include "tide.h"

#define TIDE_TIME_SHIFT 0

#define Asin_m2  1.133329e-02
#define Asin_s2  -4.847622e-04
#define Asin_k1  -2.112834e-02
#define Asin_o1  -1.130596e-02
#define Asin_p1  -6.584461e-03
#define Asin_q1  -1.558616e-03

#define Acos_m2  9.246401e-03
#define Acos_s2  4.650254e-03
#define Acos_k1  5.675759e-04
#define Acos_o1  1.100382e-02
#define Acos_p1  9.158928e-04
#define Acos_q1  3.010823e-03

REAL U_left_bndy(REAL time)
{
  REAL U = 0;
  
  U += Asin_m2 * sin((2 * PI) / (12.42 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_m2 * cos((2 * PI) / (12.42 * 3600) * (time + TIDE_TIME_SHIFT)); //M2
  U += Asin_s2 * sin((2 * PI) / (12.00 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_s2 * cos((2 * PI) / (12.00 * 3600) * (time + TIDE_TIME_SHIFT)); //S2
  U += Asin_k1 * sin((2 * PI) / (23.93 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_k1 * cos((2 * PI) / (23.93 * 3600) * (time + TIDE_TIME_SHIFT)); //K1
  U += Asin_p1 * sin((2 * PI) / (24.07 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_p1 * cos((2 * PI) / (24.07 * 3600) * (time + TIDE_TIME_SHIFT)); //P1
  U += Asin_o1 * sin((2 * PI) / (25.82 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_o1 * cos((2 * PI) / (25.82 * 3600) * (time + TIDE_TIME_SHIFT)); //O1
  U += Asin_q1 * sin((2 * PI) / (26.87 * 3600) * (time + TIDE_TIME_SHIFT)) + Acos_q1 * cos((2 * PI) / (26.87 * 3600) * (time + TIDE_TIME_SHIFT)); //Q1
  
  return U;
}

REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z) {
  //normal to face!!
  REAL u, v;
  u = (ReturnDepth(0,0)/ReturnDepth(x,y)) * U_left_bndy(0);
  v = 0;
  
  return u*n1+v*n2; //from example
}

REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno) {
  return 0;
}
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno,int nsize) {
   return 0;
}
