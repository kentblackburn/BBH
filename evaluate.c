/* Evaluate.c */

/* James Kent Blackburn
   University of Florida
   Department of Physics
   Summer 1988            */

#include <stdio.h>
#include <math.h>
#include "assert.h"

#define fmax(a,b) (((a) > (b)) ? (a) : (b))
#define fmin(a,b) (((a) < (b)) ? (a) : (b))

/* double VaryIntegral(void)  calculates the entire variational integral.
   This integral is broken up into a number of parts.
   First the hot spots near each black hole are done. Then working
   from the inside out,
         From r = 0 to the inside edge of the more massive black hole.
         If the holes themselves overlap:
           The situation is ultrarelativistic.  Might as well just
           integrate from the origin out to infinity with no messing around.
           OR, for the time being, set the maximum radius of a given hole
           at its offset from the origin.
         else if the holes radii overlap :
           From inner edge of more massive to inner edge of less massive.
           From inner edge of less massive to outer edge of more massive
           From outer edge of more massive to outer edge of less massive.
         else they don't overlap:
           Across diameter of more massive.
           Between the holes.
           Across diameter of less massive.
         endif.
         Finally, from outer edge of outer hole to infinity.

double full_VaryIntegral(double xoff, double yoff, double zoff, double rmin,
                 double rmax)
does the full variational integral centered at the offset, from rmin to rmax.

double part_VaryIntegral(double x1, double r1, double x2, double r2,
                    double rmin, double rmax);
does the variational integral centered at the origin but avoids the
hot spots around each hole.
*/

extern FILE *fP;
extern int EQ_Masses_Flag;
extern long int TIMES;
extern double U, V, M1, M2, X1, X2, omega;
extern double Rmax;

double R1, R2; /* Radii of the hot spots around the holes. */

/*
double VolumeIntegral(double (*integrand)(double x, double y, double z),
                   double rmin, double rmax,
                   double xoff, double yoff, double zoff,
                   double x1, double r1, double x2, double r2, int holes);
*/

double full_VaryIntegral(xoff,yoff,zoff,rmin,rmax)

 double xoff,yoff,zoff,rmin,rmax;

{ 
 double VolumeIntegral(),VPintegrand();

 return VolumeIntegral(VPintegrand, rmin, rmax, xoff, yoff, zoff,
                        0.0, 0.0, 0.0, 0.0, 0); /* 0 is no holes.*/
}


double part_VaryIntegral(x1,r1,x2,r2,rmin,rmax)

 double x1,r1,x2,r2,rmin,rmax;

{
 double VolumeIntegral(),VPintegrand();

 return VolumeIntegral(VPintegrand, rmin, rmax, 0.0, 0.0, 0.0,
                        x1, r1, x2, r2, 1); /* 1 is hole(s). */
}


int dcomp(d1,d2)

 double *d1,*d2;

{
 if (*d1 >= *d2) return 1;
 else return -1;
}

double VaryIntegral()

{ 
 double r[4],absX1, absX2,outer2,VParound2;
 double VP1, VP2, VPinnermost, VPinner, VPmiddle, VPouter, VPoutermost;
 double full_VaryIntegral(),part_VaryIntegral();
 
 assert(fP != NULL);
 absX1 = fabs(X1);
 absX2 = fabs(X2);

 if (EQ_Masses_Flag)
   {
    R1 = fmin(2.0*M1,fabs(X2-X1)/4.0);
    R2 = R1;
   }
 else
   {
    R1 = fmin(2.0*M1, fabs(X2-X1)/4.0);
    R2 = fmin(100.0*M2, fabs(X2-X1)/4.0);
   }

 printf("\nM1 %g  R1 %g X1 %g\nM2 %g R2 %g X2 %g\n",
           M1, R1, X1, M2, R2, X2);
 fprintf(fP, "\nM1 %g  R1 %g X1 %g\nM2 %g R2 %g X2 %g\n",
                M1, R1, X1, M2, R2, X2);

 r[0] = fabs(absX1 - R1);
 r[1] = absX1 + R1;
 r[2] = fabs(absX2 - R2);
 r[3] = absX2 + R2;

 qsort(r, 4, sizeof(double), dcomp);

 printf("Sorted Radii %g %g %g %g\n", r[0], r[1], r[2], r[3]);
 fprintf(fP, "Sorted Radii %g %g %g %g\n", r[0], r[1], r[2], r[3]);
 TIMES = 0;
 VP1 = full_VaryIntegral(X1, 0.0, 0.0, 0.001*M1, R1);
 printf("TIMES %ld  Hole #1 %g\n", TIMES, VP1);
 fprintf(fP,"TIMES %ld  Hole #1 %g\n", TIMES, VP1);
 fflush(fP);

 if (M2 < 0.99*M1)
   {
    TIMES = 0;
    outer2 = fmin(10*M2,R2);
    VParound2 = full_VaryIntegral(X2, 0.0, 0.0, 0.001*M2, outer2);
    printf("TIMES %ld  Around Hole #2 %g\n", TIMES, VParound2);
    fprintf(fP, "TIMES %ld  Around Hole #2 %g\n", TIMES, VParound2);
    fflush(fP);
   }
 else
  {
   VParound2 = 0.0;
   outer2 = 0.001*M2;
  }

 TIMES = 0;
 if (EQ_Masses_Flag)
   VP2 = VP1;
 else
   VP2 = full_VaryIntegral(X2, 0.0, 0.0, outer2, R2);
 printf("TIMES %ld  Hole #2 %g\n", TIMES, VP2);
 fprintf(fP, "TIMES %ld  Hole #2 %g\n", TIMES, VP2);
 fflush(fP);

 if (R1 < absX1 && R2 < absX2)
   {
    TIMES = 0;
    VPinnermost = full_VaryIntegral(0.0, 0.0, 0.0, 0.0, r[0]);
    fprintf(fP, "TIMES %ld innermost %g\n", TIMES, VPinnermost);
    printf("TIMES %ld innermost %g\n", TIMES, VPinnermost);
    fflush(fP);
   }
 else
   VPinnermost = 0.0; /* else there is no innermost region. */

 TIMES = 0;
 VPinner = part_VaryIntegral(X1, R1, X2, R2, r[0], r[1]);
 printf("TIMES %ld inner %g\n",TIMES, VPinner);
 fprintf(fP,"TIMES %ld  inner %g\n", TIMES, VPinner);
 fflush(fP);

 TIMES = 0;
 VPmiddle = part_VaryIntegral(X1, R1, X2, R2, r[1], r[2]);
 printf("TIMES %ld Middle %g\n", TIMES, VPmiddle);
 fprintf(fP,"TIMES %ld Middle %g\n", TIMES, VPmiddle);
 fflush(fP);

 TIMES = 0;
 VPouter = part_VaryIntegral(X1, R1, X2, R2, r[2], r[3]);
 printf("TIMES %ld Outer %g\n", TIMES, VPouter);
 fprintf(fP,"TIMES %ld Outer %g\n", TIMES, VPouter);
 fflush(fP);
 TIMES = 0;
 VPoutermost = full_VaryIntegral(0.0, 0.0, 0.0, r[3], Rmax);

 printf("TIMES %ld outermost %g\n", TIMES, VPoutermost);
 fprintf(fP,"TIMES %ld outermost %g\n", TIMES, VPoutermost);
 fflush(fP);

 return VP1 + VParound2 + VP2 + VPinnermost + VPinner + VPmiddle
            + VPouter + VPoutermost;
}


