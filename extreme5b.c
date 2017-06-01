/* Extreme5b.c */

/* Test algorithm for minimization routines.
   Calculates initial constants for trial geometry #4.
*/

/* James Kent Blackburn
   University of Florida
   Department of Physics
   Summer 1988            */

#include <stdio.h>
#include <math.h>
#include "assert.h"

#define YES          1
#define NO           0
#define FIND_MIN     YES
#define STEST        0.75
#define TOL          1.0e-5 /* change 5 to 6 */
#define INTGRLEPS    1.0e-7 /* change to -8 */
#define PI           3.14159265359
#define EPSILON      1.0e-7
#define FILENAME     "neq.data"
#define FTABULAR     "table100.data"
#define M1IRRBLE     1.00
#define M2IRRBLE     1.00
#define JDENSITY     7.75
#define INNERS     105.00
#define OUTERS     115.00 /* start at 4.90 for new unstable searches */
#define SIZES        2.00

FILE *fP;
FILE *fT;
long int TIMES;

int EQ_Masses_Flag = 0;

double Rmax   =  100000.00,
       M1irr  =   M1IRRBLE,
       M2irr  =   M2IRRBLE,
       Sbegin =     OUTERS,
       Send   =     INNERS,
       DeltaS =     -SIZES,
       Sstart =      6*(M1IRRBLE+M2IRRBLE),
       Jtotal =      JDENSITY*M1IRRBLE*M2IRRBLE;
                  /*(-) set in 1st call to initialization(s) */

double radError =   INTGRLEPS,
       thetaError = INTGRLEPS,
       phiError =   INTGRLEPS;
         
double U, V, M1, M2, X1, X2, omega; /* Set in prmtrzd */

int check_initialization(s,jtotal)

 double s,jtotal;

{
 if (  (fabs(1 - (X1 - X2)/s) > EPSILON)
    ||  (fabs(1 - V/X2/omega) > EPSILON)
    ||  (fabs(1 - U/X1/omega) > EPSILON)
    ||  (fabs(1 - (M1*U*X1/sqrt(1-U*U) + M2*V*X2/sqrt(1-V*V))/jtotal)>EPSILON)
    ||  (fabs(M1*U/sqrt(1-U*U) + M2*V/sqrt(1-V*V)) > EPSILON)
    ||  (fabs(M1irr - M1*(1+M2/(2*s))) > EPSILON)
    ||  (fabs(M2irr - M2*(1+M1/(2*s))) > EPSILON))
   {
    printf("Check Setup.\n\n  U/X1 = %f   V/X2 = %f   s = %f   X1-X2 = %f\n\n",
            omega, V/X2, s, X1 - X2);
    printf("jtotal = %f  ?= J1  + J2 =  %f  +  %f  = %f\n\n", jtotal,
            M1*U*X1/sqrt(1-U*U), M2*V*X2/sqrt(1-V*V) ,
            M1*U*X1/sqrt(1-U*U) + M2*V*X2/sqrt(1-V*V));
    printf("Ptotal = P1 + P2 = %f  +  %f = %f\n\n",
            M1*U/sqrt(1-U*U), M2*V/sqrt(1-V*V) ,
            M1*U/sqrt(1-U*U) + M2*V/sqrt(1-V*V));
    printf("M1irr = %f = M1 * (1 + M2/(2*s)) = %f\n",
            M1irr, M1 * (1 + M2/(2*s)));
    printf("M2irr = %f = M2 * (1 + M1/(2*s)) = %f\n",
            M2irr, M2 * (1 + M1/(2*s)));
    printf("U = %f   V = %f\n\n", U, V);
    printf("\nHit a Key ( <x Return> to cancel) . . .\n");
    if ('x' == getchar()) return 0;
    else return 1;
   }
 return 1;
}

int initialization(s)

 double s;

{
 int itter;
 double b,c,tmp1,tmp2;
 double mu, Eschw, Jschw, Sschw;
 double Mtotal;

 if (Jtotal <= 0.0)
   { /* Fix the total angular momentum once and for all. */

     /* For Test particle around Schwarzschild,
              isotropic radial coord. = Sstart. */

    Mtotal = M1irr + M2irr;
    mu = M1irr * M2irr /Mtotal;

    Sschw = Sstart * (1+Mtotal/(2*Sstart))*(1+Mtotal/(2*Sstart));

    Eschw = mu * (Sschw - 2*Mtotal) / sqrt(Sschw*(Sschw - 3*Mtotal));
    Jschw = mu * Sschw * sqrt(Mtotal/(Sschw - 3*Mtotal));

    Jtotal = Jschw;
    tmp1 = (1-2*Mtotal/Sschw) * (1 + pow(Jtotal/(mu*Sschw),2.0));
    tmp2 = (1-2*Mtotal/Sschw) * (1 + pow(Jtotal/(Mtotal*Sschw),2.0));

    printf("Extreme.initialization  Sstart %g   Sschw %g  Eschw %g\n",
              Sstart, Sschw, Eschw);
    printf( "Compare Eschw %g with [ ]^1/2 mu %g\n", Eschw,
          sqrt(tmp1)*mu );

    printf("  Jtotal/(mu*Mtotal) %g < 2 sqrt(3) %g\n",
                                  Jtotal/(mu*Mtotal), 2*sqrt(3.0));
    assert(fP != NULL);
    fprintf(fP,"Extreme.initialization  Sstart %g   Sschw %g  Eschw %g\n",
              Sstart, Sschw, Eschw);
    fprintf(fP,"Compare Eschw %g with [ ]^1/2 mu + M %g\n", Eschw,
          sqrt(tmp2)*mu);

    fprintf(fP,"  Jtotal/(mu*Mtotal) %g < 2 sqrt(3) %g\n",
                                  Jtotal/(mu*Mtotal), 2*sqrt(3.0));
   }

 b = 2*s - M1irr + M2irr;
 c = -2*s*M1irr;
 M1 = (-b+sqrt(b*b-4*c))/2.0;
 b = 2*s - M2irr + M1irr;
 c = -2*s*M2irr;
 M2 = (-b+sqrt(b*b-4*c))/2.0;

 U =   Jtotal / s / sqrt(M1*M1 + Jtotal*Jtotal/(s*s));
 V = - Jtotal / s / sqrt(M2*M2 + Jtotal*Jtotal/(s*s));
 X1 =   M2 * sqrt(1-U*U) * s / (M1 * sqrt(1-V*V) + M2 * sqrt(1-U*U));
 X2 = - M1 * sqrt(1-V*V) * s / (M1 * sqrt(1-V*V) + M2 * sqrt(1-U*U));
 omega = U / X1;
 if (!check_initialization(s,Jtotal)) return 0;
 else return 1;
}

main()

{
 double a,b,c,fa,fb,fc,fm,prmtrzd(),xmin,fmin,brent();
 double s,brent(),prmtrzd();

 fP = fopen(FILENAME, "a+");
 assert(fP != NULL);
 if ( fabs(M1irr - M2irr) < EPSILON ) EQ_Masses_Flag = 1;
 fmin = M1irr + M2irr;
 assert (fP != NULL);
 fprintf(fP,"\nM1irr %g M2irr %g SNewt %g Rmax %g\n DeltaS %g Error %g\n",
               M1irr, M2irr, Sstart, Rmax, DeltaS, radError);


                    /* For scanning */
 /* NOTE: Must scane hi to lo values of separation! */

 for (s = Sbegin ; s >= Send; s += DeltaS)
    {
     fm = prmtrzd(s);
     if ((fm < fmin)&&(s > STEST))
       {
        a = s;
        b = a + DeltaS/5.0;
       }
     fmin = fm;
    }

        /* For finding the extremum */
 if (FIND_MIN == YES)
   {
    mnbrak(&a,&b,&c,&fa,&fb,&fc,prmtrzd);
    printf("f(x) = %g  at  x = %g \n",fa,a);
    printf("f(x) = %g  at  x = %g \n",fb,b);
    printf("f(x) = %g  at  x = %g \n",fc,c);

    fprintf(fP,"f(x) = %g  at  x = %g \n",fa,a);
    fprintf(fP,"f(x) = %g  at  x = %g \n",fb,b);
    fprintf(fP,"f(x) = %g  at  x = %g \n",fc,c);
    fflush(fP);

    fmin = brent(a,b,fb,c,prmtrzd,TOL,&xmin);

    if ((fT = fopen(FTABULAR,"a+")) != NULL)
      {
       fprintf(fT,"%1d\n",1);
       fclose(fT);
      }

    printf("Extremized Mass = %15.9g  at  the Separation = %15.9f\n",fmin,xmin);
    fprintf(fP,"\n Extremized Mass = %15.9g  at  the Separation = %15.9f\n",fmin,xmin);
   }
 fclose(fP);
}

double prmtrzd(s)

 double s;

{
 int initialization();
 double Minf,Mmax,mVP;
 double VaryIntegral();

 if (!initialization(s)) return 0.0;

 Minf = M1 / sqrt(1-U*U) + M2 / sqrt(1-V*V);
 mVP = VaryIntegral()/(16*PI) - omega * Jtotal; /* Second term is
                                        from Version A of the V.P. */
 Mmax = Minf - mVP;
 printf("Jtotal %f Minf %13.9f - mVP %13.9f = %13.9f S = %f\n",
         Jtotal,Minf, mVP, Mmax, X1-X2);
 fprintf(fP, "Jtotal %f Minf %13.9f - mVP %13.9f = %13.9f S = %f\n",
              Jtotal,Minf, mVP, Mmax, X1-X2);
 printf("U = %g V = %g and OMEGA = %g \n",U,V,omega);
 fprintf(fP,"U = %13.9g V = %13.9g and OMEGA = %13.9g \n",U,V,omega);
 fflush(fP);
 if ((fT = fopen(FTABULAR,"a+")) != NULL)
   {
    fprintf(fT,"\n%7.4f %7.4f %11.8f %18.14f %18.14f %18.14f",
                  M1irr,M2irr,Jtotal,X1-X2,  Mmax,   omega);
    fprintf(fT,"\n   %17.14f %17.14f %17.14f %17.14f %17.14f %17.14f %1d",
                      M1,     M2,     Minf,   mVP,    U,      V,      0);
    fclose(fT);
   }
 return(Mmax);
}


