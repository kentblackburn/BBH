/* FastVol.c */
/* Has step size halving control. */

/* A quality controlled integrator over a spherical volume.  This allows for
   an offset and for the presence of two holes in the region of integration.
*/

/* James Kent Blackburn
   University of Florida
   Department of Physics
   Summer 1988            */

#include <stdio.h>
#include <math.h>
#include "assert.h"

#define M_PI      3.14159265358979323846
#define SMALL     1.0e-7
#define VERYSMALL 1.0e-10
#define JMAX     15 /* Max number halving. same as KMAX for fast but > KMAX otherwise*/
#define KMAX      7 /* Romberg of order 2K */
#define LMAX     64 /* 2 ^ (KMAX-1) */
#define KMIN      3

#define fmin(a,b) (((a) < (b)) ? (a) : (b))

extern FILE *fP;
extern long int TIMES;
extern double radError,
              thetaError,
              phiError;

int HOLES;

double local_x1,
       local_x2,
       local_r1,
       local_r2,
       local_xoff,
       local_yoff,
       local_zoff,
       local_radius,
       local_theta;

double  (*local_integrand)();


void polint(xa,ya,n,x,y,dy)

 int n;
 double *xa,*ya,x,*y,*dy;

{
 int i, m, ns = 1;
 double den, dif, dift, ho, hp, w;
 double c[JMAX+2], d[JMAX+2];

 dif = fabs(x-xa[1]);
 for (i = 1; i <= n; i++)
    {
     if ((dift = fabs(x - xa[i])) < dif)
       {
        ns = i;
        dif = dift;
       }
     c[i] = ya[i];
     d[i] = ya[i];
    }
 *y = ya[ns--];
 for (m = 1; m < n; m++)
    {
     for (i = 1; i <= n-m; i++)
        {
         ho = xa[i] - x;
         hp = xa[i+m] - x;
         w = c[i+1] - d[i];
         if ((den = ho-hp) == 0.0) assert(0);
         den = w/den;
         d[i] = hp * den;
         c[i] = ho * den;
        }
     *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
}


double trapzd(funct,a,b,n,current_s,current_it)

 int n,*current_it;
 double (*funct)(),a,b,*current_s;

{
 long int j;
 double tnm, del, sum, x, f;

 assert(a!=b);
 if (n == 1)
   {
    *current_it = 1; /* the number of points to be added on the next call. */
    return (*current_s = 0.5 * (b-a) * ((*funct)(a) + (*funct)(b)));
   }
 else
   {
    tnm = (double) *current_it;
    del = (b-a)/tnm;
    x = a + 0.5 * del;
    sum = 0.0;
    for (j = 0; j < *current_it; j++)
       {
        sum += (f = (*funct)(x));
        x = x + del;
       }
    *current_it *= 2;
    return (*current_s = 0.5 * ( *current_s + (b-a)*sum/tnm));
   }
}


double Qromb(integrand,a,b,error)

 double (*integrand)(),a,b,error;

{
 int curr_it,j, k;
 double ss, dss,s[JMAX+2], h[JMAX+2],curr_s;
 double trapzd();
 void polint();

 h[1] = 1.0;
 for (j = 1; j <= JMAX; j++)
    {
     s[j] = trapzd(integrand, a, b, j, &curr_s, &curr_it);
     if (j >= KMIN)
       {
        k = fmin(j, KMAX);
        polint(&h[j-k], &s[j-k], k, 0.0, &ss, &dss);
        if (fabs(dss) < error /* * fabs(ss) */)
          { /* success */
           return ss;
          }
/*
        print_inter = 1;
        polint(&h[j-k], &s[j-k], k, 0.0, &ss, &dss);
*/
       }
     s[j+1] = s[j];
     h[j+1] = 0.25*h[j];
    }
 assert(0); /* Too many steps. */
 return 0.0;
}

/* Qrombinf */
/* allows region of integration to extend out to infinity. */


#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))


double inf_trapzd(funk,f,a,b,n,it,s)

 int n,*it;
 double (*funk)(),*f,a,b,*s;

/* The array f[] contains all but the last needed refinement of the function. */

{
 int j,l,lstart,lplus;
 double tnm,del,sum,x;

 assert(a!=b);
 if (n == 1)
   {
    *it = 1; /* the number of points to be added on the next call. */
    return (*s = 0.5 * (b-a) * (f[0] + f[LMAX]));
   }
 else
   {
    tnm = (double) *it;
    if (n < KMAX)
      {
       lplus = LMAX / *it;
       lstart = lplus / 2;
       for (sum = 0.0, l = lstart; l < LMAX; l += lplus)
          {
           sum += f[l];
          }
      }
    else
      {
       del = (b-a)/tnm;
       x = a + 0.5 * del;
       assert(*it == LMAX/2);
       for (sum = 0.0, j = 0; j < *it; j++)
          {
           f[1+2*j] = FUNC(x);
/*         printf("trapzd 1+2*j %d  x %g  f[x] %g\n", 1+2*j, x, f[1+2*j]); */
           sum += f[1+2*j];
           x = x + del;
          }
      }
    *it *= 2;
    return (*s = 0.5 * ( *s + (b-a)*sum/tnm));
   }
}


double inf_QrombSlave(funk,f,a,b,error)

 double (*funk)(),*f,a,b,error;

{
 int curr_it,l,j;
 double curr_s,ss, dss,s[JMAX+2],h[JMAX+2],first_f[LMAX+1],first_s,second_s;
 double inf_trapzd(),inf_QrombSlave();
 void polint();

/*
 printf("enter QrombS from %g to %g\n", a, b);
 for (l = 0; l <= LMAX; l +=2)
    {
     printf("l %d  f %g\n", l, f[l]);
    }
*/
 h[1] = 1.0;
 for (j = 1; j <= JMAX; j++)
    {
     s[j] = inf_trapzd(funk, f, a, b, j, &curr_it, &curr_s);
     if (j >= KMAX)
       {
        polint(&h[j-KMAX], &s[j-KMAX], KMAX, 0.0, &ss, &dss);
/* printf("TestVol.inf_QrombSlave  polint returns ss %g error %g estError %g\n", 
                       ss, error, dss); */
        if (fabs(dss) < error /* * fabs(ss) */)
          { /* success */
           return ss;
          }
        else
          { /* halve the interval and recurse. */
           for (l = 0; l <= LMAX/2; l++)
              {
               first_f[2*l] = f[l];
              }
           first_s = inf_QrombSlave(funk, first_f, a, (b+a)/2.0, error);
           for (l = 0; l <= LMAX/2; l++)
              {
               first_f[2*l] = f[l+LMAX/2];
              }
           second_s = inf_QrombSlave(funk, first_f, (b+a)/2.0, b, error);
           return first_s + second_s;
          }
       }
     assert(j < KMAX);
     s[j+1] = s[j];
     h[j+1] = 0.25*h[j];
    }
 assert(0); /* Too many steps. */
 return 0;
}


double inf_QrombDriver(funk,aa,bb,error)

double (*funk)(),aa,bb,error;

{
 int l;
 double f[LMAX+1], del, a, b, x;
 double inf_QrombSlave();

 b = 1.0/aa;
 a = 1.0/bb;
 del = (b - a)/LMAX;
 for (l = 0; l <= LMAX; l += 2)
    {
     x = a + l * del;
     f[l] = FUNC(x);
    }
/* printf("inf_QrombDriver error %g\n", error); */
 return inf_QrombSlave(funk, f, a, b, error);
}


#undef FUNC

double New_trapzd(funct,f,a,b,n,it,s)

 int n,*it;
 double (*funct)(),*f,a,b,*s;

/* The array f[] contains all but the last needed refinement of the function. */

{
 int j, l, lstart, lplus;
 double tnm, del, sum, x;

 assert(a!=b);
 if (n == 1)
   {
    *it = 1; /* the number of points to be added on the next call. */
    return (*s = 0.5 * (b-a) * (f[0] + f[LMAX]));
   }
 else
   {
    tnm = (double) *it;
    if (n < KMAX)
      {
       lplus = LMAX / *it;
       lstart = lplus / 2;
       for (sum = 0.0, l = lstart; l < LMAX; l += lplus)
          {
           sum += f[l];
          }
      }
    else
      {
       del = (b-a)/tnm;
       x = a + 0.5 * del;
       assert(*it == LMAX/2);
       for (sum = 0.0, j = 0; j < *it; j++)
          {
           f[1+2*j] = (*funct)(x);
/*         printf("trapzd 1+2*j %d  x %g  f[x] %g\n", 1+2*j, x, f[1+2*j]); */
           sum += f[1+2*j];
           x = x + del;
          }
      }
    *it *= 2;
    return (*s = 0.5 * ( *s + (b-a)*sum/tnm));
   }
}


double QrombSlave(func,f,a,b,error)

 double (*func)(),*f,a,b,error;

{
 int j,l,curr_it;
 double curr_s,ss,dss,s[JMAX+2],h[JMAX+2],first_f[LMAX+1],first_s,second_s;
 double New_trapzd(),QrombSlave();
 void polint();
 
 h[1] = 1.0;
 for (j = 1; j <= JMAX; j++)
    {
     s[j] = New_trapzd(func, f, a, b, j, &curr_it, &curr_s);
     if (j >= KMAX)
       {
        polint(&h[j-KMAX], &s[j-KMAX], KMAX, 0.0, &ss, &dss);
        if (fabs(dss) < error /* * fabs(ss) */)
          { /* success */
           return ss;
          }
        else
          { /* halve the interval and recurse. */
           for (l = 0; l <= LMAX/2; l++)
              {
               first_f[2*l] = f[l];
              }
           first_s = QrombSlave(func, first_f, a, (b+a)/2.0, error);
           for (l = 0; l <= LMAX/2; l++)
              {
               first_f[2*l] = f[l+LMAX/2];
              }
           second_s = QrombSlave(func, first_f, (b+a)/2.0, b, error);
           return first_s + second_s;
          }
       }
     assert(j < KMAX);
     s[j+1] = s[j];
     h[j+1] = 0.25*h[j];
    }
 assert(0); /* Too many steps. */
 return 0;
}


double QrombDriver(func,a,b,error)

 double (*func)(),a,b,error;

{
 int l;
 double f[LMAX+1], del;
 double QrombSlave();

 for (l = 0; l <= LMAX; l += 2)
    {
     del = (b - a)/LMAX;
     f[l] = (*func)(a + l * del);
    }
 return QrombSlave(func, f, a, b, error);
}


double volume_element(rad,the,phi)

 double rad,the,phi;

{
 double x,y,z,Jacobian, ret;

 TIMES++;
/* if (!(TIMES % 5000)) printf("volume_element TIMES %ld\n", TIMES); */
/*
 x = rad*sin(the)*cos(phi);
 y = rad*sin(the)*sin(phi);
 z = rad*cos(the);
*/
/* We have the holes on the x-axis but the grid points are clustered near
   small values of theta.  So let (x,y,z) -> (y,z,x) to put the phi = 0 axis
   at the x-axis.
*/
 y = rad*sin(the)*cos(phi) + local_yoff;
 z = rad*sin(the)*sin(phi) + local_zoff;
 x = rad*cos(the) + local_xoff;
 Jacobian = rad*rad*sin(the);
 assert(local_integrand != NULL);
 ret = (*local_integrand)(x,y,z);
 ret = ret * Jacobian;
 return(ret);
}


double phi_integrand(phi)

 double phi;

{
 double result,volume_element();

 result = volume_element(local_radius, local_theta, phi);
 /*  printf("phi_integrand result %g\n", result); */
 return result;
}


double theta_integrand(theta)

 double theta;

{
 double result, phi_upper, phi_lower;
 double Qromb(),phi_integrand();

 local_theta = theta;
 if (sin(theta) < SMALL) return 0.0;
 phi_lower = 0.0; /* only one quarter is covered from symmetry. */
 phi_upper = M_PI/2.0;
 result = 4.0 * Qromb(phi_integrand, phi_lower, phi_upper, phiError);
 /*  printf("FastVol.theta_integrand theta %g  result %g\n",theta, result); */
 return result;
}


double radial_integrand(r)

 double r;

{
 double result,theta_upper,theta_lower,costheta;
 double QrombDriver(),theta_integrand();

 if (fabs(r) < VERYSMALL) return 0.0;
 local_radius = r;
 if (HOLES)
   {
    assert(local_x1 > local_x2); /* or else limits of costheta are wrong? */
    costheta = (local_x1 * local_x1 - local_r1 * local_r1 + r * r)
             / (2 * local_x1 * r);
    if (costheta < 1.0 && costheta > -1.0)
      {
       theta_lower = acos(costheta);
      }
    else
      {
       theta_lower = 0.0;
      }
    costheta = (local_x2 * local_x2 - local_r2 * local_r2 + r * r)
             / (2 * local_x2 * r);
    if (costheta < 1.0 && costheta > -1.0)
      {
       theta_upper = acos(costheta);
      }
    else
      {
       theta_upper = M_PI;
      }
   }
 else
   {
    theta_upper = M_PI;
    theta_lower = 0.0;
   }
 result = QrombDriver(theta_integrand, theta_lower, theta_upper, thetaError);
 /* printf("FastVol.radial_integrand r %g result %g\n", local_radius, result); */
 return result;
}


double VolumeIntegral(integrand,rmin,rmax,xoff,yoff,zoff,x1,r1,x2,r2,holes)

 int holes;
 double (*integrand)(),rmin,rmax,xoff,yoff,zoff,x1,r1,x2,r2;

{
 double result;
 double QrombDriver(),inf_QrombDriver(),radial_integrand();

 local_xoff = xoff;
 local_yoff = yoff;
 local_zoff = zoff;
 local_integrand = integrand;
 HOLES = holes;
 if (holes != 0)
   {
    local_x1 = x1;
    local_x2 = x2;
    local_r1 = r1;
    local_r2 = r2;
   }
 else
   {
    local_x1 = 0;
    local_x2 = 0;
    local_r1 = 0;
    local_r2 = 0;
   }
 printf("FastVol.VolumeIntegral from %f to %f\n", rmin, rmax); 
/*  printf("with xoff %f yoff %f zoff %f\n", local_xoff, local_yoff, local_zoff); */
/*  printf("x1 %f r1 %f x2 %f r2 %f\n", local_x1, local_r1, local_x2, local_r2);  */
 if(fabs(rmax-rmin) < VERYSMALL) return 0.0;
 if (HOLES == 0 && rmin > 0.1 && rmax > 20.0)
   { /* use inf functions. */
/*  printf("FastVol VolumeIntegral Use inf_functions.\n"); */
    result = inf_QrombDriver(radial_integrand, rmin, rmax, radError);
   }
 else
   {
    if (HOLES == 0 && local_xoff != 0.0 && rmin < SMALL) rmin = SMALL;
    result = QrombDriver(radial_integrand, rmin, rmax, radError);
   }
/* printf("FastVol.VolumeIntegral  result %g TIMES %ld\n", result, TIMES); */
 return result;
}



/*
double fun(double x, double y, double z)

{
 return x*x+y*y+z*z;
}

void main(void)

{
 double result;
 int holes = 0;
 double rmin, rmax, xoff=0.0, yoff=0.0, zoff=0.0, x1=0.0,
        r1=0.0, x2=0.0, r2=0.0;

 rmin = 1.0;
 rmax = 10.0;
 result = VolumeIntegral(fun, rmin, rmax, xoff, yoff, zoff, x1, r1, x2, r2, holes);
 printf("Volume Integral %g\n", result);
 printf(" should be %g\n", (4*M_PI)*(pow(rmax,5.0)-pow(rmin,5.0))/5.0);
}
*/

