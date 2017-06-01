/* scalar5b.c */

/* James Kent Blackburn
   University of Florida
   Department of Physics
   Summer 1988            */

#include <math.h>
#include <stdio.h>

double VPintegrand(x,y,z)
    double x,y,z;
{
  double rho1, drho1_dx, drho1_dy, drho1_dz, d2rho1_dxdx, d2rho1_dydy,
      d2rho1_dzdz, d2rho1_dxdy, d2rho1_dxdz, d2rho1_dydz, d2rho1_dydx,
      d2rho1_dzdx, d2rho1_dzdy ;
  double rho2, drho2_dx, drho2_dy, drho2_dz, d2rho2_dxdx, d2rho2_dydy,
      d2rho2_dzdz, d2rho2_dxdy, d2rho2_dxdz, d2rho2_dydz, d2rho2_dydx,
      d2rho2_dzdx, d2rho2_dzdy ;

  double common, dcommon_drho1, dcommon_drho2, dcommon_dx, dcommon_dy,
      dcommon_dz ;
  double d2common_drho1drho1, d2common_drho1drho2, d2common_drho2drho2,
      d2common_dxdx, d2common_dydy, d2common_dzdz, d2common_dxdy,
      d2common_dxdz, d2common_dydz, d2common_dydx, d2common_dzdx, d2common_dzdy
      ;

  double PSI, dPSI_drho1, dPSI_drho2, dPSI_dx, dPSI_dy, dPSI_dz ;
  double d2PSI_drho1drho1, d2PSI_drho1drho2, d2PSI_drho2drho2, d2PSI_dxdx,
      d2PSI_dydy, d2PSI_dzdz, d2PSI_dxdy, d2PSI_dxdz, d2PSI_dydz, d2PSI_dydx,
      d2PSI_dzdx, d2PSI_dzdy ;

  double psi1, dpsi1_drho1, dpsi1_drho2, dpsi1_dx, dpsi1_dy, dpsi1_dz ;
  double d2psi1_drho1drho1, d2psi1_drho1drho2, d2psi1_drho2drho2, d2psi1_dxdx,
      d2psi1_dydy, d2psi1_dzdz, d2psi1_dxdy, d2psi1_dxdz, d2psi1_dydz,
      d2psi1_dydx, d2psi1_dzdx, d2psi1_dzdy ;

  double psi2, dpsi2_drho1, dpsi2_drho2, dpsi2_dx, dpsi2_dy, dpsi2_dz ;
  double d2psi2_drho1drho1, d2psi2_drho1drho2, d2psi2_drho2drho2, d2psi2_dxdx,
      d2psi2_dydy, d2psi2_dzdz, d2psi2_dxdy, d2psi2_dxdz, d2psi2_dydz,
      d2psi2_dydx, d2psi2_dzdx, d2psi2_dzdy ;

  double QU_PSI, dQU_PSI_drho1, dQU_PSI_drho2, dQU_PSI_dx, dQU_PSI_dy,
      dQU_PSI_dz ;
  double d2QU_PSI_drho1drho1, d2QU_PSI_drho1drho2, d2QU_PSI_drho2drho2,
      d2QU_PSI_dxdx, d2QU_PSI_dydy, d2QU_PSI_dzdz, d2QU_PSI_dxdy,
      d2QU_PSI_dxdz, d2QU_PSI_dydz, d2QU_PSI_dydx, d2QU_PSI_dzdx, d2QU_PSI_dzdy
      ;

  double qu_psi1, dqu_psi1_drho1, dqu_psi1_drho2, dqu_psi1_dx, dqu_psi1_dy,
      dqu_psi1_dz ;
  double d2qu_psi1_drho1drho1, d2qu_psi1_drho1drho2, d2qu_psi1_drho2drho2,
      d2qu_psi1_dxdx, d2qu_psi1_dydy, d2qu_psi1_dzdz, d2qu_psi1_dxdy,
      d2qu_psi1_dxdz, d2qu_psi1_dydz, d2qu_psi1_dydx, d2qu_psi1_dzdx,
      d2qu_psi1_dzdy ;

  double qu_psi2, dqu_psi2_drho1, dqu_psi2_drho2, dqu_psi2_dx, dqu_psi2_dy,
      dqu_psi2_dz ;
  double d2qu_psi2_drho1drho1, d2qu_psi2_drho1drho2, d2qu_psi2_drho2drho2,
      d2qu_psi2_dxdx, d2qu_psi2_dydy, d2qu_psi2_dzdz, d2qu_psi2_dxdy,
      d2qu_psi2_dxdz, d2qu_psi2_dydz, d2qu_psi2_dydx, d2qu_psi2_dzdx,
      d2qu_psi2_dzdy ;

  double a1, da1_drho1, da1_drho2, da1_dx, da1_dy, da1_dz ;
  double d2a1_drho1drho1, d2a1_drho1drho2, d2a1_drho2drho2, d2a1_dxdx,
      d2a1_dydy, d2a1_dzdz, d2a1_dxdy, d2a1_dxdz, d2a1_dydz, d2a1_dydx,
      d2a1_dzdx, d2a1_dzdy ;

  double a2, da2_drho1, da2_drho2, da2_dx, da2_dy, da2_dz ;
  double d2a2_drho1drho1, d2a2_drho1drho2, d2a2_drho2drho2, d2a2_dxdx,
      d2a2_dydy, d2a2_dzdz, d2a2_dxdy, d2a2_dxdz, d2a2_dydz, d2a2_dydx,
      d2a2_dzdx, d2a2_dzdy ;

  double alpha1, dalpha1_drho1, dalpha1_drho2, dalpha1_dx, dalpha1_dy,
      dalpha1_dz ;
  double d2alpha1_drho1drho1, d2alpha1_drho1drho2, d2alpha1_drho2drho2,
      d2alpha1_dxdx, d2alpha1_dydy, d2alpha1_dzdz, d2alpha1_dxdy,
      d2alpha1_dxdz, d2alpha1_dydz, d2alpha1_dydx, d2alpha1_dzdx, d2alpha1_dzdy
      ;

  double alpha2, dalpha2_drho1, dalpha2_drho2, dalpha2_dx, dalpha2_dy,
      dalpha2_dz ;
  double d2alpha2_drho1drho1, d2alpha2_drho1drho2, d2alpha2_drho2drho2,
      d2alpha2_dxdx, d2alpha2_dydy, d2alpha2_dzdz, d2alpha2_dxdy,
      d2alpha2_dxdz, d2alpha2_dydz, d2alpha2_dydx, d2alpha2_dzdx, d2alpha2_dzdy
      ;

  double sq_alpha1, dsq_alpha1_drho1, dsq_alpha1_drho2, dsq_alpha1_dx,
      dsq_alpha1_dy, dsq_alpha1_dz ;
  double d2sq_alpha1_drho1drho1, d2sq_alpha1_drho1drho2,
      d2sq_alpha1_drho2drho2, d2sq_alpha1_dxdx, d2sq_alpha1_dydy,
      d2sq_alpha1_dzdz, d2sq_alpha1_dxdy, d2sq_alpha1_dxdz, d2sq_alpha1_dydz,
      d2sq_alpha1_dydx, d2sq_alpha1_dzdx, d2sq_alpha1_dzdy ;

  double sq_alpha2, dsq_alpha2_drho1, dsq_alpha2_drho2, dsq_alpha2_dx,
      dsq_alpha2_dy, dsq_alpha2_dz ;
  double d2sq_alpha2_drho1drho1, d2sq_alpha2_drho1drho2,
      d2sq_alpha2_drho2drho2, d2sq_alpha2_dxdx, d2sq_alpha2_dydy,
      d2sq_alpha2_dzdz, d2sq_alpha2_dxdy, d2sq_alpha2_dxdz, d2sq_alpha2_dydz,
      d2sq_alpha2_dydx, d2sq_alpha2_dzdx, d2sq_alpha2_dzdy ;

  double gxx, dgxx_drho1, dgxx_drho2, dgxx_dx, dgxx_dy, dgxx_dz ;
  double d2gxx_drho1drho1, d2gxx_drho1drho2, d2gxx_drho2drho2, d2gxx_dxdx,
      d2gxx_dydy, d2gxx_dzdz, d2gxx_dxdy, d2gxx_dxdz, d2gxx_dydz, d2gxx_dydx,
      d2gxx_dzdx, d2gxx_dzdy ;
  double gyy, dgyy_drho1, dgyy_drho2, dgyy_dx, dgyy_dy, dgyy_dz ;
  double d2gyy_drho1drho1, d2gyy_drho1drho2, d2gyy_drho2drho2, d2gyy_dxdx,
      d2gyy_dydy, d2gyy_dzdz, d2gyy_dxdy, d2gyy_dxdz, d2gyy_dydz, d2gyy_dydx,
      d2gyy_dzdx, d2gyy_dzdy ;
  double gzz, dgzz_drho1, dgzz_drho2, dgzz_dx, dgzz_dy, dgzz_dz ;
  double d2gzz_drho1drho1, d2gzz_drho1drho2, d2gzz_drho2drho2, d2gzz_dxdx,
      d2gzz_dydy, d2gzz_dzdz, d2gzz_dxdy, d2gzz_dxdz, d2gzz_dydz, d2gzz_dydx,
      d2gzz_dzdx, d2gzz_dzdy ;

  double beta_X,beta_Y,beta_Z,beta_z,gXX,gYY,gZZ;

  double lapse, K, K_AB_K_ab;
  double K1xx, K1yy, K1zz, K1xy, K1yx, K1yz, K1zy, K1xz, K1zx;
  double K2xx, K2yy, K2zz, K2xy, K2yx, K2yz, K2zy, K2xz, K2zx;
  double K_xx, K_yy, K_zz, K_xy, K_yx, K_yz, K_zy, K_xz, K_zx;
  double K_XX,K_YY,K_ZZ,K_XY,K_YX,K_XZ,K_ZX,K_YZ,K_ZY;

  double beta_y, dbeta_y_drho1, dbeta_y_drho2, dbeta_y_dx, dbeta_y_dy,
      dbeta_y_dz ;
  double beta_x, dbeta_x_drho1, dbeta_x_drho2, dbeta_x_dx, dbeta_x_dy,
      dbeta_x_dz ;
  double D_x_beta_y, D_y_beta_y, D_z_beta_y,
      D_x_beta_x, D_y_beta_x, D_z_beta_x, D_x_beta_z, D_y_beta_z, D_z_beta_z;


  double G_xxx,G_xxy,G_xyx,G_xxz,G_xzx,G_xyy,G_xyz,G_xzy,G_xzz,
      G_yxx,G_yxy,G_yyx,G_yxz,G_yzx,G_yyy,G_yyz,G_yzy,G_yzz,
      G_zxx,G_zxy,G_zyx,G_zxz,G_zzx,G_zyy,G_zyz,G_zzy,G_zzz,

  G_Xxx,G_Xxy,G_Xyx,G_Xxz,G_Xzx,G_Xyy,G_Xyz,G_Xzy,G_Xzz,
      G_Yxx,G_Yxy,G_Yyx,G_Yxz,G_Yzx,G_Yyy,G_Yyz,G_Yzy,G_Yzz,
      G_Zxx,G_Zxy,G_Zyx,G_Zxz,G_Zzx,G_Zyy,G_Zyz,G_Zzy,G_Zzz;

  double R, S2, S;

  extern double U,V,M1,M2,X1,X2,omega;


  {
    if (0) exit(0);
    S2 = (X1-X2)*(X1-X2);
    S = sqrt(S2);

    rho1 = sqrt( (x-X1)*(x-X1) + y*y/(1-U*U) + z*z );
    drho1_dx = (x-X1)/rho1;
    drho1_dy = y/((1-U*U)*rho1);
    drho1_dz = z/rho1;

    d2rho1_dxdx = 1/rho1 - (x-X1)/(rho1*rho1) * drho1_dx;
    d2rho1_dydy = 1/(1-U*U)/rho1 - y/(1-U*U)/(rho1*rho1) * drho1_dy;
    d2rho1_dzdz = 1/rho1 - z/(rho1*rho1) * drho1_dz;
    d2rho1_dydx = d2rho1_dxdy = -drho1_dx * drho1_dy/rho1;
    d2rho1_dzdy = d2rho1_dydz = -drho1_dy * drho1_dz/rho1;
    d2rho1_dzdx = d2rho1_dxdz = -drho1_dz * drho1_dx/rho1;

    rho2 = sqrt( (x-X2)*(x-X2) + y*y/(1-V*V) + z*z );
    drho2_dx = (x-X2)/rho2;
    drho2_dy = y/(1-V*V)/rho2;
    drho2_dz = z/rho2;

    d2rho2_dxdx = 1/rho2 - (x-X2)/(rho2*rho2) * drho2_dx;
    d2rho2_dydy = 1/(1-V*V)/rho2 - y/(1-V*V)/(rho2*rho2) * drho2_dy;
    d2rho2_dzdz = 1/rho2 - z/(rho2*rho2) * drho2_dz;
    d2rho2_dydx = d2rho2_dxdy = -drho2_dx * drho2_dy/rho2;
    d2rho2_dzdy = d2rho2_dydz = -drho2_dy * drho2_dz/rho2;
    d2rho2_dzdx = d2rho2_dxdz = -drho2_dz * drho2_dx/rho2;


    PSI = 1 + M1/(2*rho1) + M2/(2*rho2);
    dPSI_drho1 = -M1/(2*rho1*rho1);
    dPSI_drho2 = -M2/(2*rho2*rho2);
    d2PSI_drho1drho1 = M1/((rho1)*(rho1)*(rho1)) ;
    d2PSI_drho2drho2 = M2/((rho2)*(rho2)*(rho2)) ;
    d2PSI_drho1drho2 = 0;


    QU_PSI = ((PSI)*(PSI)*(PSI)*(PSI)) ;
    dQU_PSI_drho1 = 4*(QU_PSI/PSI)*dPSI_drho1;
    dQU_PSI_drho2 = 4*(QU_PSI/PSI)*dPSI_drho2;

    d2QU_PSI_drho1drho1 = 12*((PSI*dPSI_drho1)*(PSI*dPSI_drho1))
        + 4*(QU_PSI/PSI)*d2PSI_drho1drho1;
    d2QU_PSI_drho2drho2 = 12*((PSI*dPSI_drho2)*(PSI*dPSI_drho2))
        + 4*(QU_PSI/PSI)*d2PSI_drho2drho2;
    d2QU_PSI_drho1drho2 = 12 * PSI * PSI * dPSI_drho1 * dPSI_drho2
        + 4 *(QU_PSI/PSI) * d2PSI_drho1drho2;

    { dQU_PSI_dx = dQU_PSI_drho1 * drho1_dx + dQU_PSI_drho2 * drho2_dx;
      dQU_PSI_dy = dQU_PSI_drho1 * drho1_dy + dQU_PSI_drho2 * drho2_dy;
      dQU_PSI_dz = dQU_PSI_drho1 * drho1_dz + dQU_PSI_drho2 * drho2_dz;
    } ;
    {d2QU_PSI_dxdx = d2QU_PSI_drho1drho1*drho1_dx*drho1_dx +
          dQU_PSI_drho1*d2rho1_dxdx + d2QU_PSI_drho2drho2*drho2_dx*drho2_dx +
          dQU_PSI_drho2*d2rho2_dxdx + 2 * d2QU_PSI_drho1drho2 * drho1_dx *
          drho2_dx;
      d2QU_PSI_dydy = d2QU_PSI_drho1drho1*drho1_dy*drho1_dy +
          dQU_PSI_drho1*d2rho1_dydy + d2QU_PSI_drho2drho2*drho2_dy*drho2_dy +
          dQU_PSI_drho2*d2rho2_dydy + 2 * d2QU_PSI_drho1drho2 * drho1_dy *
          drho2_dy;
      d2QU_PSI_dzdz = d2QU_PSI_drho1drho1*drho1_dz*drho1_dz +
          dQU_PSI_drho1*d2rho1_dzdz + d2QU_PSI_drho2drho2*drho2_dz*drho2_dz +
          dQU_PSI_drho2*d2rho2_dzdz + 2 * d2QU_PSI_drho1drho2 * drho1_dz *
          drho2_dz;
      d2QU_PSI_dydx = d2QU_PSI_dxdy = d2QU_PSI_drho1drho1 * drho1_dx * drho1_dy
          + d2QU_PSI_drho2drho2 * drho2_dx * drho2_dy + d2QU_PSI_drho1drho2 *
          (drho1_dx*drho2_dy+drho2_dx*drho1_dy) + dQU_PSI_drho1 * d2rho1_dxdy +
          dQU_PSI_drho2 * d2rho2_dxdy;
      d2QU_PSI_dzdy = d2QU_PSI_dydz = d2QU_PSI_drho1drho1 * drho1_dy * drho1_dz
          + d2QU_PSI_drho2drho2 * drho2_dy * drho2_dz + d2QU_PSI_drho1drho2 *
          (drho1_dy*drho2_dz+drho2_dy*drho1_dz) + dQU_PSI_drho1 * d2rho1_dydz +
          dQU_PSI_drho2 * d2rho2_dydz;
      d2QU_PSI_dxdz = d2QU_PSI_dzdx = d2QU_PSI_drho1drho1 * drho1_dz * drho1_dx
           + d2QU_PSI_drho2drho2 * drho2_dz * drho2_dx  + d2QU_PSI_drho1drho2 *
          (drho1_dz*drho2_dx+drho2_dz*drho1_dx)  + dQU_PSI_drho1 * d2rho1_dxdz
          + dQU_PSI_drho2 * d2rho2_dxdz;
    } ;


    psi1 = 1 + M1/(2*rho1) + M2*M1*M1/(2*S*(M1*M1+4*rho1*rho1));
    dpsi1_drho1 = -M1/(2*rho1*rho1)
        - M2*M1*M1*8*rho1/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1+4*rho1*rho1)) );
    dpsi1_drho2 = 0;
    d2psi1_drho1drho1 = 2*M1/(2*rho1*rho1*rho1)
        - M2*M1*M1*8/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1+4*rho1*rho1)) )
        + M2*M1*M1*8*2*4*2*rho1*rho1/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1
        +4*rho1*rho1)*(M1*M1+4*rho1*rho1)) );
    d2psi1_drho2drho2 = 0;
    d2psi1_drho1drho2 = 0;


    qu_psi1 = ((psi1)*(psi1)*(psi1)*(psi1)) ;
    dqu_psi1_drho1 = 4*(qu_psi1/psi1)*dpsi1_drho1;
    dqu_psi1_drho2 = 4*(qu_psi1/psi1)*dpsi1_drho2;

    d2qu_psi1_drho1drho1 = 12*((psi1*dpsi1_drho1)*(psi1*dpsi1_drho1))
        + 4*(qu_psi1/psi1)*d2psi1_drho1drho1;
    d2qu_psi1_drho2drho2 = 12*((psi1*dpsi1_drho2)*(psi1*dpsi1_drho2))
        + 4*(qu_psi1/psi1)*d2psi1_drho2drho2;
    d2qu_psi1_drho1drho2 = 12 * psi1 * psi1 * dpsi1_drho1 * dpsi1_drho2
        + 4 *(qu_psi1/psi1) * d2psi1_drho1drho2;
    { dqu_psi1_dx = dqu_psi1_drho1 * drho1_dx + dqu_psi1_drho2 * drho2_dx;
      dqu_psi1_dy = dqu_psi1_drho1 * drho1_dy + dqu_psi1_drho2 * drho2_dy;
      dqu_psi1_dz = dqu_psi1_drho1 * drho1_dz + dqu_psi1_drho2 * drho2_dz;
    } ;


    psi2 = 1 + M2/(2*rho2) + M1*M2*M2/(2*S*(M2*M2+4*rho2*rho2));
    dpsi2_drho1 = 0;
    dpsi2_drho2 = -M2/(2*rho2*rho2)
        - M1*M2*M2*8*rho2/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2+4*rho2*rho2)) );
    d2psi2_drho1drho1 = 0;
    d2psi2_drho2drho2 = 2*M2/(2*rho2*rho2*rho2)
        - M1*M2*M2*8/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2+4*rho2*rho2)) )
        + M1*M2*M2*8*2*4*2*rho2*rho2/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2
        +4*rho2*rho2)*(M2*M2+4*rho2*rho2)) );
    d2psi2_drho1drho2 = 0;


    qu_psi2 = ((psi2)*(psi2)*(psi2)*(psi2)) ;
    dqu_psi2_drho1 = 4*(qu_psi2/psi2)*dpsi2_drho1;
    dqu_psi2_drho2 = 4*(qu_psi2/psi2)*dpsi2_drho2;

    d2qu_psi2_drho1drho1 = 12*((psi2*dpsi2_drho1)*(psi2*dpsi2_drho1))
        + 4*(qu_psi2/psi2)*d2psi2_drho1drho1;
    d2qu_psi2_drho2drho2 = 12*((psi2*dpsi2_drho2)*(psi2*dpsi2_drho2))
        + 4*(qu_psi2/psi2)*d2psi2_drho2drho2;
    d2qu_psi2_drho1drho2 = 12 * psi2 * psi2 * dpsi2_drho1 * dpsi2_drho2
        + 4 *(qu_psi2/psi2) * d2psi2_drho1drho2;
    { dqu_psi2_dx = dqu_psi2_drho1 * drho1_dx + dqu_psi2_drho2 * drho2_dx;
      dqu_psi2_dy = dqu_psi2_drho1 * drho1_dy + dqu_psi2_drho2 * drho2_dy;
      dqu_psi2_dz = dqu_psi2_drho1 * drho1_dz + dqu_psi2_drho2 * drho2_dz;
    } ;


    {
      double f, df_drho1, df_drho2, d2f_drho1drho1, d2f_drho1drho2,
          d2f_drho2drho2 ;
      double g, dg_drho1, dg_drho2, d2g_drho1drho1, d2g_drho1drho2,
          d2g_drho2drho2 ;
      double h, dh_drho1, dh_drho2, d2h_drho1drho1, d2h_drho1drho2,
          d2h_drho2drho2 ;

      g = 1+M2/(2*rho2);
      dg_drho2 = -M2/(2*rho2*rho2);
      d2g_drho2drho2 = M2/(rho2*rho2*rho2);

      h = M1/(2*rho1*g);
      dh_drho1 = -h/rho1;
      d2h_drho1drho1 = 2*h/(rho1*rho1);
      dh_drho2 = -h*dg_drho2/g;
      d2h_drho2drho2 = h*(-d2g_drho2drho2
          + 2*dg_drho2*dg_drho2/g)/g;
      d2h_drho1drho2 = -dh_drho1*dg_drho2/g;

      a1 = (1-h)/(1+h);
      da1_drho1 = -2*dh_drho1/((1+h)*(1+h)) ;
      d2a1_drho1drho1 = -2*d2h_drho1drho1/((1+h)*(1+h))
          +4*((dh_drho1)*(dh_drho1)) /((1+h)*(1+h)*(1+h)) ;
      da1_drho2 = -2*dh_drho2/((1+h)*(1+h)) ;
      d2a1_drho2drho2 = -2*d2h_drho2drho2/((1+h)*(1+h))
          +4*((dh_drho2)*(dh_drho2)) /((1+h)*(1+h)*(1+h)) ;
      d2a1_drho1drho2 = -2*d2h_drho1drho2/((1+h)*(1+h))
          +4*dh_drho1*dh_drho2/((1+h)*(1+h)*(1+h)) ;

      if (0) {

        f = 1.0;
        df_drho1 = d2f_drho1drho1 = df_drho2
            = d2f_drho2drho2 = d2f_drho1drho2 = 0.0;
      } else {
        f = 1 + M2*M1*M1/(2*S*(M1*M1+4*rho1*rho1));
        df_drho1 = - M2*M1*M1*8*rho1/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1
            +4*rho1*rho1)) );
        d2f_drho1drho1 = - M2*M1*M1*8/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1
            +4*rho1*rho1)) )
            + M2*M1*M1*8*2*8*rho1*rho1/(2*S*((M1*M1+4*rho1*rho1)*(M1*M1
            +4*rho1*rho1)*(M1*M1+4*rho1*rho1)) );
        df_drho2 = d2f_drho2drho2 = d2f_drho1drho2 = 0;
      }

      alpha1 = a1*f*f;
      dalpha1_drho1 = da1_drho1*f*f + 2*a1*f*df_drho1;
      d2alpha1_drho1drho1 = d2a1_drho1drho1*f*f + 4*f*da1_drho1*df_drho1
          + 2*a1*df_drho1*df_drho1 + 2*f*a1*d2f_drho1drho1;

      dalpha1_drho2 = da1_drho2*f*f + 2*a1*f*df_drho2;
      d2alpha1_drho2drho2 = d2a1_drho2drho2*f*f + 4*f*da1_drho2*df_drho2
          + 2*a1*df_drho2*df_drho2 + 2*f*a1*d2f_drho2drho2;

      d2alpha1_drho1drho2 = d2a1_drho1drho2*f*f + 2*f*da1_drho1*df_drho2
          + 2*f*da1_drho2*df_drho1 + 2*a1*df_drho1*df_drho2
          + 2*f*a1*d2f_drho1drho2;

      sq_alpha1 = ((alpha1)*(alpha1)) ;
      dsq_alpha1_drho1 = 2*alpha1*dalpha1_drho1;
      d2sq_alpha1_drho1drho1 = 2*((dalpha1_drho1)*(dalpha1_drho1))
          + 2*alpha1*d2alpha1_drho1drho1;
      dsq_alpha1_drho2 = 2*alpha1*dalpha1_drho2;
      d2sq_alpha1_drho2drho2 = 2*((dalpha1_drho2)*(dalpha1_drho2))
          + 2*alpha1*d2alpha1_drho2drho2;
      d2sq_alpha1_drho1drho2 = 2*dalpha1_drho1*dalpha1_drho2
          + 2*alpha1*d2alpha1_drho1drho2;
    }


    {
      double f, df_drho1, df_drho2, d2f_drho1drho1, d2f_drho1drho2,
          d2f_drho2drho2 ;
      double g, dg_drho1, dg_drho2, d2g_drho1drho1, d2g_drho1drho2,
          d2g_drho2drho2 ;
      double h, dh_drho1, dh_drho2, d2h_drho1drho1, d2h_drho1drho2,
          d2h_drho2drho2 ;

      g = 1+M1/(2*rho1);
      dg_drho1 = -M1/(2*rho1*rho1);
      d2g_drho1drho1 = M1/(rho1*rho1*rho1);

      h = M2/(2*rho2*g);
      dh_drho2 = -h/rho2;
      d2h_drho2drho2 = 2*h/(rho2*rho2);
      dh_drho1 = -h*dg_drho1/g;
      d2h_drho1drho1 = h*(-d2g_drho1drho1
          + 2*dg_drho1*dg_drho1/g)/g;
      d2h_drho1drho2 = -dh_drho2*dg_drho1/g;

      a2 = (1-h)/(1+h);
      da2_drho1 = -2*dh_drho1/((1+h)*(1+h)) ;
      d2a2_drho1drho1 = -2*d2h_drho1drho1/((1+h)*(1+h))
          +4*((dh_drho1)*(dh_drho1)) /((1+h)*(1+h)*(1+h)) ;
      da2_drho2 = -2*dh_drho2/((1+h)*(1+h)) ;
      d2a2_drho2drho2 = -2*d2h_drho2drho2/((1+h)*(1+h))
          +4*((dh_drho2)*(dh_drho2)) /((1+h)*(1+h)*(1+h)) ;
      d2a2_drho1drho2 = -2*d2h_drho1drho2/((1+h)*(1+h))
          +4*dh_drho1*dh_drho2/((1+h)*(1+h)*(1+h)) ;

      if (0) {

        f = 1.0;
        df_drho1 = d2f_drho1drho1 = df_drho2
            = d2f_drho2drho2 = d2f_drho1drho2 = 0.0;
      } else {
        f = 1 + M1*M2*M2/(2*S*(M2*M2+4*rho2*rho2));
        df_drho2 = - M1*M2*M2*8*rho2/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2
            +4*rho2*rho2)) );
        d2f_drho2drho2 = - M1*M2*M2*8/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2
            +4*rho2*rho2)) )
            + M1*M2*M2*8*2*8*rho2*rho2/(2*S*((M2*M2+4*rho2*rho2)*(M2*M2
            +4*rho2*rho2)*(M2*M2+4*rho2*rho2)) );
        df_drho1 = d2f_drho1drho1 = d2f_drho1drho2 = 0;
      }

      alpha2 = a2*f*f;
      dalpha2_drho1 = da2_drho1*f*f + 2*a2*f*df_drho1;
      d2alpha2_drho1drho1 = d2a2_drho1drho1*f*f + 4*f*da2_drho1*df_drho1
          + 2*a2*df_drho1*df_drho1 + 2*f*a2*d2f_drho1drho1;

      dalpha2_drho2 = da2_drho2*f*f + 2*a2*f*df_drho2;
      d2alpha2_drho2drho2 = d2a2_drho2drho2*f*f + 4*f*da2_drho2*df_drho2
          + 2*a2*df_drho2*df_drho2 + 2*f*a2*d2f_drho2drho2;

      d2alpha2_drho1drho2 = d2a2_drho1drho2*f*f + 2*f*da2_drho1*df_drho2
          + 2*f*da2_drho2*df_drho1 + 2*a2*df_drho1*df_drho2
          + 2*f*a2*d2f_drho1drho2;

      sq_alpha2 = ((alpha2)*(alpha2)) ;
      dsq_alpha2_drho1 = 2*alpha2*dalpha2_drho1;
      d2sq_alpha2_drho1drho1 = 2*((dalpha2_drho1)*(dalpha2_drho1))
          + 2*alpha2*d2alpha2_drho1drho1;
      dsq_alpha2_drho2 = 2*alpha2*dalpha2_drho2;
      d2sq_alpha2_drho2drho2 = 2*((dalpha2_drho2)*(dalpha2_drho2))
          + 2*alpha2*d2alpha2_drho2drho2;
      d2sq_alpha2_drho1drho2 = 2*dalpha2_drho1*dalpha2_drho2
          + 2*alpha2*d2alpha2_drho1drho2;
    }


    { dalpha1_dx = dalpha1_drho1 * drho1_dx + dalpha1_drho2 * drho2_dx;
      dalpha1_dy = dalpha1_drho1 * drho1_dy + dalpha1_drho2 * drho2_dy;
      dalpha1_dz = dalpha1_drho1 * drho1_dz + dalpha1_drho2 * drho2_dz;
    } ;
    { dalpha2_dx = dalpha2_drho1 * drho1_dx + dalpha2_drho2 * drho2_dx;
      dalpha2_dy = dalpha2_drho1 * drho1_dy + dalpha2_drho2 * drho2_dy;
      dalpha2_dz = dalpha2_drho1 * drho1_dz + dalpha2_drho2 * drho2_dz;
    } ;
    { dsq_alpha1_dx = dsq_alpha1_drho1 * drho1_dx + dsq_alpha1_drho2 * drho2_dx;
      dsq_alpha1_dy = dsq_alpha1_drho1 * drho1_dy + dsq_alpha1_drho2 * drho2_dy;
      dsq_alpha1_dz = dsq_alpha1_drho1 * drho1_dz + dsq_alpha1_drho2 * drho2_dz;
    } ;
    { dsq_alpha2_dx = dsq_alpha2_drho1 * drho1_dx + dsq_alpha2_drho2 * drho2_dx;
      dsq_alpha2_dy = dsq_alpha2_drho1 * drho1_dy + dsq_alpha2_drho2 * drho2_dy;
      dsq_alpha2_dz = dsq_alpha2_drho1 * drho1_dz + dsq_alpha2_drho2 * drho2_dz;
    } ;


    common = QU_PSI - U*U*(sq_alpha1-qu_psi1)/(1-U*U)
        - V*V*(sq_alpha2-qu_psi2)/(1-V*V);

    dcommon_drho1 = dQU_PSI_drho1
        - (U*U/(1-U*U))*(dsq_alpha1_drho1 - dqu_psi1_drho1)
        - (V*V/(1-V*V))*(dsq_alpha2_drho1 - dqu_psi2_drho1);

    dcommon_drho2 = dQU_PSI_drho2
        - (U*U/(1-U*U))*(dsq_alpha1_drho2 - dqu_psi1_drho2)
        - (V*V/(1-V*V))*(dsq_alpha2_drho2 - dqu_psi2_drho2);

    d2common_drho1drho1 = d2QU_PSI_drho1drho1
        - (U*U/(1-U*U))*(d2sq_alpha1_drho1drho1 - d2qu_psi1_drho1drho1)
        - (V*V/(1-V*V))*(d2sq_alpha2_drho1drho1 - d2qu_psi2_drho1drho1);

    d2common_drho2drho2 = d2QU_PSI_drho2drho2
        - (U*U/(1-U*U))*(d2sq_alpha1_drho2drho2 - d2qu_psi1_drho2drho2)
        - (V*V/(1-V*V))*(d2sq_alpha2_drho2drho2 - d2qu_psi2_drho2drho2);

    d2common_drho1drho2 = d2QU_PSI_drho1drho2
        - (U*U/(1-U*U))*(d2sq_alpha1_drho1drho2 - d2qu_psi1_drho1drho2)
        - (V*V/(1-V*V))*(d2sq_alpha2_drho1drho2 - d2qu_psi2_drho1drho2);


    { dcommon_dx = dcommon_drho1 * drho1_dx + dcommon_drho2 * drho2_dx;
      dcommon_dy = dcommon_drho1 * drho1_dy + dcommon_drho2 * drho2_dy;
      dcommon_dz = dcommon_drho1 * drho1_dz + dcommon_drho2 * drho2_dz;
    } ;
    {d2common_dxdx = d2common_drho1drho1*drho1_dx*drho1_dx +
          dcommon_drho1*d2rho1_dxdx + d2common_drho2drho2*drho2_dx*drho2_dx +
          dcommon_drho2*d2rho2_dxdx + 2 * d2common_drho1drho2 * drho1_dx *
          drho2_dx;
      d2common_dydy = d2common_drho1drho1*drho1_dy*drho1_dy +
          dcommon_drho1*d2rho1_dydy + d2common_drho2drho2*drho2_dy*drho2_dy +
          dcommon_drho2*d2rho2_dydy + 2 * d2common_drho1drho2 * drho1_dy *
          drho2_dy;
      d2common_dzdz = d2common_drho1drho1*drho1_dz*drho1_dz +
          dcommon_drho1*d2rho1_dzdz + d2common_drho2drho2*drho2_dz*drho2_dz +
          dcommon_drho2*d2rho2_dzdz + 2 * d2common_drho1drho2 * drho1_dz *
          drho2_dz;
      d2common_dydx = d2common_dxdy = d2common_drho1drho1 * drho1_dx * drho1_dy
          + d2common_drho2drho2 * drho2_dx * drho2_dy + d2common_drho1drho2 *
          (drho1_dx*drho2_dy+drho2_dx*drho1_dy) + dcommon_drho1 * d2rho1_dxdy +
          dcommon_drho2 * d2rho2_dxdy;
      d2common_dzdy = d2common_dydz = d2common_drho1drho1 * drho1_dy * drho1_dz
          + d2common_drho2drho2 * drho2_dy * drho2_dz + d2common_drho1drho2 *
          (drho1_dy*drho2_dz+drho2_dy*drho1_dz) + dcommon_drho1 * d2rho1_dydz +
          dcommon_drho2 * d2rho2_dydz;
      d2common_dxdz = d2common_dzdx = d2common_drho1drho1 * drho1_dz * drho1_dx
           + d2common_drho2drho2 * drho2_dz * drho2_dx  + d2common_drho1drho2 *
          (drho1_dz*drho2_dx+drho2_dz*drho1_dx)  + dcommon_drho1 * d2rho1_dxdz
          + dcommon_drho2 * d2rho2_dxdz;
    } ;

    lapse = a1 * a2 * ( (((1+M2/(2*S))*(1+M2/(2*S*(1+M1/(2*rho1)))))*((1
        +M2/(2*S))*(1+M2/(2*S*(1+M1/(2*rho1))))))
        * M1*M1*rho2*rho2/(M1*M1*rho2*rho2+4*rho1*rho1*rho1*rho1)
        + (((1+M1/(2*S))*(1+M1/(2*S*(1+M2/(2*rho2)))))*((1+M1/(2*S))*(1
        +M1/(2*S*(1+M2/(2*rho2))))))
        * M2*M2*rho1*rho1/(M2*M2*rho1*rho1+4*rho2*rho2*rho2*rho2)
        + (1 - M1*M1*rho2*rho2/(M1*M1*rho2*rho2+4*rho1*rho1*rho1*rho1)
        - M2*M2*rho1*rho1/(M2*M2*rho1*rho1+4*rho2*rho2*rho2*rho2))
        ) / sqrt(common/QU_PSI);


    {gxx = QU_PSI;
      dgxx_drho1 = dQU_PSI_drho1;
      dgxx_drho2 = dQU_PSI_drho2;
      d2gxx_drho1drho1 = d2QU_PSI_drho1drho1;
      d2gxx_drho2drho2 = d2QU_PSI_drho2drho2;
      d2gxx_drho1drho2 = d2QU_PSI_drho1drho2;
      dgxx_dx = dQU_PSI_dx;
      dgxx_dy = dQU_PSI_dy;
      dgxx_dz = dQU_PSI_dz;
      d2gxx_dxdx = d2QU_PSI_dxdx;
      d2gxx_dxdy = d2gxx_dydx = d2QU_PSI_dxdy;
      d2gxx_dzdx = d2gxx_dxdz = d2QU_PSI_dxdz;
      d2gxx_dydy = d2QU_PSI_dydy;
      d2gxx_dzdy = d2gxx_dydz = d2QU_PSI_dydz;
      d2gxx_dzdz = d2QU_PSI_dzdz;
    } ;
    {gyy = common;
      dgyy_drho1 = dcommon_drho1;
      dgyy_drho2 = dcommon_drho2;
      d2gyy_drho1drho1 = d2common_drho1drho1;
      d2gyy_drho2drho2 = d2common_drho2drho2;
      d2gyy_drho1drho2 = d2common_drho1drho2;
      dgyy_dx = dcommon_dx;
      dgyy_dy = dcommon_dy;
      dgyy_dz = dcommon_dz;
      d2gyy_dxdx = d2common_dxdx;
      d2gyy_dxdy = d2gyy_dydx = d2common_dxdy;
      d2gyy_dzdx = d2gyy_dxdz = d2common_dxdz;
      d2gyy_dydy = d2common_dydy;
      d2gyy_dzdy = d2gyy_dydz = d2common_dydz;
      d2gyy_dzdz = d2common_dzdz;
    } ;
    {gzz = QU_PSI;
      dgzz_drho1 = dQU_PSI_drho1;
      dgzz_drho2 = dQU_PSI_drho2;
      d2gzz_drho1drho1 = d2QU_PSI_drho1drho1;
      d2gzz_drho2drho2 = d2QU_PSI_drho2drho2;
      d2gzz_drho1drho2 = d2QU_PSI_drho1drho2;
      dgzz_dx = dQU_PSI_dx;
      dgzz_dy = dQU_PSI_dy;
      dgzz_dz = dQU_PSI_dz;
      d2gzz_dxdx = d2QU_PSI_dxdx;
      d2gzz_dxdy = d2gzz_dydx = d2QU_PSI_dxdy;
      d2gzz_dzdx = d2gzz_dxdz = d2QU_PSI_dxdz;
      d2gzz_dydy = d2QU_PSI_dydy;
      d2gzz_dzdy = d2gzz_dydz = d2QU_PSI_dydz;
      d2gzz_dzdz = d2QU_PSI_dzdz;
    } ;


    gXX = 1/gxx;
    gYY = 1/gyy;
    gZZ = 1/gzz;

    { double f, df_drho1, df_drho2, df_dx, df_dy, df_dz;

      f = 1 + U*U*M1*M1*rho2*rho2/((1-U*U)*(M1*M1*rho2*rho2 +
          4*((rho1)*(rho1)*(rho1)*(rho1)) ))
          + V*V*M2*M2*rho1*rho1/((1-V*V)*(M2*M2*rho1*rho1 +
          4*((rho2)*(rho2)*(rho2)*(rho2)) ));

      df_drho1 = U*U/(1-U*U) * M1*M1*rho2*rho2/(M1*M1*rho2*rho2
          +4*((rho1)*(rho1)*(rho1)*(rho1)) )
          * (-16*((rho1)*(rho1)*(rho1)) /(M1*M1*rho2*rho2
          +4*((rho1)*(rho1)*(rho1)*(rho1)) ))
          + V*V/(1-V*V) * M2*M2*rho1*rho1/(M2*M2*rho1*rho1
          +4*((rho2)*(rho2)*(rho2)*(rho2)) )
          *(2/rho1 - 2*rho1*M2*M2/(M2*M2*rho1*rho1
          +4*((rho2)*(rho2)*(rho2)*(rho2)) ));

      df_drho2 = V*V/(1-V*V) * M2*M2*rho1*rho1/(M2*M2*rho1*rho1
          +4*((rho2)*(rho2)*(rho2)*(rho2)) )
          * (-16*((rho2)*(rho2)*(rho2)) /(M2*M2*rho1*rho1
          +4*((rho2)*(rho2)*(rho2)*(rho2)) ))
          + U*U/(1-U*U) * M1*M1*rho2*rho2/(M1*M1*rho2*rho2
          +4*((rho1)*(rho1)*(rho1)*(rho1)) )
          *(2/rho2 - 2*rho2*M1*M1/(M1*M1*rho2*rho2
          +4*((rho1)*(rho1)*(rho1)*(rho1)) ));

      { df_dx = df_drho1 * drho1_dx + df_drho2 * drho2_dx;
        df_dy = df_drho1 * drho1_dy + df_drho2 * drho2_dy;
        df_dz = df_drho1 * drho1_dz + df_drho2 * drho2_dz;
      } ;

      beta_X = -f*omega*y;
      beta_Y = omega*x +
          (U*(sq_alpha1 - qu_psi1)/(1-U*U) + (V*(sq_alpha2 - qu_psi2)/(1
          -V*V)))/gyy;
      beta_Z = 0.0;

      beta_x = gxx * beta_X;
      beta_y = gyy * beta_Y;
      beta_z = gzz * beta_Z;

      dbeta_x_dx = dgxx_dx * beta_X - df_dx * omega * y * gxx;
      dbeta_x_dy = dgxx_dy * beta_X - f * omega * gxx - df_dy * omega * y * gxx;
      dbeta_x_dz = dgxx_dz * beta_X - df_dz * omega * y * gxx;

      dbeta_y_dx = omega * x * dgyy_dx + omega * gyy
          + U*(dsq_alpha1_dx - dqu_psi1_dx)/(1-U*U)
          + V*(dsq_alpha2_dx - dqu_psi2_dx)/(1-V*V);
      dbeta_y_dy = dgyy_dy * omega * x
          + U*(dsq_alpha1_dy - dqu_psi1_dy)/(1-U*U)
          + V*(dsq_alpha2_dy - dqu_psi2_dy)/(1-V*V);
      dbeta_y_dz = dgyy_dz * omega * x
          + U*(dsq_alpha1_dz - dqu_psi1_dz)/(1-U*U)
          + V*(dsq_alpha2_dz - dqu_psi2_dz)/(1-V*V);
    }
  }


  {


    G_xxx =  dgxx_dx / 2;
    G_xxy =  dgxx_dy / 2;
    G_xyx = G_xxy;
    G_xxz =  dgxx_dz / 2;
    G_xzx = G_xxz;
    G_xyy = -dgyy_dx / 2;
    G_xyz =  0.0;
    G_xzy = 0.0;
    G_xzz = -dgzz_dx / 2;

    G_yxx = -dgxx_dy / 2;
    G_yxy =  dgyy_dx / 2;
    G_yyx = G_yxy;
    G_yxz =  0.0;
    G_yzx = 0.0;
    G_yyy =  dgyy_dy / 2;
    G_yyz =  dgyy_dz / 2;
    G_yzy = G_yyz;
    G_yzz = -dgzz_dy / 2;

    G_zxx = -dgxx_dz / 2;
    G_zxy =  0.0;
    G_zyx = 0.0;
    G_zxz =  dgzz_dx / 2;
    G_zzx = G_zxz;
    G_zyy = -dgyy_dz / 2;
    G_zyz =  dgzz_dy / 2;
    G_zzy = G_zyz;
    G_zzz =  dgzz_dz / 2;


    G_Xxx = gXX * G_xxx;
    G_Xxy = gXX * G_xxy;
    G_Xyx = G_Xxy;
    G_Xxz = gXX * G_xxz;
    G_Xzx = G_Xxz;
    G_Xyy = gXX * G_xyy;
    G_Xyz = gXX * G_xyz;
    G_Xzy = G_Xyz;
    G_Xzz = gXX * G_xzz;

    G_Yxx = gYY * G_yxx;
    G_Yxy = gYY * G_yxy;
    G_Yyx = G_Yxy;
    G_Yxz = gYY * G_yxz;
    G_Yzx = G_Yxz;
    G_Yyy = gYY * G_yyy;
    G_Yyz = gYY * G_yyz;
    G_Yzy = G_Yyz;
    G_Yzz = gYY * G_yzz;

    G_Zxx = gZZ * G_zxx;
    G_Zxy = gZZ * G_zxy;
    G_Zyx = G_Zxy;
    G_Zxz = gZZ * G_zxz;
    G_Zzx = G_Zxz;
    G_Zyy = gZZ * G_zyy;
    G_Zyz = gZZ * G_zyz;
    G_Zzy = G_Zyz;
    G_Zzz = gZZ * G_zzz;
  }


  {


    K1xx = U *alpha1*drho1_dy*dgxx_drho1   /(2*sqrt(QU_PSI*gyy));
    K1zz = U *alpha1*drho1_dy*dgzz_drho1   /(2*sqrt(QU_PSI*gyy));
    K2xx = V *alpha2*drho2_dy*dgxx_drho2   /(2*sqrt(QU_PSI*gyy));
    K2zz = V *alpha2*drho2_dy*dgzz_drho2   /(2*sqrt(QU_PSI*gyy));


    K1yy =    U *(-alpha1*drho1_dy*dgyy_drho1   + 4*gyy*drho1_dy*dalpha1_drho1)
          /(2*sqrt(QU_PSI*gyy)) ;
    K2yy =    V *(-alpha2*drho2_dy*dgyy_drho2   + 4*gyy*drho2_dy*dalpha2_drho2)
          /(2*sqrt(QU_PSI*gyy)) ;


    K1yx = K1xy =    U *sqrt(gyy/QU_PSI)*drho1_dx    * (dalpha1_drho1
        -alpha1*dgyy_drho1/(2*gyy));
    K1yz = K1zy =    U *sqrt(gyy/QU_PSI)*drho1_dz    * (dalpha1_drho1
        -alpha1*dgyy_drho1/(2*gyy));
    K2yx = K2xy =    V *sqrt(gyy/QU_PSI)*drho2_dx    * (dalpha2_drho2
        -alpha2*dgyy_drho2/(2*gyy));
    K2yz = K2zy =    V *sqrt(gyy/QU_PSI)*drho2_dz    * (dalpha2_drho2
        -alpha2*dgyy_drho2/(2*gyy));


    K_xx= K1xx+ K2xx;
    K_yy= K1yy+ K2yy;
    K_zz= K1zz+ K2zz;
    K_xy= K1xy+ K2xy;
    K_yz= K1yz+ K2yz;
    K_yx= K1yx+ K2yx;
    K_zy= K1zy+ K2zy;


    K_xz = 0;
    K_zx = 0;

    K = gXX * K_xx + gYY * K_yy + gZZ * K_zz;

    K_XX = gXX * gXX * K_xx;
    K_YY = gYY * gYY * K_yy;
    K_ZZ = gZZ * gZZ * K_zz;
    K_XY = gXX * gYY * K_xy;
    K_YX = K_XY;
    K_XZ = gXX * gZZ * K_xz;
    K_ZX = K_XZ;
    K_YZ = gYY * gZZ * K_yz;
    K_ZY = K_YZ;
    K_AB_K_ab = K_XX * K_xx + K_XY * K_xy + K_XZ * K_xz
        + K_YX * K_yx + K_YY * K_yy + K_YZ * K_yz
        + K_ZX * K_zx + K_ZY * K_zy + K_ZZ * K_zz;


    D_x_beta_y = dbeta_y_dx - G_Yxy*beta_y - G_Xxy*beta_x;
    D_y_beta_y = dbeta_y_dy - G_Yyy*beta_y - G_Xyy*beta_x;
    D_z_beta_y = dbeta_y_dz - G_Yzy*beta_y - G_Xzy*beta_x;
    D_x_beta_x = dbeta_x_dx - G_Yxx*beta_y - G_Xxx*beta_x;
    D_y_beta_x = dbeta_x_dy - G_Yyx*beta_y - G_Xyx*beta_x;
    D_z_beta_x = dbeta_x_dz - G_Yzx*beta_y - G_Xzx*beta_x;
    D_x_beta_z =  - G_Yxz*beta_y - G_Xxz*beta_x;
    D_y_beta_z =  - G_Yyz*beta_y - G_Xyz*beta_x;
    D_z_beta_z =  - G_Yzz*beta_y - G_Xzz*beta_x;

  }


  {
    R = - d2gxx_dydy / (gxx * gyy) - d2gxx_dzdz / (gxx * gzz)
        - d2gyy_dxdx / (gxx * gyy) - d2gyy_dzdz / (gyy * gzz)
        - d2gzz_dxdx / (gxx * gzz) - d2gzz_dydy / (gyy * gzz)

    + (dgxx_dy)*(dgxx_dy) / (2*gxx*gxx*gyy) + (dgxx_dz)*(dgxx_dz) /
        (2*gxx*gxx*gzz)
        + (dgyy_dx)*(dgyy_dx) / (2*gxx*gyy*gyy) + (dgyy_dz)*(dgyy_dz) /
        (2*gyy*gyy*gzz)
        + (dgzz_dx)*(dgzz_dx) / (2*gxx*gzz*gzz) + (dgzz_dy)*(dgzz_dy) /
        (2*gyy*gzz*gzz)

    + (dgxx_dx)*(dgyy_dx) / (2*gxx*gxx*gyy) + (dgxx_dx)*(dgzz_dx) /
        (2*gxx*gxx*gzz)

    + (dgxx_dy)*(dgyy_dy) / (2*gxx*gyy*gyy) - (dgxx_dy)*(dgzz_dy) /
        (2*gxx*gyy*gzz)

    - (dgxx_dz)*(dgyy_dz) / (2*gxx*gyy*gzz) + (dgxx_dz)*(dgzz_dz) /
        (2*gxx*gzz*gzz)

    - (dgyy_dx)*(dgzz_dx) / (2*gxx*gyy*gzz)

    + (dgyy_dy)*(dgzz_dy) / (2*gyy*gyy*gzz)

    + (dgyy_dz)*(dgzz_dz) / (2*gyy*gzz*gzz);
  }


  {
    double D_b_K_XB,D_b_K_YB,D_b_K_ZB,D_x_K_XX,D_y_K_XY,D_z_K_XZ,
        D_x_K_YX,D_y_K_YY,D_z_K_YZ,D_x_K_ZX,D_y_K_ZY,D_z_K_ZZ,
        D_X_K,D_Y_K,D_Z_K,measure,scalar;


    scalar = lapse * (R + K*K - K_AB_K_ab)
        + 2 * (  D_x_beta_x * (K_XX - gXX*K)
        +  D_y_beta_y * (K_YY - gYY*K)
        +  D_z_beta_z * (K_ZZ - gZZ*K)
        + (D_x_beta_y + D_y_beta_x) * K_XY
        + (D_x_beta_z + D_z_beta_x) * K_XZ
        + (D_z_beta_y + D_y_beta_z) * K_ZY);


    measure = sqrt( gxx * gyy * gzz );
    if (0) {
      printf(" x %g, y %g, z %g\n", x, y, z);
      printf("  R  %g    measure               %g\n", R, measure);
      printf("  lapse (R + K*K - K_AB_K_ab)    %g   %g\n",
          lapse, (R + K*K - K_AB_K_ab));
      printf("D_x_beta_x * (K_XX - gXX*K) %g\n",
          D_x_beta_x * (K_XX - gXX*K));
      printf("D_y_beta_y * (K_YY - gYY*K) %g\n",
          D_y_beta_y * (K_YY - gYY*K));
      printf("D_z_beta_z * (K_ZZ - gZZ*K) %g\n",
          D_z_beta_z * (K_ZZ - gZZ*K));
      printf("(D_x_beta_y + D_y_beta_x) * K_XY %g \n",
          (D_x_beta_y + D_y_beta_x) * K_XY);
      printf("(D_x_beta_z + D_z_beta_x) * K_XZ %g \n",
          (D_x_beta_z + D_z_beta_x) * K_XZ);
      printf("(D_z_beta_y + D_y_beta_z) * K_ZY %g \n",
          (D_z_beta_y + D_y_beta_z) * K_ZY);

      printf("  returning integrand %g   NO r^2 sin(theta)\n", scalar *
          measure);
    }
    return( (measure * scalar) );
  }


}

