/*
C FUNCTIONS - DESCRIPTIONS TO BE ADDED
TO COMPILE: gcc -shared -fPIC -o C_Kernel_WAD-jack.so C_Kernel_WAD-jack.c -lm
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <limits.h>
# include <float.h>

# define EPS DBL_EPSILON

double phi(double x) {return(exp(-0.5*pow(x,2))/sqrt(2.*M_PI));};
double Dphi(double x) {return(-x*phi(x));};

/*
double w(double x, double trim) {return(exp(-pow(x,2)/pow(pow(trim,2)-pow(x,2),1/4)));};
double Dw(double x, double trim) {return((-2.*x/pow(pow(trim,2)-pow(x,2),1/4) -0.5*pow(x,3)/pow(pow(trim,2)-pow(x,2),5/4))*w(x,trim));};
*/

double w(double x, double trim) {return(exp(-pow(x,4)/(pow(trim,4)*(pow(trim,4)-pow(x,4)))));};
double Dw(double x, double trim) {return(4*( -pow(x,7)/(pow(trim,4)*pow(pow(trim,4)-pow(x,4),2)) - pow(x,3)/(pow(trim,4)*(pow(trim,4)-pow(x,4))) ) * w(x,trim));};

double K(double x, int P){
    double Pol = 0.;
    switch (P) {
    case 2: 
        Pol=1.; break;
    case 4:
        Pol=(3.-pow(x,2))/2.; break;
    case 6:
        Pol=(pow(x,4)-10.*pow(x,2)+15.)/8.; break;
    case 8:
        Pol=(-pow(x,6)+21.*pow(x,4)-105.*pow(x,2)+105.)/48.; break;
    case 10:
        Pol=(pow(x,8)-36.*pow(x,6)+378.*pow(x,4)-1260.*pow(x,2)+945.)/384.; break;
    };
    return(Pol * phi(x));
};

double DK(double x, int P){
    double Pol = 0.;
    switch (P) {
    case 2: 
        Pol=1.; break;
    case 4:
        Pol=(5.-pow(x,2))/2.; break;
    case 6:
        Pol=(pow(x,4)-14.*pow(x,2)+35.)/8.; break;
    case 8:
        Pol=(-pow(x,6)+27.*pow(x,4)-189.*pow(x,2)+315.)/48.; break;
    case 10:
        Pol=(pow(x,8)-44.*pow(x,6)+594.*pow(x,4)-2772.*pow(x,2)+3465.)/384.; break;
    };
    return(Pol * Dphi(x));
};

int v_ij(int i, int j, int n) {
/*    if (i<j) {return(i*n-(i*(i+1))/2+j-i);} else {return(j*n-(j*(j+1))/2+i-j);};*/
    return(i*n-(i*(i+1))/2+j-i-1);
};

void wad(double *theta1_hat_LI,
         double *sigma1_hat_LI,
         double *K_ij, double *DK_ij, double *w_i, double *Dw_i, int *used,
         double *f_LI_i, double *Df_LI_i,
         double *e_LI_i, double *De_LI_i,
         double *y, double *x, int *d, double *trim,
         int *P, int *n, double *h) { 

    int i, j, k; double temp, s_LI_i;

    double K_ii  = pow(K(0.,P[0]),d[0]); double DK_ii = DK(0.,P[0])*pow(K(0.,P[0]),d[0]-1);

    double theta2_hat_LI=0., theta3_hat_LI=0.;

    for (i=0; i<n[0]; i++) {

        for (k=0; k<d[0]; k++) if (-trim[k]<x[k*n[0]+i] && x[k*n[0]+i]<trim[k]) used[i] *= 1; else used[i] *= 0;

        if (used[i] == 1) {
            temp = 1.; for (k=1; k<d[0]; k++) {temp *= w(x[k*n[0]+i],trim[k]);}
            w_i[i] = w(x[i],trim[0])*temp; Dw_i[i] = Dw(x[i],trim[0])*temp;
        }

        for (j=i+1; j<n[0]; j++) {
            temp = 1.; for (k=1; k<d[0]; k++) {temp *= K((x[k*n[0]+i]-x[k*n[0]+j])/h[k],P[0]);}
            K_ij[v_ij(i,j,n[0])] = K((x[i]-x[j])/h[0],P[0])*temp; DK_ij[v_ij(i,j,n[0])] = DK((x[i]-x[j])/h[0],P[0])*temp;
        }
    }

    theta2_hat_LI=0.; theta3_hat_LI=0.;

    for (i=0; i<n[0]; i++) if (used[i] == 1) {
        f_LI_i[i]=0.; Df_LI_i[i]=0.; e_LI_i[i]=0.; De_LI_i[i]=0.;

        for (j=0; j<n[0]; j++) {
            if (i < j) {
                f_LI_i[i] += K_ij[v_ij(i,j,n[0])]      ; Df_LI_i[i] += DK_ij[v_ij(i,j,n[0])];
                e_LI_i[i] += y[j]*K_ij[v_ij(i,j,n[0])] ; De_LI_i[i] += y[j]*DK_ij[v_ij(i,j,n[0])];
            } else if (i == j) {
                f_LI_i[i] += K_ii      ; Df_LI_i[i] += DK_ii;
                e_LI_i[i] += y[i]*K_ii ; De_LI_i[i] += y[i]*DK_ii;
            } else if (i > j) {
                f_LI_i[i] += K_ij[v_ij(j,i,n[0])]      ; Df_LI_i[i] += -DK_ij[v_ij(j,i,n[0])];
                e_LI_i[i] += y[j]*K_ij[v_ij(j,i,n[0])] ; De_LI_i[i] += y[j]*(-DK_ij[v_ij(j,i,n[0])]);
            }
        }

        s_LI_i = -Dw_i[i] - w_i[i] * Df_LI_i[i]/(h[0]*f_LI_i[i]);

        theta1_hat_LI[0] += y[i] * s_LI_i;

        sigma1_hat_LI[0] += pow((w_i[i]*De_LI_i[i]/h[0] + Dw_i[i]*e_LI_i[i])/f_LI_i[i] + y[i]*s_LI_i,2);

        theta2_hat_LI += w_i[i]*(De_LI_i[i]/f_LI_i[i] - e_LI_i[i]*Df_LI_i[i]/pow(f_LI_i[i],2));

        theta3_hat_LI += e_LI_i[i]/f_LI_i[i] * s_LI_i;
    }

    theta1_hat_LI[0] /= n[0]; theta2_hat_LI /= n[0]*h[0]; theta3_hat_LI /= n[0];

    sigma1_hat_LI[0] = sigma1_hat_LI[0]/n[0] - pow(theta1_hat_LI[0],2)
                       - 2.*theta1_hat_LI[0]*theta2_hat_LI + 2.*theta1_hat_LI[0]*theta3_hat_LI;

};


