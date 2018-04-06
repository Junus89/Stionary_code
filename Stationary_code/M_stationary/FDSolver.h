#ifndef FDSolver_H_
#include<complex.h> 
#include<math.h>
#define FDSolver_H_
//#define PI atan(1,1);
#define PI	3.14159265358979323846


//double Tint,TR,fR,OmegaR;
double DT;
int DSNum, BNum, FNum, TNum;
double R,Rou_0,C_0,RMaTip,Omega;
//double VectorX,VectorY,VectorZ,VuN,VvN,VUN;
double OX,OY,OZ, MaX, MaY, MaZ;
double Time[64];

double OmegaM, OmegaR, Gamma;
double complex z1 = 0.0 + 1.0*I;


/*---------- Memory allocation function declaretaion -----------------*/
float *vector(long nl,long nh);
void free_vector(float *inputAr,long nl,long nh);

double *dvector(long nl,long nh);
void free_dvector(double *inputAr,long nl,long nh);

float **matrix(long nrl,long nrh, long ncl,long nch);
void free_matrix(float **inputAr,long nrl,long nrh,long ncl,long nch);

double **dmatrix(long nrl,long nrh, long ncl,long nch);
void free_dmatrix(double **inputAr,long nrl,long nrh,long ncl,long nch);

float ***f3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh);
void free_f3dmatrix(float ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh);

double ***d3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh);
void free_d3dmatrix(double ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh);

float ****f4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh);
void free_f4dmatrix(float ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh);

double ****d4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh);
void free_d4dmatrix(double ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh);












// ----fifthloop prepreation-----



//void fifthLoop();
//void fifthLoop(double OmegaR, double Omega, double MaX, double MaY, double MaZ, double ka, double **pF,double complex SP1,double complex SP2,double complex SP3);
void fifthLoop(double OmegaR, double Omega, double MaX, double MaY, double MaZ, double ka);


double *Ones(int n);// Ones function


#endif


