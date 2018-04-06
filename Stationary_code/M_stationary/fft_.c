/************************************************
* FFT code from the book Numerical Recipes in C *
* Visit www.nr.com for the licence.             *
************************************************/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

/*
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

 Inputs:
	data[] : array of complex* data points of size 2*NFFT+1.
		data[0] is unused,
		* the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
			data[2*n+1] = real(x(n))
			data[2*n+2] = imag(x(n))
		if length(Nx) < NFFT, the remainder of the array must be padded with zeros

	nn : FFT order NFFT. This MUST be a power of 2 and >= length(x).
	isign:  if set to 1, 
				computes the forward FFT
			if set to -1, 
				computes Inverse FFT - in this case the output values have
				to be manually normalized by multiplying with 1/NFFT.
 Outputs:
	data[] : The FFT or IFFT results are stored in data, overwriting the input.
*/

void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    //for (int i = 1; i < n; i += 2) {
	for (int i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

/********************************************************
* The following is a test routine that generates a ramp *
* with 10 elements, finds their FFT, and then finds the *
* original sequence using inverse FFT                   *
********************************************************/

int main(int argc, char * argv[])
{
	int i;
	int Nx=0;
	int NFFT;
	double *x;
	double *X;

	/* generate a ramp with 10 numbers */
	printf("Nx = %d\n", Nx);
	x = (double *) malloc(Nx * sizeof(double));
	for(i=0; i<Nx; i++)
	{
		x[i] = i;
	}
	
	//double test[]={1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0};
	
	/* reading data from a file */

	
  /* reading WriteObserverSMesh.dat data  */
	
  double real_temp;
  double imag_temp;
  double *R = NULL;
  double *Im = NULL;
  int index, lineNum = 0, j = 0, ii = 1;
	
	FILE *file;
	//file=fopen("pXOYM150000.txt","r");
	//file=fopen("pF_NewTest.txt","r");
	file=fopen("pF.txt","r");
  //storing the data into the structure
  while(fscanf(file,"%lf %lf",&real_temp, &imag_temp)!=EOF)
    {
      if(R == NULL && Im == NULL)
	{

	  R = malloc(sizeof(real_temp));
	  Im = malloc(sizeof(imag_temp));
	  *R = real_temp;
	  *Im = imag_temp;

	}
      else
	{
	  ii++;
	  R = realloc(R,sizeof(R)*ii);
	  Im = realloc(Im,sizeof(Im)*ii);
	  index = ii-1;
	  *(R+index)=real_temp;
	  *(Im+index)=imag_temp;
	}
      
    }
  real_temp = 0.0;
  imag_temp = 0.0;
  while(index>=0)
    {
      printf("[%d]: R = %lf, Im = %lf\n",j, R[j], Im[j]);
      index--;
      lineNum = 1+j;
      j++;
      real_temp++;
      imag_temp++;
   
    }
  printf("----- checking ----> lineNum = %d\n",lineNum);
  printf("----- checking ----> [53]: point[53] = %lf	+%lfi\n", R[53], Im[53]);
  printf("----- checking ----> [4095]: Point[149] = %lf		+%lfi\n", R[4095], Im[4095]);
	
	free(R);free(Im);
  fclose(file); 
	

	
	
	
	
	
	Nx=lineNum;
	//Nx=8;	
		/* calculate NFFT as the next higher power of 2 >= Nx */
	NFFT = (int)pow(2.0, ceil(log((double)Nx)/log(2.0)));
	printf("NFFT = %d\n", NFFT);

	/* allocate memory for NFFT complex numbers (note the +1) */
	X = (double *) malloc((2*NFFT+1) * sizeof(double));

	/* Storing x(n) in a complex array to make it work with four1. 
	This is needed even though x(n) is purely real in this case. */
	//double R[]={4,1,0,1,0,1,0,1};
	//double Im[]={0,2.41,0,0.41,0,-0.41,0,-2.41};
	
	
	for(i=0; i<Nx; i++)
	{
		//X[2*i+1] = x[i];
		X[2*i+1] = R[i]; // real part
		//X[2*i+2] = 0.0;
		X[2*i+2] = -Im[i];// imaginary part
	}
	/* pad the remainder of the array with zeros (0 + 0 j) */
	for(i=Nx; i<NFFT; i++)
	{
		X[2*i+1] = 0.0;
		X[2*i+2] = 0.0;
	}

	printf("\nInput complex sequence (padded to next highest power of 2):\n");
	for(i=0; i<NFFT; i++)
	{
		printf("x[%d] = (%lf + j %lf)\n", i, X[2*i+1], X[2*i+2]);
	}

	/* calculate FFT *
	four1(X, NFFT, 1);

	printf("\nFFT:\n");
	for(i=0; i<NFFT; i++)
	{
		printf("X[%d] = (%lf + j %lf)\n", i, X[2*i+1], X[2*i+2]);
	}*/

	/*calculate IFFT */
	four1(X, NFFT, -1);

	/* normalize the IFFT */
	for(i=0; i<NFFT; i++)
	{
		X[2*i+1] /= NFFT;
		X[2*i+2] /= NFFT;
	}

	printf("\nComplex sequence reconstructed by IFFT:\n");
	FILE *fwriteIFFT;
	fwriteIFFT = fopen("IFFT.txt","w");
	for(i=0; i<NFFT; i++)
	{
		printf("x[%d] = (%lf + j %lf)\n", i, X[2*i+1], X[2*i+2]);
		fprintf(fwriteIFFT,"%12.10f		%12.10f\n", X[2*i+1], X[2*i+2]);
	}
	fclose(fwriteIFFT);

	return 0;
	
  
}


