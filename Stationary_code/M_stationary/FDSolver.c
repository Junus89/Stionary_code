#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */

#include"FDSolver.h"

/* ---------- function definitions --------- */

/* Memory allocation function for different dimension arrays */
//***************************1D*******************************************/
float *vector(long nl,long nh)
/*	
 */
{
    float *outputAr; 

    outputAr = (float *)malloc((nh+1)*sizeof(float));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
    return outputAr;
}

void free_vector(float *inputAr,long nl,long nh)
{
    free(inputAr);
    return;
}

/*----------------------------------------------------------- */
double *dvector(long nl,long nh)
/*	
 */
{
    double *outputAr; 

    outputAr = (double *)malloc((nh+1)*sizeof(double));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
    return outputAr;
}


void free_dvector(double *inputAr,long nl,long nh)
{
    free(inputAr);
    return;
}

//***************************2D*******************************************/
float **matrix(long nrl,long nrh, long ncl,long nch)
/*	
 */
{
    int i;
    float **outputAr;
    outputAr = (float **)malloc((nrh+1)*sizeof(float *));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (float *)malloc((nch+1)*sizeof(float));
        if (outputAr[i] == NULL)
		{
	       printf("No enough memory space.");
	       exit(1);
		}
    }		
    return outputAr;
}

void free_matrix(float **inputAr,long nrl,long nrh,long ncl,long nch)
/*	
 */
{
    int i;
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


/*-----------------------------------------------------------------------*/

double **dmatrix(long nrl,long nrh, long ncl,long nch)
/*	
 */
{
    int i;
    double **outputAr;
    outputAr = (double **)malloc((nrh+1)*sizeof(double *));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (double *)malloc((nch+1)*sizeof(double));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}		
    return outputAr;
}

void free_dmatrix(double **inputAr,long nrl,long nrh,long ncl,long nch)
/*	
 */
{
    int i;
    for (i = 0; i <= nrh; i++)  free(inputAr[i]);
    free(inputAr);		
    return;
}

//***************************3D*******************************************/

float ***f3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh)
/*
 */
{
    int i,j;
    float ***outputAr;
    outputAr = (float ***)malloc((nrh+1)*sizeof(float**));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (float **)malloc((nch+1)*sizeof(float*));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}

    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{

	        outputAr[i][j] = (float*)malloc((ndh+1)*sizeof(float));
	        if (outputAr[i][j] == NULL)
			{
	            printf("No enough memory space.");
	            exit(1);
			}
		}
    }	
    return outputAr;
}



void free_f3dmatrix(float ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
/*	
 */
{
    int i,j;
	
    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{
	        free(inputAr[i][j]);
		}
	}
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


/*-----------------------------------------------------------------------*/

double ***d3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh)
/*
 */
{
    int i,j;
    double ***outputAr;
    outputAr = (double ***)malloc((nrh+1)*sizeof(double **));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (double **)malloc((nch+1)*sizeof(double*));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}

    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{

	        outputAr[i][j] = (double *)malloc((ndh+1)*sizeof(double));
	        if (outputAr[i][j] == NULL)
			{
	            printf("No enough memory space.");
	            exit(1);
			}
		}
	}	
    return outputAr;
}

void free_d3dmatrix(double ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
/*	
 */
{
    int i,j;
	
    for (i = 0; i < nrh+1; i++)	
    {
        for (j = 0; j < nch+1; j++)
		{
	        free(inputAr[i][j]);
		}
	}
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


//***************************4D*******************************************/


float ****f4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh)
{
	int i,j, k;
	float ****outputAr;
	outputAr = (float ****)malloc((nrh+1)*sizeof(float ***));
	if(outputAr == NULL)
	{
		printf("No enough memory space !");
		exit(1);
	}
	for(i=0;i< nrh+1;i++)
	{
		outputAr[i] = (float ***)malloc((nch+1)*sizeof(float **));
		if(outputAr == NULL)
		{
			printf("No enough memory space!");
			exit(1);
		}
	}
	for (i=0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			outputAr[i][j] = (float **)malloc((ndh+1)*sizeof(float*));
			if(outputAr[i][j] == NULL)
			{
				printf("No enough memory space!");
				exit(1);
			}
		}
	}
	
	for(i = 0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			for(k=0;k<ndh+1;k++)
			{
				outputAr[i][j][k] = (float *)malloc((nfh+1)*sizeof(float));
				if(outputAr[i][j][k]==NULL)
				{
					printf("No enough memory space!");
					exit(1);
				}
			}
		}	
	}
	
	return outputAr;
}

void free_f4dmatrix(float ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh)
{
	int i,j,k;

	for(i=0;i<nrh+1;i++)
		{
		for(j=0;j<nch+1;j++)
			{
				for(k=0;k<ndh+1;k++)
					{
						free(inputAr[i][j][k]);
					}
					free(inputAr[i][j]);
			}
			free(inputAr[i]);
		}
/*
	for(i=0;i<nrh+1;i++)
		{
			for(j=0;j<nch+1;j++)
				{
					free(inputAr[i][j]);
				}
		}

	for(i=0;i<nrh+1;i++) free(inputAr[i]);*/
	free(inputAr);

}

/*-----------------------------------------------------------------------*/

double ****d4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh)
{
	int i,j, k;
	double ****outputAr;
	outputAr = (double ****)malloc((nrh+1)*sizeof(double ***));
	if(outputAr == NULL)
	{
		printf("No enough memory space !");
		exit(1);
	}
	for(i=0;i< nrh+1;i++)
	{
		outputAr[i] = (double ***)malloc((nch+1)*sizeof(double **));
		if(outputAr == NULL)
		{
			printf("No enough memory space!");
			exit(1);
		}
	}
	for (i=0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			outputAr[i][j] = (double **)malloc((ndh+1)*sizeof(double*));
			if(outputAr[i][j] == NULL)
			{
				printf("No enough memory space!");
				exit(1);
			}
		}
	}
	
	for(i = 0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			for(k=0;k<ndh+1;k++)
			{
				outputAr[i][j][k] = (double *)malloc((nfh+1)*sizeof(double));
				if(outputAr[i][j][k]==NULL)
				{
					printf("No enough memory space!");
					exit(1);
				}
			}
		}	
	}
	
	return outputAr;
}

void free_d4dmatrix(double ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh)
{
	int i,j,k;

	for(i=0;i<nrh+1;i++)
		{
		for(j=0;j<nch+1;j++)
			{
				for(k=0;k<ndh+1;k++)
					{
						free(inputAr[i][j][k]);
					}
					free(inputAr[i][j]);
			}
			free(inputAr[i]);
		}
/*
	for(i=0;i<nrh+1;i++)
		{
			for(j=0;j<nch+1;j++)
				{
					free(inputAr[i][j]);
				}
		}

	for(i=0;i<nrh+1;i++) free(inputAr[i]);*/
	free(inputAr);

}



double *Ones(int Tnum)
{
  double *One;
  One = dvector(1,Tnum);
  for(int i=0;i<Tnum;i++)
  {
	  One[i] = 1.0;
  }
  return One;
}




