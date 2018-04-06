#include<stdio.h>
#include<math.h>

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
	int i;
	double **outputAr;
	outputAr = (double **)malloc((nrh+1)*sizeof(double *));
	if(outputAr == NULL)
	{
		printf("No enough memory space!");
		exit(1);
	}
	for(i=0;i<nrh+1;i++)
	{
		outputAr[i]=(double *)malloc((nch+1)*sizeof(double));
		if(outputAr[i]==NULL)
		{
			printf("No enough memory space!\n");
			exit(1);
		}
	}
	return outputAr;
}

void free_dmatrix(double **inputAr,long nrl,long nrh, long ncl, long nch)
{
	int i;
	for(i=0;i<=nrh;i++)
		free(inputAr[i]);
	free(inputAr);
}


double anorm2(double **a, int n,int m)
{
	int i,j;
	double anorm;

	for (i=1;i<=n;i++)
		for (j=1;j<=m;j++)
			anorm=0.0;
			anorm += fabs(a[i][j]);
	return sqrt(anorm);
}


int main()
{
	int i,j,n=40,m=30;
	double **a,anorm,aanorm,Anorm;
	a = dmatrix(1,n,1,n);
		
	for(i=1;i<=n;i++){
		for(j=1;j<=m;j++){
			anorm=0.0;Anorm=0.0;
			a[i][j] = 1.0/(1+i+j);
			
			printf("a[%d][%d] = %lf\n",i,j,a[i][j]);
			//printf("b[%d][%d] = %lf\n",i,j,b[i][j]);
			anorm += fabs(a[i][j]); 
			aanorm = sqrt(anorm);
			anorm2(a,n,m);
			//Anorm = anorm2(a,n,m);
		}
	}
	printf("the norm of the matrix a is : %lf\n",anorm);
	//Anorm = anorm2(a,n,m);
	printf("the anorm2 of the matrix a is: %lf\n",Anorm);
	//printf("the anorm2 of the matrix b is: %lf\n",bnorm);
	free_dmatrix(a,1,n,1,n);
	
	return 0;
	
}
	