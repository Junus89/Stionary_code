#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */


/* ----- including -----header defined functions ------*/
#include "FDSolver.h"
#include "FDSolver.c"



int main()
{
	
	/*reading ConstantVariable file including file names */
  FILE *fp;
  char fdName[128], sdName[128], odName[128];
  int ObsrSthetaNum, ObsrSFaiNum;
	//double MaX, MaY, MaZ;
  fp = fopen("ConstantV.txt","r");
  while(1 == fscanf(fp,"%s%*[^\n] %s%*[^\n] %s%*[^\n] %d%*[^\n] %lf%*[^\n] %d%*[^\n] %d%*[^\n] %lf%*[^\n]\
		 %lf%*[^\n] %lf%*[^\n] %d%*[^\n] %d%*[^\n]\n %lf%*[^\n] %lf%*[^\n] %lf%*[^\n] %lf%*[^\n]",fdName, sdName, odName,\
    	 &BNum,&R,&TNum,&FNum,&MaX,&MaY,&MaZ,&ObsrSthetaNum,&ObsrSFaiNum,&Rou_0,&C_0,&RMaTip,&Omega)){
    }
  fclose(fp);
  printf("Mx = %f My = %f, Mz = %f\n",MaX,MaY,MaZ);
  double HalfNum = (ObsrSthetaNum-1)*(ObsrSFaiNum+1)/2;
  printf("\n-------->Checking:\n ObsrSthetaNum = %d\n ObsrSFeinum = %d\n HalfNum = %lf\n",ObsrSthetaNum, ObsrSFaiNum, HalfNum);
  Gamma = sqrt(1/(1-(pow(MaX,2)+pow(MaY,2)+pow(MaZ,2)))); 
  printf("Gamma = %lf\n",Gamma);
	
	double T; T = 2*PI/Omega;
	
	double A; A = 1.0; /* source strength */
  DT = T/(TNum-1);          //...discrete time or sampling time? or sampling time step?


  //double Time[TNum];
  for(int i=0;i<TNum;i++)
  {
	  Time[i]=DT*i;
  }

 
  /* reading WriteObserverSMesh.dat data  */
  double XO_temp;
  double YO_temp;
  double ZO_temp;
  double *XO = NULL;
  double *YO = NULL;
  double *ZO = NULL;
  int OSNum = 0;
  int indexC, jC = 0, iC = 1;
  FILE *fp_ObMesh;
  fp_ObMesh = fopen("ObserverGeo.txt","r");
  //storing the data into the structure
  while(fscanf(fp_ObMesh,"%lf %lf %lf",&XO_temp, &YO_temp, &ZO_temp)!=EOF)
    {
      if(XO == NULL && YO == NULL && ZO == NULL)
	{

	  XO = malloc(sizeof(XO_temp));
	  YO = malloc(sizeof(YO_temp));
	  ZO = malloc(sizeof(ZO_temp));
	  *XO = XO_temp;
	  *YO = YO_temp;
	  *ZO = ZO_temp;
	}/* in case the memory is not enough, it hase to be reallocated */
      else
	{
	  iC++;
	  XO = realloc(XO,sizeof(XO)*iC);
	  YO = realloc(YO,sizeof(YO)*iC);
	  ZO = realloc(ZO,sizeof(ZO)*iC);
	  indexC = iC-1;
	  *(XO+indexC)=XO_temp;
	  *(YO+indexC)=YO_temp;
	  *(ZO+indexC)=ZO_temp;
	}
      /* showing and checking the result */
    }
  XO_temp = 0.0;
  YO_temp = 0.0;
  ZO_temp = 0.0;
  while(indexC>=0)
    {
     // printf("[%d]: XO = %lf, YO = %lf, ZO = %lf\n",jC, XO[jC], YO[jC], ZO[jC]);
      indexC--;
      OSNum = 1+jC;
      jC++;
      XO_temp++;
      YO_temp++;
      ZO_temp++;
    }
  printf("----- checking ----> OSNum = %d\n",OSNum);
  printf("----- checking ----> [53]: XO[0] = %lf, YO[0] = %lf, ZO[0] = %lf\n", XO[0], YO[0], ZO[0]);
  
  fclose(fp_ObMesh);
  /* freeing the allocated memory for 'WritesObserverSMeshData.dat */
 
	double *r,*RStar,*R,*VPFCoefX1,*VPFCoefY1,*VPFCoefZ1; 
	double complex *VPFCoefX2,*VPFCoefY2,*VPFCoefZ2;
	r = dvector(1,OSNum); RStar = dvector(1,OSNum); r= dvector(1,OSNum); R = dvector(1,OSNum);
	VPFCoefX1 = dvector(1,OSNum);VPFCoefX2 = dvector(1,OSNum);VPFCoefY1 = dvector(1,OSNum);
	VPFCoefY2 = dvector(1,OSNum);VPFCoefZ1 = dvector(1,OSNum);VPFCoefZ2 = dvector(1,OSNum);
	
	for(int i=0;i<OSNum;i++)
	{
		r[i] = sqrt(pow(XO[i],2)+pow(YO[i],2)+pow(ZO[i],2));
		RStar[i] = sqrt(pow(r[i],2)+pow(Gamma,2)*pow((MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i]),2))/Gamma;
		R[i] = pow(Gamma,2)*(RStar[i]-(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i]));
		VPFCoefX1[i] = (XO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaX)/(pow(Gamma,2)*RStar[i]);
		VPFCoefY1[i] = (YO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaY)/(pow(Gamma,2)*RStar[i]);
		VPFCoefZ1[i] = (ZO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaZ)/(pow(Gamma,2)*RStar[i]);

		VPFCoefX2[i] = (z1*Omega/C_0)*(XO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaX-pow(Gamma,2)*RStar[i]*MaX);
		VPFCoefY2[i] = (z1*Omega/C_0)*(YO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaY-pow(Gamma,2)*RStar[i]*MaY);
		VPFCoefZ2[i] = (z1*Omega/C_0)*(ZO[i]+pow(Gamma,2)*(MaX*XO[i]+MaY*YO[i]+MaZ*ZO[i])*MaZ-pow(Gamma,2)*RStar[i]*MaZ);
	}
	free_dvector(r,1,OSNum);free_dvector(RStar,1,OSNum);free_dvector(R,1,OSNum);free_dvector(VPFCoefX1,1,OSNum);
	free_dvector(VPFCoefX2,1,OSNum);free_dvector(VPFCoefY1,1,OSNum);free_dvector(VPFCoefY2,1,OSNum);
	free_dvector(VPFCoefZ1,1,OSNum);free_dvector(VPFCoefZ2,1,OSNum);
	
	double complex **VPF,**VPFCoef,**VPFX,**VPFY,**VPFZ,**p,**rou;
	double **realu,**realv,**realw,**realp,**realrou,**pXOY;
	VPF = dmatrix(1,OSNum,1,TNum);
	VPFCoef=dmatrix(1,OSNum,1,TNum);
	VPFX=dmatrix(1,OSNum,1,TNum);
	VPFY=dmatrix(1,OSNum,1,TNum);
	VPFZ=dmatrix(1,OSNum,1,TNum);
	realu = dmatrix(1,OSNum,1,TNum);
	realv = dmatrix(1,OSNum,1,TNum);
	realw = dmatrix(1,OSNum,1,TNum);
	realp = dmatrix(1,OSNum,1,TNum);
	realrou = dmatrix(1,OSNum,1,TNum);
	realp = dmatrix(1,OSNum,1,TNum);
	realrou = dmatrix(1,OSNum,1,TNum);
	pXOY = dmatrix(1,OSNum,1,TNum);
	
	for(int i=0;i<OSNum;i++)
	{
		for(int j=0;j<TNum;j++)
		{
	    VPF[i][j] = A*cexp(z1*Omega*(Time[j]-R[i]/C_0))/(4*PI*RStar[i]);
	    VPFCoef[i][j] = -A*cexp(z1*Omega*(Time[j]-R[i]/C_0))/(4*PI*pow(RStar[i],2));
    
	    VPFX[i][j] = VPFCoef[i][j]*(VPFCoefX1[i]+VPFCoefX2[i]);
	    VPFY[i][j] = VPFCoef[i][j]*(VPFCoefY1[i]+VPFCoefY2[i]);
	    VPFZ[i][j] = VPFCoef[i][j]*(VPFCoefZ1[i]+VPFCoefZ2[i]);
			
			realu[i][j] = real(VPFX[i][j]);
			realv[i][j] = real(VPFY[i][j]);
			realw[i][j] = real(VPFZ[i][j]);
			p[i][j] = -Rou_0*(z1*Omega*VPF[i][j]+C_0*(MaX*VPFX[i][j]+MaY*VPFY[i][j]+MaZ*VPFZ[i][j]));
			rou[i][j] = p[i][j]/pow(C_0,2);
			realp[i][j] = real(p[i][j]);
			realrou[i][j] = real(rou[i][j]);
			pXOY[i][j] = realp[i][j];
			printf("pXOY[%d][%d] = %lf\n",i,j,pXOY[i][j]);
			
		}
	}
	
	printf("------ Pressure Prediction Is Done!-------");
	
	double *RMSpXOY,*ObserverSTheta=NULL,anorm;
	RMSpXOY = dvector(1,OSNum+1);

	for(int j=0;j<OSNum;j++)
	{
		for(int k=0;k<TNum;k++)
		{
			anorm = 0.0;
			//RMSpXOY[j] = norm(pXOY[j][k])/sqrt(TNum[k]);
			anorm += fabs(pXOY[j][k]);
			RMSpXOY[j] = sqrt(anorm)/sqrt(TNum);
		}
	}
	for(int i=0;i<OSNum+1;i+=2*PI)
	{
		ObserverSTheta[i] = i;
	}	
	RMSpXOY[OSNum+1] = RMSpXOY[1];
	
	
	
	
	
	
	free_dmatrix(VPF,1,OSNum,1,TNum);free_dmatrix(VPFCoef,1,OSNum,1,TNum);free_dmatrix(VPFX,1,OSNum,1,TNum);
	free_dmatrix(VPFY,1,OSNum,1,TNum);free_dmatrix(VPFZ,1,OSNum,1,TNum);free_dmatrix(realu,1,OSNum,1,TNum);
	free_dmatrix(realv,1,OSNum,1,TNum);free_dmatrix(realw,1,OSNum,1,TNum);free_dmatrix(p,1,OSNum,1,TNum);
	free_dmatrix(rou,1,OSNum,1,TNum);free_dmatrix(realp,1,OSNum,1,TNum);free_dmatrix(realrou,1,OSNum,1,TNum);
	free_dmatrix(pXOY,1,OSNum,1,TNum); free_dvector(RMSpXOY,1,OSNum+1);free_dvector(ObserverSTheta,1,OSNum+1);


  return 0;
}


