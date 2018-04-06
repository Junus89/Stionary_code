#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */


/* ----- including -----header defined functions ------*/
#include "xmalloc.h"
#include "xmalloc.c"
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
	double T; T = 2*PI/Omega;
	

  

  /* Mathematical relations of TNum data */
  //double Tint = 27*TR; /* Noninteger ratio */
  double Tint = TR; /* Integer ratio */
	//double Tint =1.0;
  DT = Tint/(TNum-1);          //...discrete time or sampling time? or sampling time step?


  //double Time[TNum];
  for(int i=0;i<TNum;i++)
  {
	  Time[i]=i/(TNum*1.0-1.0)*Tint;
	  //printf("Time[%d] = %4.6f\n",i,(double)Time[i]);
  }
    
  
  int NFFT = (int)pow(2.0,ceil(log((double)FNum)/log(2.0))); 
  double ODT =Tint/NFFT;
  printf("the next power of FNum->NFFT = %d\n",NFFT);
  printf("ODT = %lf\n",ODT);
  /* for construction of OTime */
  //double *OTime;r
  
  OTime = make_vector(OTime,NFFT-1); /* calling forOTime function and initializing it */
  for(int i=0;i<NFFT-1;i++)
    {
		OTime[i] = ODT*i;
      printf("OTime[%d] = %4.4f\n",i,OTime[i]);
    }
  printf("----checking ----> OTime[5] = %4.4f\n",OTime[5]); /*checking for OTime */
  //double DF = ((1/ODT)/2)*(1/(1.0*FNum/2)); /* here as FNum is int type, should multiply with 1.0 to get the double type DF result */
  //double DF = 1/Tint;
  double DF=1.0;
  printf("----checking ----> DF = %4.9f\n",DF);
  
  /* reading bladesource.dat data */

  double Dx_temp;
  double Dy_temp;
  double Dz_temp;
  double Dp_temp;
  double Drou_temp;
  double Du_temp;
  double Dv_temp;
  double Dw_temp;
	
  double *DSX = NULL; double *DSp = NULL; double *DSv = NULL;
  double *DSY = NULL; double *DSrou = NULL; double *DSw = NULL;
  double *DSZ = NULL; double *DSu = NULL;
  int index, DSNum = 0, j = 0, i = 1;

  FILE *fp_BladeS;
  fp_BladeS = fopen("bladesource.txt","r");
  //storing the data into the structure
  while(fscanf(fp_BladeS,"%lf %lf %lf %lf %lf %lf %lf %lf",&Dx_temp, &Dy_temp, &Dz_temp,&Dp_temp, &Drou_temp, &Du_temp,&Dv_temp, &Dw_temp)!=EOF)
    {
      if(DSX == NULL && DSY == NULL && DSX == NULL && DSp == NULL && DSrou == NULL && DSu == NULL && DSv == NULL && DSw == NULL)
	{

	  DSX = malloc(sizeof(Dx_temp)); DSp = malloc(sizeof(Dp_temp)); 		DSv = malloc(sizeof(Dv_temp));
	  DSY = malloc(sizeof(Dy_temp)); DSrou = malloc(sizeof(Drou_temp)); DSw = malloc(sizeof(Dw_temp));
	  DSZ = malloc(sizeof(Dz_temp)); DSu = malloc(sizeof(Du_temp));
	  *DSX = Dx_temp; *DSp = Dp_temp; 		*DSv = Dv_temp;
	  *DSY = Dy_temp; *DSrou = Drou_temp; *DSw = Dw_temp;
	  *DSZ = Dz_temp; *DSu = Du_temp;
	}/* in case the memory is not enough, it hase to be reallocated */
      else
	{
	  i++;
	  DSX = realloc(DSX,sizeof(DSX)*i); DSp = realloc(DSp,sizeof(DSp)*i); 			DSv = realloc(DSv,sizeof(DSv)*i);
	  DSY = realloc(DSY,sizeof(DSY)*i); DSrou = realloc(DSrou,sizeof(DSrou)*i); DSw = realloc(DSw,sizeof(DSw)*i);
	  DSZ = realloc(DSZ,sizeof(DSZ)*i); DSu = realloc(DSu,sizeof(DSu)*i);
	  index = i-1;
	  *(DSX+index)=Dx_temp; *(DSp+index)=Dp_temp; 		*(DSv+index)=Dv_temp;
	  *(DSY+index)=Dy_temp; *(DSrou+index)=Drou_temp; *(DSw+index)=Dw_temp;
	  *(DSZ+index)=Dz_temp; *(DSu+index)=Du_temp;
	}
      /* showing and checking the result */
    }
  Dx_temp = 0.0; Dp_temp = 0.0; 	Dv_temp = 0.0;
  Dy_temp = 0.0; Drou_temp = 0.0; Dw_temp = 0.0;
  Dz_temp = 0.0; Du_temp = 0.0;
  while(index>=0)
    {
      printf("[%d]: DSX= %lf, DSY = %lf, DSZ = %lf, DSp= %lf, DSrou = %lf, DSu = %lf, DSv= %lf, DSw = %lf\n",j, DSX[j], DSY[j], DSZ[j],\
				DSp[j], DSrou[j], DSu[j], DSv[j], DSw[j]);
      index--;
      DSNum = 1+j;
      j++;
      Dx_temp++; Dp_temp++; 	Dv_temp++;
      Dy_temp++; Drou_temp++; Dw_temp++;
      Dz_temp++; Du_temp++;
    }
  printf("----- checking ----> DSNum = %d\n",DSNum);
  /*printf("----- checking ----> [53]: XO[53] = %lf, YO[53] = %lf, ZO[53] = %lf\n", XO[53], YO[53], ZO[53]);
  printf("----- checking ----> [99]: XO[99] = %lf, YO[99] = %lf, ZO[99] = %lf\n", XO[99], YO[99], ZO[99]);
  printf("----- checking ----> [149]: XO[149] = %lf, YO[1499] = %lf, ZO[149] = %lf\n", XO[149], YO[149], ZO[149]);*/
  fclose(fp_BladeS);
	
	



  /* reading bladenormals.dat data  */
  double nx_temp;
  double ny_temp;
  double nz_temp;
  double *Nx = NULL;
  double *Ny = NULL;
  double *Nz = NULL;
  int normalLNum=0;
	int indexN, jN = 0, iN = 1;
  FILE *fp_normals;
  fp_normals = fopen("bladenormals.txt","r");
  //storing the data into the structure
  while(fscanf(fp_normals,"%lf %lf %lf",&nx_temp, &ny_temp, &nz_temp)!=EOF)
    {
      if(Nx == NULL && Ny == NULL && Nz == NULL)
	{

	  Nx = malloc(sizeof(nx_temp));
	  Ny = malloc(sizeof(ny_temp));
	  Nz = malloc(sizeof(nz_temp));
	  *Nx = nx_temp;
	  *Ny = ny_temp;
	  *Nz = nz_temp;
	}/* in case the memory is not enough, it hase to be reallocated */
      else
	{
	  iN++;
	  Nx = realloc(Nx,sizeof(Nx)*iN);
	  Ny = realloc(Ny,sizeof(Ny)*iN);
	  Nz = realloc(Nz,sizeof(Nz)*iN);
	  indexN = iN-1;
	  *(Nx+indexN)=nx_temp;
	  *(Ny+indexN)=ny_temp;
	  *(Nz+indexN)=nz_temp;
	}
      /* showing and checking the result */
    }
  nx_temp = 0.0;
  ny_temp = 0.0;
  nz_temp = 0.0;
  while(indexN>=0)
    {
      printf("[%d]: Nx = %e, Ny = %e, Nz = %e\n",jN, Nx[jN], Ny[jN], Nz[jN]);
      indexN--;
      normalLNum = 1+jN;
      jN++;
      nx_temp++;
      ny_temp++;
      nz_temp++;
    }
  printf("----- checking ----> normalLLNum = %d\n",normalLNum);
  printf("----- checking ----> [53]: Nx[53] = %lf, Ny[53] = %lf, Nz[53] = %lf\n", Nx[53], Ny[53], Nz[53]);
  printf("----- checking ----> [990]: Nx[990] = %lf, YO[990] = %lf, Nz[990] = %lf\n", Nx[990], Ny[990], Nz[990]);
  printf("----- checking ----> [149]: Nz[10049] = %lf, Nz[10049] = %lf, Nz[10049] = %lf\n", Nx[10049], Ny[10049], Nz[10049]);
  fclose(fp_normals);






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
  fp_ObMesh = fopen("WriteObserverSMesh.txt","r");
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
  


  Gamma = sqrt(1/(1-(pow(MaX,2)+pow(MaY,2)+pow(MaZ,2)))); 
  printf("Gamma = %lf\n",Gamma);

 
  double *Area, **UnitNormal;
  Area = make_vector(Area,DSNum);
  UnitNormal = make_dmatrix(DSNum,3);
  for(int m=0;m<DSNum;m++)
    {
      Area[m] = sqrt(pow(Nx[m],2)+pow(Ny[m],2)+pow(Nz[m],2));
      printf("Area[%d]  = %lf\n",m,Area[m]);
      UnitNormal[m][0] = Nx[m]/Area[m];
      UnitNormal[m][1] = Ny[m]/Area[m];
      UnitNormal[m][2] = Nz[m]/Area[m];
      printf("UnitNormal[%d][0] = %g\n Area[%d] = %g\n",m,UnitNormal[m][0],m,Area[m]);
      printf("UnitNormal[%d][1] = %g\n Area[%d] = %g\n",m,UnitNormal[m][1],m,Area[m]); 
      printf("UnitNormal[%d][2] = %g\n Area[%d] = %g\n",m,UnitNormal[m][2],m,Area[m]);   
    }
  /*double *DataSArea, **DataSVector;
  r = Area;
  DataSVector = UnitNormal;*/
  
  double complex ****PBD1,****PBD2,****PBD3; 
	double complex ***PB1,***PB2,***PB3;
	double complex **P1,**P2,**P3;
	double complex **pF, **pXOYM;
	double  **Op, *OF;

  PBD1 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD2 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  PBD3 = make_4Ddmatrix(FNum,OSNum,BNum,DSNum);
  

  PB1 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB2 = make_3Ddmatrix(FNum,OSNum,BNum);
  PB3 = make_3Ddmatrix(FNum,OSNum,BNum);
  

  P1 = make_dmatrix(FNum,OSNum);
  P2 = make_dmatrix(FNum,OSNum);
  P3 = make_dmatrix(FNum,OSNum);
  
  


  DataXR=make_vector(DataXR,TNum);
  DataYR=make_vector(DataYR,TNum);
  DataZR=make_vector(DataZR,TNum);

  
  DOrX=make_vector(DOrX,TNum);
  DOrY=make_vector(DOrY,TNum);
  DOrZ=make_vector(DOrZ,TNum);
  DOr=make_vector(DOr,TNum);
  
  DORStar=make_vector(DORStar,TNum);
  DOR=make_vector(DOR,TNum);
  
  DORStarX=make_vector(DORStarX,TNum);
  DORStarY=make_vector(DORStarY,TNum);
  DORStarZ=make_vector(DORStarZ,TNum);
  
  DORX=make_vector(DORX,TNum);
  DORY=make_vector(DORY,TNum);
  DORZ=make_vector(DORZ,TNum);
  RGamma=make_vector(RGamma,TNum);
  
  Vx=make_vector(Vx,TNum);
  Vy=make_vector(Vy,TNum);
  Vz=make_vector(Vz,TNum);
  Q=make_vector(Q,TNum);
  
  Lx=make_vector(Lx,TNum);
  Ly=make_vector(Ly,TNum);
  Lz=make_vector(Lz,TNum);
  FxM=make_vector(FxM,TNum);
  FyM=make_vector(FyM,TNum); 
  FzM=make_vector(FzM,TNum); 
  FxP=make_vector(FxP,TNum); 
  FyP=make_vector(FyP,TNum); 
  FzP=make_vector(FzP,TNum);
  FRStarM=make_vector(FRStarM,TNum); 
  FRStarP=make_vector(FRStarP,TNum);
  FRM=make_vector(FRM,TNum); 
  FRP=make_vector(FRP,TNum);
  //double *One;
  pF = make_dmatrix(FNum,OSNum);
  
  //double **pXOYM,**Op,*OF;
  //double complex SP1; double complex SP2; double complex SP3;
	/*SP1=0.0 + 0.0*I;
	SP2=0.0 + 0.0*I;
	SP3=0.0 + 0.0*I;*/
  FILE *fwriteSP1,*fpPB1,*fpPB2;
  fwriteSP1 = fopen("PBD1.txt","w");	
	fpPB1 = fopen("PB1.txt","w");	
	fpPB2 = fopen("PB2.txt","w");	

	
  for(int n =0; n<FNum;n++)
  //for(int n =0; n< FNum;n=n+16)
  //for(int n =0; n< 2;n++)
  {

	  double Omega =0.0, ka=0.0;
	  Omega= OmegaM+BNum*n*OmegaR; /* first omega formula */
	  //Omega = (n)*2*PI*DF; /* second omega formula */
	  ka = Omega/C_0;
	  //printf("Omega[%d]=%4.4g   ka[%d] = %4.4g\n",n,Omega[n],n,ka[n]);
  

  
	  
	  for(int m=20;m<21;m++)//OSNum=
	  {
		  OX = XO[m];//the first value of XO vector
		  OY = YO[m];// the first value of YO vector
		  OZ = ZO[m];// the first value of ZO vector

		  /*
		  double complex Transa1 = 0.0 + 0.0*I;
		  double complex Transb1 = 0.0 + 0.0*I;
		  double complex Transc1 = 0.0 + 0.0*I;
*/
		  double complex Transa2 = 0.0 + 0.0*I;
		  double complex Transb2 = 0.0 + 0.0*I;
		  double complex Transc2 = 0.0 + 0.0*I;
	
			
			
		  for(int k=0;k<BNum;k++)
		  {
			  Theta = 2*((k+1)-1)*PI/BNum; // Theta = 2*(k-1)*PI/BNum;

			  double complex Transa1 = 0.0 + 0.0*I;
			  double complex Transb1 = 0.0 + 0.0*I;
			  double complex Transc1 = 0.0 + 0.0*I;
		  
		  		for(int j=0;j<DSNum;j++)
		  	  	{
		  				printf("Harmonic number %d, Observer number %d, Blade number %d, Source number %d\n",n,m+1,k+1,j+1); // since m, k, j starts from zero, add it with 1
					    	printf("------checking-----: k = %d\n",k);
							DX=DSX[j];
							DY=DSY[j];
							DZ=DSZ[j];
							DR=sqrt(pow(DX,2)+pow(DY,2));
              						DVX = UnitNormal[j][0];
							DVY = UnitNormal[j][1];
							DVZ = UnitNormal[j][2];
							A = Area[j];
							Dp = DSp[j];
							Drou = DSrou[j];
							Du = DSu[j];
							Dv = DSv[j];
							Dw = DSw[j];
					
						
					// calling fifthLoop function
							//fifthLoop(OmegaR, Omega, MaX, MaY, MaZ, ka, pF,SP1,SP2,SP3);
							fifthLoop(OmegaR, Omega, MaX, MaY, MaZ, ka);
						
					
							PBD1[n][m][k][j] = SP1;
							PBD2[n][m][k][j] = SP2;
							PBD3[n][m][k][j] = SP3;
							
							fprintf(fwriteSP1,"%12.10f	%12.10f\n",creal(PBD1[n][m][k][j]),cimag(PBD1[n][m][k][j]));
					
							Transa1 += PBD1[n][m][k][j];
							Transb1 += PBD2[n][m][k][j];
							Transc1 += PBD3[n][m][k][j];
				
							printf("\n------checking-----: SP1 = %g + %gi\n",creal(SP1),cimag(SP1));
							printf(" PBD1[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD1[n][m][k][j]),cimag(PBD1[n][m][k][j]));
          		printf(" PBD2[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD2[n][m][k][j]),cimag(PBD2[n][m][k][j]));
          		printf(" PBD3[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD3[n][m][k][j]),cimag(PBD3[n][m][k][j]));
							/*PB1[n][m][k] = Transa1;
							PB2[n][m][k] = Transb1;
							PB3[n][m][k] = Transc1;*/

					
		  	  	}
					//	printf("the output is PBD1[%d][%d][%d][%d] = %g + %gi\n",n,m,k,j,creal(PBD1[n][m][k][j]),cimag(PBD1[n][m][k][j]));
					    printf(" Transa1 = %g + %gi\n",creal(Transa1),cimag(Transa1));
							PB1[n][m][k] = Transa1;
							PB2[n][m][k] = Transb1;
							PB3[n][m][k] = Transc1;
							

							
						 	Transa2 += PB1[n][m][k];
							Transb2 += PB2[n][m][k];
							Transc2 += PB3[n][m][k]; 
						 
							
							/*
						  Transa2 = PB1[n][m][k];
							Transb2 = PB2[n][m][k];
							Transc2 = PB3[n][m][k]; */
						
						
						printf("PB1[%d][%d][%d] = %g + %g\n",n,m,k,creal(PB1[n][m][k]),cimag(PB1[n][m][k]));
        		printf("PB2[%d][%d][%d] = %g + %g\n",n,m,k,creal(PB2[n][m][k]),cimag(PB2[n][m][k]));
        		printf("PB3[%d][%d][%d] = %g + %g\n",n,m,k,creal(PB3[n][m][k]),cimag(PB3[n][m][k]));
					
				    fprintf(fpPB1,"%12.10f	%12.10f\n",creal(PB1[n][m][k]),cimag(PB1[n][m][k]));
						fprintf(fpPB2,"%12.10f	%12.10f\n",creal(PB2[n][m][k]),cimag(PB2[n][m][k]));
			}
			//printf(" Transa2 = %g + %gi\n",creal(Transa2),cimag(Transa2));
			
		  P1[n][m] = Transa2;
			P2[n][m] = Transb2;
			P3[n][m] = Transc2;
			
			
			
			printf("P1[%d][%d] = %g + %g\n",n,m,creal(P1[n][m]),cimag(P1[n][m])); // m+1 just for show its index not zero but 1
			printf("P2[%d][%d] = %g + %g\n",n,m,creal(P2[n][m]),cimag(P2[n][m]));
			printf("P3[%d][%d] = %g + %g\n",n,m,creal(P3[n][m]),cimag(P3[n][m]));
			
			// pressure final
		     	// OSNum = 1;
		  pF[n][m] = (P1[n][m]+P2[n][m]+P3[n][m])/(4*PI*Tint);
		 // printf("pF[%d][%d]= %g + %g\n",n,m,creal(pF[n][m]),cimag(pF[n][m]));// m+1 just for show its index not zero but 1
			
			
	  }
	  

  }
	fclose(fwriteSP1);fclose(fpPB1);fclose(fpPB2);
	
	/*
	for(int n=2;n<FNum;n++){
		for(int m = 20;m<21;m++){
			pF[n][m] -=pF[1][m] ;
		}
	}*/
	
	
  /* Writing the Pressure spectrum into a Data */
  FILE *fwritepF;
  fwritepF = fopen("pF.txt","w");	
	  
	for(int j=0;j<FNum;j++)
	  {
			//for(int i=0;i<OSNum;i++)
            for(int i=20;i<21;i++)
			{
									
				fprintf(fwritepF,"%12.10f	%12.10f\n",creal(pF[j][i]),cimag(pF[j][i]));
			}
	  }
  fclose(fwritepF); 
	
	FILE *fRgama;
	fRgama = fopen("RGamma.txt","w");
	
	for(int t = 0; t<TNum;t++)
	{
		fprintf(fRgama,"%12.10f\n",RGamma[t]);
	}
	fclose(fRgama);
	
	/*
  FILE *fwriteSP1;
  fwritepF = fopen("PBD1.txt","w");	
	  
	for(int j=0;j<FNum;j++)
	  {
			//for(int i=0;i<OSNum;i++)
            for(int i=20;i<21;i++)
						{
							for(int s=0;s<BNum;s++)
							{
								for(int r=0;r<DSNum;r++)
								{
									fprintf(fwriteSP1,"%12.10f	%12.10f\n",creal(PBD1[j][i][s][r]),cimag(PBD1[j][i][s][r]));
								}
							}
						}
	  }
  fclose(fwriteSP1); */

	
	

 Op = make_dmatrix(FNum,OSNum);
 //pXOYM = make_dmatrix(FNum,OSNum);
  
		 
    for(int i=0;i<FNum;i++)
	{
            for(int j=20;j<21;j++)
		{
				printf("pF[%d][%d]= %g + %g\n",i,j,creal(pF[i][j]),cimag(pF[i][j]));
		}
	}
	
  /* Writing the Pressure spectrum into a Data */
 /* FILE *fwrite33;
  fwrite33 = fopen("pXOYM.txt","w");	
	  
	for(int j=0;j<FNum;j++)
	  {
			for(int i=20;i<21;i++)
			{
			  fprintf(fwrite33,"%12.9f	%12.9f\n",creal(pXOYM[j][i]),cimag(pXOYM[j][i]));
			}
	  }
  fclose(fwrite33);*/
	
  //double Opreal=0;
  //double Opimag=0;

  for(int i=0;i<FNum;i++)
  {
        //for(int i=0;i<OSNum;i++)
		for(int j=20;j<21;j++)
			{

				Op[i][j] = 2*sqrt(pow(creal(pF[i][j]),2)+pow(cimag(pF[i][j]),2));

	  			printf("Op[%d][%d] = %12.10f\n",i,j,Op[i][j]);
	  
			}

    }




  
  
  int rr;
  OF = make_vector(OF,FNum);
  for(rr=0;rr<FNum;rr++)
  {
	  //OF[rr] = (OmegaM+BNum*(rr-5)*OmegaR)/(2*PI);
	  OF[rr] = (OmegaM+BNum*rr*OmegaR)/(2*PI);
	  //OF[rr] = ((rr)*2*PI*DF)/(2*PI);
		
	  printf("OF[%d] = %g\n",rr,OF[rr]);
  }


  
  FILE *fwrite1;
  fwrite1 = fopen("FDPressureSpectrum.txt","w");
  //for(int k=0;k<OSNum;k++)
  for(int k=20;k<21;k++)
  {
	  
	  for(int j=0;j<FNum;j++)
	  {
		  fprintf(fwrite1,"%12.10f\t %12.10f\n",OF[j],Op[j][k]);
	  }
	  
  }
  fclose(fwrite1);
  
  
  
  
  
  
  
  
  
  
  
  
  //------checking variable values-------
  printf("----- checking ----: fR: %g\n",fR);
  printf("PBD1[0][0][0][0] = %g + %gi\n",creal(PBD1[0][0][0][0]),cimag(PBD1[0][0][0][0]));
  printf("----- checking ----: XO = %g, YO = %g, ZO = %g\n", OX, OY, OZ);
  printf("----- checking ----: Theta: %g\n",Theta);
  printf("------checking-----: DataXR[3] = %g\n",DataXR[3]);
  printf("------checking-----: DataYR[32] = %g\n",DataYR[32]);
  printf("------checking-----: DataZR[3] = %g\n",DataZR[34]);
  printf("----- checking ----: DX = %g, DY = %g, DZ = %g, DR = %g\n",DX,DY,DZ,DR);
  printf("------checking-----: DORStar[0] = %g\n",DORStar[0]);
	printf("------checking-----: DORStar[TNum-1] = %g\n",DORStar[TNum-1]);
  printf("------checking-----: DOrX[0] = %g\n",DOrX[0]);
  printf("------checking-----: DOrX[TNum-1] = %g\n",DOrX[TNum-1]);
  printf("------checking-----: DORY[0] = %g\n",DOrY[0]);
  printf("------checking-----: DOrY[TNum-1] = %g\n",DOrY[TNum-1]);
  printf("------checking-----: DOrZ[43] = %g\n",DOrZ[0]);
  printf("------checking-----: DORZ[43] = %g\n",DOrZ[TNum-1]);
  printf("------checking-----: DOr[43] = %g\n",DOr[143]);
  printf("------checking-----: DOR[43] = %g\n",DOR[143]);
  printf("------checking-----: DORStarX[0] = %g\n",DORStarX[0]);
	printf("------checking-----: DORStarX[TNum-1] = %g\n",DORStarX[TNum-1]);
  printf("------checking-----: DORStarY[0] = %g\n",DORStarY[0]);
	printf("------checking-----: DORStarY[TNum-1] = %g\n",DORStarX[TNum-1]);
  printf("------checking-----: DORStarZ[0] = %g\n",DORStarZ[0]);
	printf("------checking-----: DORStarZ[TNum-1] = %g\n",DORStarX[TNum-1]);
  printf("------checking-----: DORX[143] = %g\n",DORX[0]);
  printf("------checking-----: DORZ[343] = %g\n",DORX[TNum-1]);
  printf("------checking-----: DORY[0] = %g\n",DORY[0]);
  printf("------checking-----: DORY[TNum-1] = %g\n",DORY[TNum-1]);
  printf("------checking-----: DORZ[0] = %g\n",DORZ[0]);
	printf("------checking-----: DORZ[TNum-1] = %g\n",DORZ[TNum-1]);
  printf("------checking-----: RGamma[0] = %g\n",RGamma[0]);
  printf("------checking-----: RGamma[TNum-1] = %g\n",RGamma[TNum-1]);
  printf("------checking-----: Vx[0] = %g\n",Vx[0]);
  printf("------checking-----: Vx[TNum-1] = %g\n",Vx[TNum-1]);
  printf("------checking-----: Vy[0] = %g\n",Vy[0]);
  printf("------checking-----: Vy[TNum-1] = %g\n",Vy[TNum-1]);
  printf("------checking-----: Vz[0] = %g\n",Vz[0]);
  printf("------checking-----: Vz[TNum-1] = %g\n",Vz[TNum-1]);
  printf("------checking-----: Q[0] = %g\n",Q[0]);
  printf("------checking-----: Q[TNum-1] = %g\n",Q[TNum-1]);
  printf("------checking-----: OmegaM = %g\n",OmegaM);
  printf("------checking-----: Lx[0] = %g\n",Lx[0]);
  printf("------checking-----: Lx[TNum-1] = %g\n",Lx[TNum-1]);
  printf("------checking-----: Ly[0] = %g\n",Ly[0]);
  printf("------checking-----: Ly[TNum-1] = %g\n",Ly[TNum-1]);
  printf("------checking-----: Lz[0] = %g\n",Lz[0]);
  printf("------checking-----: Lz[TNum-1] = %g\n",Lz[TNum-1]);
  printf("------checking-----: FxM[0] = %g\n",FxM[0]);
  printf("------checking-----: FxM[TNum-1] = %g\n",FxM[TNum-1]);
  printf("------checking-----: FyM[0] = %g\n",FyM[0]);
  printf("------checking-----: FyM[TNum-1] = %g\n",FyM[TNum-1]);
  printf("------checking-----: FzM[0] = %g\n",FzM[0]);
  printf("------checking-----: FzM[TNum-1] = %g\n",FzM[TNum-1]);
  printf("------checking-----: FRM[0] = %g\n",FRM[0]);
  printf("------checking-----: FRM[TNum-1] = %g\n",FRM[TNum-1]);
  printf("------checking-----: FRStarM[0] = %g\n",FRStarM[0]);
  printf("------checking-----: FRStarM[TNum-1] = %g\n",FRStarM[TNum-1]);
	
  printf("------checking-----: FxP[5] = %g\n",FxP[5]); // here MaX, MaY, and MaZ are zero, is given from the file so the remaining result is zero
  printf("------checking-----: FyP[311] = %g\n",FyP[311]);
  printf("------checking-----: FzP[265] = %g\n",FzP[265]);
  printf("------checking-----: FRStarM[185] = %g\n",FRStarM[185]);
  printf("------checking-----: FRStarP[355] = %g\n",FRStarP[355]);
  printf("------checking-----: FRM[171] = %g\n",FRM[191]);
  printf("------checking-----: FRP[335] = %g\n",FRP[335]);

  printf("------checking-----: Gamma = %g\n",Gamma);
  printf("------checking-----: DT = %g\n",DT);
  printf("------checking-----: Tint = %g\n",Tint);
  printf("------checking-----: ODT = %g\n",ODT);
  printf("------checking-----: Q[0] = %g\n",Q[0]);
  printf("------checking-----: Q[TNum-1] = %g\n",Q[TNum-1]);
  printf("------checking-----: A = %g\n",A);
  printf("------checking-----: Drou = %g\n",Drou);
  printf("------checking-----: Dp = %g\n",Dp);
  printf("------checking-----: Du = %g\n",Du);
  printf("------checking-----: Dv = %g\n",Dv);
  printf("------checking-----: Dw = %g\n",Dw);
  printf("------checking-----: Gu = %g\n",Gu);
  printf("------checking-----: Gv = %g\n",Gv);
  printf("------checking-----: Gw = %g\n",Gw);
  printf("------checking-----: VectorX = %g\n",VectorX);
  printf("------checking-----: VectorY = %g\n",VectorY);
  printf("------checking-----: VectorZ = %g\n",VectorZ);
  printf("------checking-----: VUN = %g\n",VUN);
  printf("------checking-----: VvN = %g\n",VvN);
  printf("------checking-----: VuN = %g\n",VuN);
  printf("------checking-----: DVX = %g\n",DVX);
  printf("------checking-----: DVY = %g\n",DVY);
  printf("------checking-----: DVZ = %g\n",DVZ);
  printf("------checking-----: Area[0] = %g\n",Area[0]);
  printf("------checking-----: Area[0] = %g\n",Area[DSNum-1]);
  printf("------checking-----: UnitNormal[0][0] = %g\n",UnitNormal[0][0]);
  printf("------checking-----: UnitNormal[DSNum-1][0] = %g\n",UnitNormal[DSNum-1][0]);
  printf("------checking-----: UnitNormal[0][1] = %g\n",UnitNormal[0][0]);
  printf("------checking-----: UnitNormal[DSNum-1][1] = %g\n",UnitNormal[DSNum-1][1]);
  printf("------checking-----: UnitNormal[0][2] = %g\n",UnitNormal[0][0]);
  printf("------checking-----: UnitNormal[DSNum-1][2] = %g\n",UnitNormal[DSNum-1][2]);


  
 /* printf("------checking-----: PBD1[34][1][1][1] = %g\n",PBD1[34][0][0][0]);
  printf("------checking-----: PBD1[155][1][1][1] = %g\n",PBD2[155][0][0][0]);
  printf("------checking-----: PBD1[238][1][1][1] = %g\n",PBD1[238][0][0][0]);
  
  printf("------checking-----: PB1[34][1][1] = %g\n",PB1[34][0][0]);
  printf("------checking-----: PB1[155][1][1] = %g\n",PB2[155][0][0]);
  printf("------checking-----: PB1[238][1][1] = %g\n",PB1[238][0][0]);
  
  printf("------checking-----: P1[34][1] = %g\n",P1[34][0]);
  printf("------checking-----: P1[155][1] = %g\n",P2[155][0]);
  printf("------checking-----: P1[238][1] = %g\n",P3[238][0]);
  printf("------checking-----: pF[78][1] = %g\n",pF[78][0]);*/

  printf("RMaTip = %4.4f,C_0 = %4.4f, R = %4.4f, OmegaR = %4.4f [rad/s], fR = %4.4f [hz], TR = %4.4f [s], OmegaM = %4.4f [rad/s]\n\n",RMaTip, C_0, R, OmegaR, fR, TR, OmegaM);
  printf("\n------------ Pressure Prediciton is DONE!----------------\n");
  
  
  

 
  
  /*--------------------freeing all used memories--------------------*/


  free(DSX); 
  free(DSY);  
  free(DSZ); 
  free(DSp); 
  free(DSrou);
  free(DSu);
  free(DSv);
  free(DSw);


  free(Nx);
  free(Ny);
  free(Nz);

  free(XO);
  free(YO);
  free(ZO);

  free_vector(OTime);

  
  free_vector(DataXR);
  free_vector(DataYR);
  free_vector(DataZR);
  free_vector(DOrX);
  free_vector(DOrY);
  free_vector(DOrZ);
  free_vector(DOr);
  free_vector(DORStar);
  free_vector(DOR);
  free_vector(DORStarX);
  free_vector(DORStarY);
  free_vector(DORStarZ);
  free_vector(DORX);
  free_vector(DORY);
  free_vector(DORZ);
  free_vector(RGamma);
  free_vector(Vx); 
  free_vector(Vy); 
  free_vector(Vz); 
  free_vector(Q); 
  free_vector(Lx);
  free_vector(Ly);
  free_vector(Lz);
  free_vector(FxM);
  free_vector(FyM); 
  free_vector(FzM); 
  free_vector(FxP); 
  free_vector(FyP); 
  free_vector(FzP);
  free_vector(FRStarM); 
  free_vector(FRStarP);
  free_vector(FRM); 
  free_vector(FRP);
  free_vector(OF);
  free_vector(Area);
  
  free_4Ddmatrix(PBD1,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD2,OSNum,BNum,DSNum);
  free_4Ddmatrix(PBD3,OSNum,BNum,DSNum);
  free_3Ddmatrix(PB1,OSNum,BNum);
  free_3Ddmatrix(PB2,OSNum,BNum);
  free_3Ddmatrix(PB3,OSNum,BNum); 
  free_dmatrix(P1,OSNum);
  free_dmatrix(P2,OSNum);
  free_dmatrix(P3,OSNum);
  free_dmatrix(pF,FNum);
  free_dmatrix(Op,FNum);
  //free_dmatrix(pXOYM,FNum);
  free_dmatrix(UnitNormal,DSNum);
  //free_dmatrix(DataSVector,DSNum);
  

  


  return 0;
}


