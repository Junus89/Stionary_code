clear
close all
clc
tic
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FWHFD Codes for Stationary Sources in a Mean Flow %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultLineLineWidth',1.1)
set(0,'DefaultaxesLineWidth',1)
set(0,'DefaultaxesFontSize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlowData = importdata('FlowParameter.dat');

Rou0 = FlowData(1);                                                        %[kg/m^3]
C0 = FlowData(2);                                                          %[m/s]

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

ComData = importdata('ComputationParameter.dat');

Omega = ComData(1);
T = 2*pi/Omega;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

SamData = importdata('SampleParameter.dat');

TNum = SamData(1);
DT = T/TNum;
Time = DT*(0:TNum-1);

TimeHalf(1:TNum/2) = Time(2*(1:TNum/2)-1);
OmegaHalf = zeros(TNum/2,1);
KHalf = zeros(TNum/2,1);

for j = 1:TNum/2
    
    OmegaHalf(j,1) = Omega*(j-1);
    if(OmegaHalf(j,1)==0)
        
        OmegaHalf(j,1)=10^-12;
        
    end
    
    KHalf(j,1)=OmegaHalf(j)/C0;
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

MMa = importdata('MachParameter.dat');

MaX = MMa(1);
MaY = MMa(2);
MaZ = MMa(3);
Gama = sqrt(1/(1-(MaX^2+MaY^2+MaZ^2)));

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

ObserverS = importdata('ObserverGeo.dat');

XO = ObserverS(:,1);
YO = ObserverS(:,2);
ZO = ObserverS(:,3);

OSNum = length(XO);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

DataS = importdata('SourceGeo.dat');
DataX = DataS(:,1);
DataY = DataS(:,2);
DataZ = DataS(:,3);
DataNX = DataS(:,4);
DataNY = DataS(:,5);
DataNZ = DataS(:,6);

DSNum = length(DataX);

DataSArea = zeros(DSNum,1);
DataSVector = zeros(DSNum,3);

for m = 1:DSNum
    
    DataSArea(m) = sqrt((DataNX(m))^2+(DataNY(m))^2+(DataNZ(m))^2) ;
    DataSVector(m,1) = DataNX(m)/DataSArea(m) ;
    DataSVector(m,2) = DataNY(m)/DataSArea(m) ;
    DataSVector(m,3) = DataNZ(m)/DataSArea(m) ;

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

Nvar=5;                                                                    % Number of flow varaibles
        
SourceData = importdata('SourceData.dat');
SourceData3d = reshape(SourceData',[Nvar TNum DSNum]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOrX = zeros(DSNum,OSNum);
DOrY = zeros(DSNum,OSNum);
DOrZ = zeros(DSNum,OSNum);

for j = 1:DSNum
    
    for k = 1:OSNum
        
        DOrX(j,k) = XO(k)-DataX(j);
        DOrY(j,k) = YO(k)-DataY(j);
        DOrZ(j,k) = ZO(k)-DataZ(j);
        
    end
    
end

DOr = sqrt(DOrX.^2+DOrY.^2+DOrZ.^2);
DORStar = sqrt(DOr.^2+Gama^2*(MaX*DOrX+MaY*DOrY+MaZ*DOrZ).^2)/Gama;
DOR = Gama^2*(DORStar-(MaX*DOrX+MaY*DOrY+MaZ*DOrZ));
DORX = (DOr./DORStar).*(DOrX./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaX)-Gama^2*MaX;
DORY = (DOr./DORStar).*(DOrY./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaY)-Gama^2*MaY;
DORZ = (DOr./DORStar).*(DOrZ./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaZ)-Gama^2*MaZ;
DORStarX = (DOr./DORStar).*(DOrX./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaX)/Gama^2;
DORStarY = (DOr./DORStar).*(DOrY./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaY)/Gama^2;
DORStarZ = (DOr./DORStar).*(DOrZ./DOr+Gama^2*(MaX*DOrX./DOr+MaY*DOrY./DOr+MaZ*DOrZ./DOr)...
        *MaZ)/Gama^2; 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

QT = zeros(TNum,DSNum);
FTX = zeros(TNum,DSNum);
FTY = zeros(TNum,DSNum);
FTZ = zeros(TNum,DSNum);
%
for j = 1:DSNum
    
    SourceData2d = SourceData3d(:,:,j)';
    
%==========================================================================
% First column: pressure fluctuation.
% Second column: density fluctuation.
% Third column: x velocity fluctuation.
% Fourth column: y velocity fluctuation.
% Fifth column: z velocity fluctuation.
%==========================================================================

    for k = 1:TNum
        
        QT(k,j) = (SourceData2d(k,2)+Rou0)*((SourceData2d(k,3)+MaX*C0)...
                *DataSVector(j,1)+(SourceData2d(k,4)+MaY*C0)*...
                DataSVector(j,2)+(SourceData2d(k,5)+MaZ*C0)*DataSVector(j,3));
        FTX(k,j) = (SourceData2d(k,2)+Rou0)*SourceData2d(k,3)*((SourceData2d(k,3)+MaX*C0)...
                *DataSVector(j,1)+(SourceData2d(k,4)+MaY*C0)*...
                DataSVector(j,2)+(SourceData2d(k,5)+MaZ*C0)*DataSVector(j,3))...
                +SourceData2d(k,1)*DataSVector(j,1)-MaX*C0*QT(k,j);
        FTY(k,j) = (SourceData2d(k,2)+Rou0)*SourceData2d(k,4)*((SourceData2d(k,3)+MaX*C0)...
                *DataSVector(j,1)+(SourceData2d(k,4)+MaY*C0)*...
                DataSVector(j,2)+(SourceData2d(k,5)+MaZ*C0)*DataSVector(j,3))...
                +SourceData2d(k,1)*DataSVector(j,2)-MaY*C0*QT(k,j);
        FTZ(k,j) = (SourceData2d(k,2)+Rou0)*SourceData2d(k,5)*((SourceData2d(k,3)+MaX*C0)...
                *DataSVector(j,1)+(SourceData2d(k,4)+MaY*C0)*...
                DataSVector(j,2)+(SourceData2d(k,5)+MaZ*C0)*DataSVector(j,3))...
                +SourceData2d(k,1)*DataSVector(j,3)-MaZ*C0*QT(k,j);
        
    end
    
    QT(:,j) = QT(:,j)-mean(QT(:,j));
    FTX(:,j) = FTX(:,j)-mean(FTX(:,j));
    FTY(:,j) = FTY(:,j)-mean(FTY(:,j));
    FTZ(:,j) = FTZ(:,j)-mean(FTZ(:,j));
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

QF = fft(QT);
FFX = fft(FTX);
FFY = fft(FTY);
FFZ = fft(FTZ);
QFH = QF(1:TNum/2,:);
FFXH = FFX(1:TNum/2,:);
FFYH = FFY(1:TNum/2,:);
FFZH = FFZ(1:TNum/2,:);

GQ = zeros(TNum/2,DSNum);
GFRX = zeros(TNum/2,DSNum);
GFRY = zeros(TNum/2,DSNum);
GFRZ = zeros(TNum/2,DSNum);
GFRStarX = zeros(TNum/2,DSNum);
GFRStarY = zeros(TNum/2,DSNum);
GFRStarZ = zeros(TNum/2,DSNum);
pMonoF = zeros(TNum/2,OSNum);
pDiFR = zeros(TNum/2,OSNum);
pDiFRStar = zeros(TNum/2,OSNum);

for j = 1:OSNum
    
    EKR = exp(-1i*KHalf*(DOR(:,j)'));
    
    for k = 1:TNum/2
        
        GQ(k,:) = (1i*OmegaHalf(k)*EKR(k,:).*DataSArea(:,1)')./(4*pi*DORStar(:,j)');
        GFRX(k,:) = ((1i*KHalf(k)*EKR(k,:).*DORX(:,j)').*DataSArea(:,1)')./(4*pi*DORStar(:,j)');
        GFRY(k,:) = ((1i*KHalf(k)*EKR(k,:).*DORY(:,j)').*DataSArea(:,1)')./(4*pi*DORStar(:,j)');
        GFRZ(k,:) = ((1i*KHalf(k)*EKR(k,:).*DORZ(:,j)').*DataSArea(:,1)')./(4*pi*DORStar(:,j)');
        GFRStarX(k,:) = ((EKR(k,:).*DORStarX(:,j)').*DataSArea(:,1)')./(4*pi*(DORStar(:,j)').^2);
        GFRStarY(k,:) = ((EKR(k,:).*DORStarY(:,j)').*DataSArea(:,1)')./(4*pi*(DORStar(:,j)').^2);
        GFRStarZ(k,:) = ((EKR(k,:).*DORStarZ(:,j)').*DataSArea(:,1)')./(4*pi*(DORStar(:,j)').^2);
        
    end
    
    pMonoF_M = QFH.*GQ;
    pDiFR_M = FFXH.*GFRX+FFYH.*GFRY+FFZH.*GFRZ;
    pDiFRStar_M = FFXH.*GFRStarX+FFYH.*GFRStarY+FFZH.*GFRStarZ;
    
    for k = 1:TNum/2
        
        pMonoF(k,j) = sum(pMonoF_M(k,:));
        pDiFR(k,j) = sum(pDiFR_M(k,:));
        pDiFRStar(k,j) = sum(pDiFRStar_M(k,:));
        
    end
    
end

pMonoT = real(ifft(pMonoF));
pDiTR = real(ifft(pDiFR));
pDiTRStar = real(ifft(pDiFRStar));
pObser = pMonoT+pDiTR+pDiTRStar;
pXOY = pObser;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Pressure Prediction is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMSpXOY = zeros(OSNum+1,1);

for j = 1:OSNum
    
    RMSpXOY(j) = norm(pXOY(:,j))/sqrt(TNum/2);
    
end

ObserverSTheta = linspace(0,2*pi,OSNum+1);     %[rad]
RMSpXOY(OSNum+1) = RMSpXOY(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Postprocess is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
grid on
hold on
box on

plot(TimeHalf/T,pXOY(:,1),'r-o');

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

figure(2)
grid on
box on

polar(ObserverSTheta',RMSpXOY,'r-o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Plot is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:OSNum
   
    File = ['FWH_Result/FWHTimePressure',num2str(k)];
    fid = fopen(File,'w');
    
    for j = 1:TNum/2
        
        fprintf(fid,'%e\t %e\n',TimeHalf(j)/T,pXOY(j,k));
        
    end
    
    fclose(fid);
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

fid = fopen('FWHDirectivity.dat','w');

for j = 1:OSNum+1
    
    fprintf(fid,'%e\t %e\n',ObserverSTheta(j),RMSpXOY(j));
    
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FWHFD Codes for Stationary Sources in a Mean Flow %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
