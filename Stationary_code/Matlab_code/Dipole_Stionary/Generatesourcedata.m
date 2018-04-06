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

A = 100;    %Source strengh [m^2/s]
d = 0.01;   %[m]

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

SamData = importdata('SampleParameter.dat');

TNum = SamData(1);
DT = T/TNum;
Time = DT*(0:TNum-1);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

MMa = importdata('MachParameter.dat');

MaX = MMa(1);
MaY = MMa(2);
MaZ = MMa(3);
Gama = sqrt(1/(1-(MaX^2+MaY^2+MaZ^2)));

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

DataS = importdata('SourceGeo.dat');

XD = DataS(:,1);
YD = DataS(:,2);
ZD = DataS(:,3);

DSNum = length(XD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = sqrt(XD.^2+(YD+0.5*d).^2+ZD.^2);
r2 = sqrt(XD.^2+(YD-0.5*d).^2+ZD.^2);
RStar1 = sqrt(r1.^2+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD).^2)/Gama;
RStar2 = sqrt(r2.^2+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD).^2)/Gama;
R1 = Gama^2*(RStar1-(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD));
R2 = Gama^2*(RStar2-(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD));

VPFCoefX11 = (XD+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaX)./(Gama^2*RStar1);
VPFCoefY11 = ((YD+0.5*d)+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaY)./(Gama^2*RStar1);
VPFCoefZ11 = (ZD+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaZ)./(Gama^2*RStar1);
VPFCoefX21 = (XD+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaX)./(Gama^2*RStar2);
VPFCoefY21 = ((YD-0.5*d)+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaY)./(Gama^2*RStar2);
VPFCoefZ21 = (ZD+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaZ)./(Gama^2*RStar2);

VPFCoefX12=(1i*Omega/C0)*(XD+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaX...
           -Gama^2*RStar1*MaX);
VPFCoefY12=(1i*Omega/C0)*((YD+0.5*d)+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaY...
           -Gama^2*RStar1*MaY);
VPFCoefZ12=(1i*Omega/C0)*(ZD+Gama^2*(MaX*XD+MaY*(YD+0.5*d)+MaZ*ZD)*MaZ...
           -Gama^2*RStar1*MaZ);
VPFCoefX22=(1i*Omega/C0)*(XD+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaX...
           -Gama^2*RStar2*MaX);
VPFCoefY22=(1i*Omega/C0)*((YD-0.5*d)+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaY...
           -Gama^2*RStar2*MaY);
VPFCoefZ22=(1i*Omega/C0)*(ZD+Gama^2*(MaX*XD+MaY*(YD-0.5*d)+MaZ*ZD)*MaZ...
           -Gama^2*RStar2*MaZ);
       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

VPF1 = zeros(DSNum,TNum);
VPF2 = zeros(DSNum,TNum);
VPFCoef1 = zeros(DSNum,TNum);
VPFCoef2 = zeros(DSNum,TNum);
VPFX1 = zeros(DSNum,TNum);
VPFY1 = zeros(DSNum,TNum);
VPFZ1 = zeros(DSNum,TNum);
VPFX2 = zeros(DSNum,TNum);
VPFY2 = zeros(DSNum,TNum);
VPFZ2 = zeros(DSNum,TNum);

for j = 1:TNum
    
    VPF1(:,j) = A*exp(1i*Omega*(Time(j)-R1(:)/C0))./(4*pi*RStar1(:));
    VPF2(:,j) = A*exp(1i*Omega*(Time(j)-R2(:)/C0))./(4*pi*RStar2(:));
    
    VPFCoef1(:,j)=-A*exp(1i*Omega*(Time(j)-R1(:)/C0))./(4*pi*RStar1(:).^2);
    VPFCoef2(:,j)=-A*exp(1i*Omega*(Time(j)-R2(:)/C0))./(4*pi*RStar2(:).^2);    

    VPFX1(:,j)=VPFCoef1(:,j).*(VPFCoefX11(:)+VPFCoefX12(:));
    VPFY1(:,j)=VPFCoef1(:,j).*(VPFCoefY11(:)+VPFCoefY12(:));
    VPFZ1(:,j)=VPFCoef1(:,j).*(VPFCoefZ11(:)+VPFCoefZ12(:));
    VPFX2(:,j)=VPFCoef2(:,j).*(VPFCoefX21(:)+VPFCoefX22(:));
    VPFY2(:,j)=VPFCoef2(:,j).*(VPFCoefY21(:)+VPFCoefY22(:));
    VPFZ2(:,j)=VPFCoef2(:,j).*(VPFCoefZ21(:)+VPFCoefZ22(:));
    
end

realu = real(-VPFX1+VPFX2);
realv = real(-VPFY1+VPFY2);
realw = real(-VPFZ1+VPFZ2);
p1 = -Rou0*(1i*Omega*VPF1+C0*(MaX*VPFX1+MaY*VPFY1+MaZ*VPFZ1));
p2 = -Rou0*(1i*Omega*VPF2+C0*(MaX*VPFX2+MaY*VPFY2+MaZ*VPFZ2));
p = -p1+p2;
rou = p/C0^2;
realp = real(p);
realrou = real(rou);
pXOY=realp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Source Calculation is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('SourceData.dat','w');

for j = 1:DSNum
    
    for k = 1:TNum
        
        fprintf(fid,'%e\t %e\t %e\t %e\t %e\t\n',realp(j,k),realrou(j,k),realu(j,k),realv(j,k),realw(j,k));
        
    end
    
        fprintf(fid,'\n');
        
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Output is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FWHFD Codes for Stationary Sources in a Mean Flow %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
