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

A=1;                                                                       %Source strengh [m^2/s]

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

r = sqrt(XD.^2+YD.^2+ZD.^2);
RStar = sqrt(r.^2+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD).^2)/Gama;
R = Gama^2*(RStar-(MaX*XD+MaY*YD+MaZ*ZD));

VPFCoefX1 = (XD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaX)./(Gama^2*RStar);
VPFCoefY1 = (YD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaY)./(Gama^2*RStar);
VPFCoefZ1 = (ZD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaZ)./(Gama^2*RStar);

VPFCoefX2 = (1i*Omega/C0)*(XD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaX...
           -Gama^2*RStar*MaX);
VPFCoefY2 = (1i*Omega/C0)*(YD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaY...
           -Gama^2*RStar*MaY);
VPFCoefZ2 = (1i*Omega/C0)*(ZD+Gama^2*(MaX*XD+MaY*YD+MaZ*ZD)*MaZ...
           -Gama^2*RStar*MaZ);
       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

VPF=zeros(DSNum,TNum);
VPFCoef=zeros(DSNum,TNum);
VPFX=zeros(DSNum,TNum);
VPFY=zeros(DSNum,TNum);
VPFZ=zeros(DSNum,TNum);

for j = 1:TNum
    
    VPF(:,j) = A*exp(1i*Omega*(Time(j)-R(:)/C0))./(4*pi*RStar(:));
    VPFCoef(:,j)=-A*exp(1i*Omega*(Time(j)-R(:)/C0))./(4*pi*RStar(:).^2);
    
    VPFX(:,j) = VPFCoef(:,j).*(VPFCoefX1(:)+VPFCoefX2(:));
    VPFY(:,j) = VPFCoef(:,j).*(VPFCoefY1(:)+VPFCoefY2(:));
    VPFZ(:,j) = VPFCoef(:,j).*(VPFCoefZ1(:)+VPFCoefZ2(:));
    
end

realu = real(VPFX);
realv = real(VPFY);
realw = real(VPFZ);
p = -Rou0*(1i*Omega*VPF+C0*(MaX*VPFX+MaY*VPFY+MaZ*VPFZ));
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
