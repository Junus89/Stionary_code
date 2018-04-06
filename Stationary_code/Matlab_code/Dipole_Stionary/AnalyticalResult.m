clear
close all
clc
tic

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Analytical Codes for Stationary Sources in a Mean Flow %%%%%%%
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

ObserverS = importdata('ObserverGeo.dat');

XO = ObserverS(:,1);
YO = ObserverS(:,2);
ZO = ObserverS(:,3);

OSNum = length(XO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = sqrt(XO.^2+(YO+0.5*d).^2+ZO.^2);
r2 = sqrt(XO.^2+(YO-0.5*d).^2+ZO.^2);
RStar1 = sqrt(r1.^2+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO).^2)/Gama;
RStar2 = sqrt(r2.^2+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO).^2)/Gama;
R1 = Gama^2*(RStar1-(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO));
R2 = Gama^2*(RStar2-(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO));

VPFCoefX11 = (XO+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaX)./(Gama^2*RStar1);
VPFCoefY11 = ((YO+0.5*d)+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaY)./(Gama^2*RStar1);
VPFCoefZ11 = (ZO+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaZ)./(Gama^2*RStar1);
VPFCoefX21 = (XO+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaX)./(Gama^2*RStar2);
VPFCoefY21 = ((YO-0.5*d)+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaY)./(Gama^2*RStar2);
VPFCoefZ21 = (ZO+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaZ)./(Gama^2*RStar2);

VPFCoefX12=(1i*Omega/C0)*(XO+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaX...
           -Gama^2*RStar1*MaX);
VPFCoefY12=(1i*Omega/C0)*((YO+0.5*d)+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaY...
           -Gama^2*RStar1*MaY);
VPFCoefZ12=(1i*Omega/C0)*(ZO+Gama^2*(MaX*XO+MaY*(YO+0.5*d)+MaZ*ZO)*MaZ...
           -Gama^2*RStar1*MaZ);
VPFCoefX22=(1i*Omega/C0)*(XO+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaX...
           -Gama^2*RStar2*MaX);
VPFCoefY22=(1i*Omega/C0)*((YO-0.5*d)+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaY...
           -Gama^2*RStar2*MaY);
VPFCoefZ22=(1i*Omega/C0)*(ZO+Gama^2*(MaX*XO+MaY*(YO-0.5*d)+MaZ*ZO)*MaZ...
           -Gama^2*RStar2*MaZ);
       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

VPF1 = zeros(OSNum,TNum);
VPF2 = zeros(OSNum,TNum);
VPFCoef1 = zeros(OSNum,TNum);
VPFCoef2 = zeros(OSNum,TNum);
VPFX1 = zeros(OSNum,TNum);
VPFY1 = zeros(OSNum,TNum);
VPFZ1 = zeros(OSNum,TNum);
VPFX2 = zeros(OSNum,TNum);
VPFY2 = zeros(OSNum,TNum);
VPFZ2 = zeros(OSNum,TNum);

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

disp('Pressure Prediction is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMSpXOY = zeros(OSNum+1,1);

for j = 1:OSNum
    
    RMSpXOY(j) = norm(pXOY(j,:))/sqrt(TNum);
    
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

plot(Time/T,pXOY(11,:),'r-o');

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

for j = 1:OSNum
    
    File = ['Ana_Result/AnaTimePressure',num2str(j)];
    fid = fopen(File,'w');
    for k = 1:TNum
        
        fprintf(fid,'%e\t %e\n',Time(k)/T,pXOY(j,k));
        
    end
    
    fclose(fid);
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

fid = fopen('AnaDirectivity.dat','w');

for j = 1:OSNum+1
    
    fprintf(fid,'%e\t %e\n',ObserverSTheta(j),RMSpXOY(j));
    
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Output is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Analytical Codes for Stationary Sources in a Mean Flow %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
