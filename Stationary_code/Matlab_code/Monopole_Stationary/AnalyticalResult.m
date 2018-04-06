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

ObserverS = importdata('ObserverGeo.dat');

XO = ObserverS(:,1);
YO = ObserverS(:,2);
ZO = ObserverS(:,3);

OSNum = length(XO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = sqrt(XO.^2+YO.^2+ZO.^2);
RStar = sqrt(r.^2+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO).^2)/Gama;
R = Gama^2*(RStar-(MaX*XO+MaY*YO+MaZ*ZO));

VPFCoefX1 = (XO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaX)./(Gama^2*RStar);
VPFCoefY1 = (YO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaY)./(Gama^2*RStar);
VPFCoefZ1 = (ZO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaZ)./(Gama^2*RStar);

VPFCoefX2 = (1i*Omega/C0)*(XO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaX...
           -Gama^2*RStar*MaX);
VPFCoefY2 = (1i*Omega/C0)*(YO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaY...
           -Gama^2*RStar*MaY);
VPFCoefZ2 = (1i*Omega/C0)*(ZO+Gama^2*(MaX*XO+MaY*YO+MaZ*ZO)*MaZ...
           -Gama^2*RStar*MaZ);
       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

VPF=zeros(OSNum,TNum);
VPFCoef=zeros(OSNum,TNum);
VPFX=zeros(OSNum,TNum);
VPFY=zeros(OSNum,TNum);
VPFZ=zeros(OSNum,TNum);

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

plot(Time/T,pXOY(1,:),'r-o');

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
