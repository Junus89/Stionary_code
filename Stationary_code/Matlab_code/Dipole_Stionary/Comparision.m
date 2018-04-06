clear
close all
clc
tic

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultLineLineWidth',1.1)
set(0,'DefaultaxesLineWidth',1)
set(0,'DefaultaxesFontSize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MMa = importdata('MachParameter.dat');

MaX = MMa(1);
MaY = MMa(2);
MaZ = MMa(3);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

ComData = importdata('ComputationParameter.dat');

Omega = ComData(1);
T = 2*pi/Omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%Comparision of pressure time history.
ANA = importdata('Ana_Result/AnaTimePressure11');
FWH = importdata('FWH_Result/FWHTimePressure11');
%
figure(1)
grid on
hold on
box on

plot(ANA(:,1),ANA(:,2),'k-')
plot(FWH(:,1),FWH(:,2),'ro')

legend('Analytic','Predicted')

xlabel('t/T')
ylabel('Acoustic Pressure [Pa]')

axis([0 1 -0.003 0.003])

set(gcf, 'PaperPositionMode','Auto')   % Use screen size

Filename1 = ['timehistory_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename1,'-depsc'); 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

DT = importdata('AnaDirectivity.dat');
DP = importdata('FWHDirectivity.dat');
DF = 1.2*max(DT(:,2))*ones(length(DT),1);

figure(3)
grid on

p = polar(DT(:,1),DF,'w');
set(p,'visible','off')
hold on
polar(DT(:,1),DT(:,2),'k-')
hold on
polar(DP(:,1),DP(:,2),'ro')
hold off

% legend('Time Domain','Frequency Domain','location','southeast')

set(gcf, 'PaperPositionMode','Auto')   % Use screen size
Filename3 = ['directivity_',num2str(10*MaX),num2str(10*MaY),num2str(10*MaZ)];
print(Filename3,'-depsc'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

