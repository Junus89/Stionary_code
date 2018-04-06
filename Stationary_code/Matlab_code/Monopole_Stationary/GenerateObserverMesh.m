clear
close all
clc
tic

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate the observer points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 20;

ThetaNum = 1;
PhiNum = 41;


% Theta = linspace(0,0,ThetaNum);   %Polar anlge
Phi = linspace(0,2*pi,PhiNum);  %Azimuthal angle

Theta = pi/2;
% Phi = 0;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

SCoordX = zeros(PhiNum-1,ThetaNum);
SCoordY = zeros(PhiNum-1,ThetaNum);
SCoordZ = zeros(PhiNum-1,ThetaNum);

for j = 1:ThetaNum 
    
    for k = 1:PhiNum-1
        
        SCoordX(k,j) = R*sin(Theta(j))*cos(Phi(k));
        SCoordY(k,j) = R*sin(Theta(j))*sin(Phi(k));
        SCoordZ(k,j) = R*cos(Theta(j));
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate the observer points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Coordinates Generation is done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Output the observer points %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('ObserverGeo.dat','w');

for j = 1:ThetaNum
    
    for k = 1:PhiNum-1
        
        fprintf(fid,'%e\t %e\t %e\n',SCoordX(k,j),SCoordY(k,j),SCoordZ(k,j));
        
    end
    
    fprintf(fid,'\n');
        
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Output the observer points %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

