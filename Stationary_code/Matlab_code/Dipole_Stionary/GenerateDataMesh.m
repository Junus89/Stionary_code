clear
close all
clc
tic

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate data surface points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 1;

ThetaNum = 51;
PhiNum = 81;

Theta = linspace(0,pi,ThetaNum);   %Polar anlge
Phi = linspace(0,2*pi,PhiNum);  %Azimuthal angle

DCoordXP = zeros(PhiNum,ThetaNum);
DCoordYP = zeros(PhiNum,ThetaNum);
DCoordZP = zeros(PhiNum,ThetaNum);

for j = 1:ThetaNum 
    
    for k = 1:PhiNum
        
        DCoordXP(k,j) = R*sin(Theta(j))*cos(Phi(k));
        DCoordYP(k,j) = R*sin(Theta(j))*sin(Phi(k));
        DCoordZP(k,j) = R*cos(Theta(j));
        
    end
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

kj = 0;

for j = 1:ThetaNum-1 
    
    for k = 1:PhiNum-1
        
        kj = kj+1;
        
        DCoordXC(kj) = (DCoordXP(k,j)+DCoordXP(k,j+1)+DCoordXP(k+1,j)+DCoordXP(k+1,j+1))/4;
        DCoordYC(kj) = (DCoordYP(k,j)+DCoordYP(k,j+1)+DCoordYP(k+1,j)+DCoordYP(k+1,j+1))/4;
        DCoordZC(kj) = (DCoordZP(k,j)+DCoordZP(k,j+1)+DCoordZP(k+1,j)+DCoordZP(k+1,j+1))/4;
        
        crossx1 = DCoordXP(k+1,j+1)-DCoordXP(k,j);
        crossy1 = DCoordYP(k+1,j+1)-DCoordYP(k,j);
        crossz1 = DCoordZP(k+1,j+1)-DCoordZP(k,j);

        crossx2 = DCoordXP(k+1,j)-DCoordXP(k,j+1);
        crossy2 = DCoordYP(k+1,j)-DCoordYP(k,j+1);
        crossz2 = DCoordZP(k+1,j)-DCoordZP(k,j+1);
        
        NXC(kj) = (crossy1*crossz2-crossy2*crossz1)/2;
        NYC(kj) = -(crossx1*crossz2-crossx2*crossz1)/2;
        NZC(kj) = (crossx1*crossy2-crossx2*crossy1)/2;
        
    end
    
end

DSNum = kj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate data surface points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Coordinates and Normals Generation is done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate data surface points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
grid off
axis off
surf(DCoordXP,DCoordYP,DCoordZP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Plot is Done !\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('SourceGeo.dat','w');

for j = 1:DSNum
    
    fprintf(fid,'%e\t %e\t %e\t %e\t %e\t %e\n',DCoordXC(j),DCoordYC(j),DCoordZC(j),NXC(j),NYC(j),NZC(j));
        
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Generate data surface points %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc