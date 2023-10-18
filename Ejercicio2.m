%% Ejercicio2 TP HT
% 1D t(i+1) -2Ti+T(i-1) = 0
% 2D -T(i+1,j)-T(i-1,j)-T(i,j+1)-T(i,j-1)+4*T(i,j) = 0
clear; clc; close all

%elección de refinado
nVolumesLength = 30; %Volumenes en longitud(minimo 3)
nVolumesHeight = 30; %Volumenes en altura(minimo 3)
nVolumes = nVolumesLength*nVolumesHeight;

%declaración de variables
L = 1; W = 1; %m
K = 1; %W/mK
qin = 30; %W/m2
T1 = 10; %°C

dx = L/nVolumesLength;
dy = L/nVolumesHeight;

Qt = sparse(nVolumes,nVolumes); %temperatures equation matrixcc

% numeración de los vols de control(creo que es así)
% 3| 6| 9
% 2| 5| 8
% 1| 4| 7

volumesCoordinates(nVolumes,2) = 0; % matriz de coordenadas de los volumenes
for iColumn = 1:nVolumesLength
    volumes = 1+nVolumesHeight*(iColumn-1):nVolumesHeight+nVolumesHeight*(iColumn-1);
    volumesCoordinates(volumes,2) = 1:nVolumesHeight;
    volumesCoordinates(volumes,1) = iColumn;
end

for iVolume = 1:nVolumes
    ix = volumesCoordinates(iVolume,1);
    iy = volumesCoordinates(iVolume,2);       
    if ix == 1 %pared Oeste Tfija Ti-1,j = T1
        if iy == nVolumesHeight %cara norte-flujo fijo

            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume+1,iVolume) = -dy^2;

        elseif iy == 1 %cara sur Ti,j-1 fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume+1,iVolume) = -dy^2;

        else
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume+1,iVolume) = -dy^2;      
        end

    elseif ix == nVolumesLength %pared Oeste Tfija Ti+1,j = T1

        if iy == nVolumesHeight %cara norte-flujo fijo

            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;

        elseif iy == 1 %cara sur Ti,j-1 fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;

        else
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;      
        end

    elseif iy == 1
        if ix == nVolumesHeight %cara este-T fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;

        elseif ix == 1 %cara oeste Ti,j-1 fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume+1,iVolume) = -dy^2;

        else
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume+1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;
            Qt(iVolume+1,iVolume) = -dy^2;      
        end

    elseif iy == nVolumesHeight

        if ix == nVolumesHeight %cara este-T fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;

        elseif ix == 1 %cara oeste Ti,j-1 fija                
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume+1,iVolume) = -dy^2;

        else
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,iVolume-1) = -dx^2;
            Qt(iVolume-1,iVolume) = -dy^2;
            Qt(iVolume+1,iVolume) = -dy^2;      
        end

    else
        %caso normal
        Qt(iVolume,iVolume) = 2*dx^2+2*dy^2;
        Qt(iVolume,iVolume+1) = -dx^2;
        Qt(iVolume,iVolume-1) = -dx^2;
        Qt(iVolume+1,iVolume) = -dy^2;
        Qt(iVolume-1,iVolume) = -dy^2;            
    end
end

% full(Qt)
% error('corte y confeccion')
%boundary conditions
B = sparse(nVolumes,1); %el elemento B(i) corresponde a la ecuacion del vol de control en Qt(i,i)
sideWest = 1:nVolumesHeight;
B(sideWest) = B(sideWest)+dy^2*T1; %pared Oeste T fija
sideSouth = 1:nVolumesHeight:nVolumes;
B(sideSouth) = B(sideSouth)+dx^2*T1; %pared Sur T fija
sideNorth = nVolumesHeight:nVolumesHeight:nVolumes;
B(sideNorth) = B(sideNorth)-qin*dx^2*dy/K; %pared Norte flujo fijo
sideEast = nVolumes-nVolumesHeight+1:nVolumes;
B(sideEast) = B(sideEast)+dy^2*T1; %pared Este T fija
% error('corte y confeccion')
%solver
T = Qt\B;
Tfield = reshape(full(T),nVolumesHeight, nVolumesLength);
figure
contourf(Tfield)
colormap('jet')
colorbar
title('Campo de Temperaturas numérico')

%% Resolucion analitica

[X, Y] = meshgrid(0:0.01:L, 0:0.01:W);
Tan = 0;
for n = 1:100
    Tan = Tan+2*qin*L/(K*pi^2)*((1+(-1)^(n+1))/(n^2*cosh(n*pi*W/L))).*sin(n*pi*X/L).*sinh(n*pi*Y/L); 
end
Tan = Tan+T1;
figure
contourf(X,Y,Tan)
title('Campo de Temperaturas Analítico')
colormap('jet')
colorbar