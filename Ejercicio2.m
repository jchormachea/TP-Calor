%% Ejercicio2 TP HT Hormachea 61439 - Nieto 61459
% 1D t(i+1) -2Ti+T(i-1) = 0
% 2D -T(i+1,j)-T(i-1,j)-T(i,j+1)-T(i,j-1)+4*T(i,j) = 0
clear; clc; close all
tic
%elección de refinado
nVolumesLength = 30; %Volumenes en longitud(minimo 3)
nVolumesHeight = 30; %Volumenes en altura(minimo 3)
nVolumes = nVolumesLength*nVolumesHeight;

%declaración de variables
L = 1; W = 1; %m
K = 1; %W/mK
qin = 30; %W/m2
T1 = 20; %°C

dx = L/nVolumesLength;
dy = W/nVolumesHeight;

Qt = sparse(nVolumes,nVolumes); %temperatures equation matrix

% numeración de los vols de control
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
    volumeSouth = iVolume-1;
    volumeNorth = iVolume+1;
    volumeWest = iVolume-nVolumesHeight;
    volumeEast = iVolume+nVolumesHeight;
    if ix == 1 %pared Oeste Tfija Ti-1,j = T1
        if iy == nVolumesHeight %cara norte-flujo fijo
            Qt(iVolume,iVolume) = dx^2+3*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeEast) = -dy^2;

        elseif iy == 1 %cara sur Ti,j-1 fija                
            Qt(iVolume,iVolume) = 3*dx^2+3*dy^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeEast) = -dy^2;

        else
            Qt(iVolume,iVolume) = 2*dx^2+3*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeEast) = -dy^2;      
        end

    elseif ix == nVolumesLength %pared Este Tfija Ti+1,j = T1

        if iy == nVolumesHeight %cara norte-flujo fijo

            Qt(iVolume,iVolume) = dx^2+3*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;

        elseif iy == 1 %cara sur Ti,j-1 fija                
            Qt(iVolume,iVolume) = 3*dx^2+3*dy^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;

        else
            Qt(iVolume,iVolume) = 2*dx^2+3*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;      
        end

    elseif iy == 1
        if ix == nVolumesHeight %cara este-T fija                
            Qt(iVolume,iVolume) = 3*dx^2+3*dy^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;

        elseif ix == 1 %cara oeste Ti,j-1 fija                
            Qt(iVolume,iVolume) = 3*dx^2+3*dy^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeEast) = -dy^2;

        else
            Qt(iVolume,iVolume) = 3*dx^2+2*dy^2;
            Qt(iVolume,volumeNorth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;
            Qt(iVolume,volumeEast) = -dy^2;      
        end

    elseif iy == nVolumesHeight

        if ix == nVolumesHeight %cara este-T fija                
            Qt(iVolume,iVolume) = 2*dx^2+2*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;

        elseif ix == 1 %cara oeste Ti,j-1 fija                
            Qt(iVolume,iVolume) = 2*dx^2+2*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeEast) = -dy^2;

        else
            Qt(iVolume,iVolume) = dx^2+2*dy^2;
            Qt(iVolume,volumeSouth) = -dx^2;
            Qt(iVolume,volumeWest) = -dy^2;
            Qt(iVolume,volumeEast) = -dy^2;      
        end

    else
        %caso normal
        
        Qt(iVolume,iVolume) = 2*dx^2+2*dy^2;
        Qt(iVolume,volumeNorth) = -dx^2;
        Qt(iVolume,volumeSouth) = -dx^2;
        Qt(iVolume,volumeEast) = -dy^2;
        Qt(iVolume,volumeWest) = -dy^2;            
    end
end

% full(Qt)
%boundary conditions
B = sparse(nVolumes,1); %el elemento B(i) corresponde a la ecuacion del vol de control en Qt(i,i)
sideWest = 1:nVolumesHeight;
B(sideWest) = B(sideWest)+dy^2*T1*2; %pared Oeste T fija
sideSouth = 1:nVolumesHeight:nVolumes;
B(sideSouth) = B(sideSouth)+dx^2*T1*2; %pared Sur T fija
sideNorth = nVolumesHeight:nVolumesHeight:nVolumes;
B(sideNorth) = B(sideNorth)+qin*dx^2*dy/K; %pared Norte flujo fijo
sideEast = nVolumes-nVolumesHeight+1:nVolumes;
B(sideEast) = B(sideEast)+dy^2*T1*2; %pared Este T fija

%solver
T = Qt\B;
T = full(T);

solucionNumerica = toc;

%% Resolucion analitica


[X, Y] = meshgrid(0.5*dx:dx:(L-0.5*dx),0.5*dy:dy:(W-0.5*dy));

Tan = 0; %temperatura analítica
for n = 1:100
    Tan = Tan+2*qin*L/(K*pi^2)*((1+(-1)^(n+1))/(n^2*cosh(n*pi*W/L))).*sin(n*pi*X/L).*sinh(n*pi*Y/L); 
end
Tan = Tan+T1;
solucionAnalitica = toc;

%% convergencia 
 
% numeración de los vols de control
% 3| 6| 9
% 2| 5| 8
% 1| 4| 7

ncol = ceil(nVolumesLength/2); %fix(nVolumesLength/2)+mod(nVolumesLength,2);
centerLine = (ncol-1)*nVolumesHeight+1:1:((ncol-1)*nVolumesHeight+nVolumesHeight); %volumenes de la col del centro
ECM = 0;
TECM = T(centerLine); %Temperaturas de MATLAB
sz = size(Tan);
y = fix(sz(1)/2)+mod(sz(1),2);%col del medio de Tan
TanECM = Tan(1:sz(2),y);
ECM = sqrt(sum((TECM-TanECM).^2)/nVolumesHeight);

%% flujos en extremos

%flujo sur
southLine = 1:nVolumesHeight:nVolumes; %volumenes cara sur
westLine = 1:nVolumesHeight;
qsouth = K/(0.5*dy)*(T(southLine)-T1);
qwest = K/(0.5*dx)*(T(westLine)-T1);

%% post procesado para plots
Tfield = reshape(full(T),nVolumesHeight, nVolumesLength); %para el contourf de temperatura
Tsup = qin*dy*0.5+K*T(nVolumesHeight:nVolumesHeight:nVolumes); %Tsuperior MATLAB
TsupAn = qin*0.5*dy+K*Tan(end,:);%Tsuperior Analitica
T1vec = ones(1,nVolumesLength)*T1; 
Tfield = [T1vec;Tfield];
Tfield = [Tfield;Tsup'];
Tfield = [[T1vec T1 T1]',Tfield];
Tfield = [Tfield, [T1vec T1 T1]'];
Tan = [T1vec;Tan];
Tan = [Tan;TsupAn];
Tan = [[T1vec T1 T1]',Tan];
Tan = [Tan, [T1vec T1 T1]'];
[Xplot, Yplot] = meshgrid([0 0.5*dx:dx:(L-0.5*dx) L],[0 0.5*dy:dy:(W-0.5*dy) W]);

%% plots


%error cuadrático medio
ECM
nVolumes

ECMval = [9 0.4203; 16 0.2521; 25 0.1563; 36 0.1112; 49 0.0805; 81 0.0488; 144 0.0277; 225 0.0177; 400 0.01; 625 0.0064]; %900 0.0044];
plot(ECMval(:,1),ECMval(:,2))
grid on; hold on
xlabel('Número de Volúmenes')
ylabel('ECM')
title('Error cuadrático Medio')
xticks(ECMval(:,1))
xlim([0 ECMval(end,1)])

%flujo de calor Oeste
figure
plot(linspace(0,W,nVolumesHeight),qwest,'r')
grid on
xlabel('posicion [m]')
ylabel('flujo de calor [W/m^2]')
title('flujo de calor en la pared Oeste')

%flujo de calor Sur
figure
plot(linspace(0,L,nVolumesLength),qsouth,'r')
grid on
xlabel('posicion [m]')
ylabel('flujo de calor [W/m^2]')
title('flujo de calor en la pared Sur')


%campo de temperatura numérico
figure
contourf(Xplot,Yplot,Tfield,'linecolor','none')
colormap('jet')
c = colorbar;
ylabel(c,'Temp (°C)','Rotation', 270)
c.Label.Position(1) = 4;
title('Campo de Temperaturas numérico')

%campo de temperatura analítico
figure
contourf(Xplot,Yplot,Tan,'linecolor','none')
title('Campo de Temperaturas Analítico')
colormap('jet')
c = colorbar;
ylabel(c,'Temp (°C)','Rotation', 270)
c.Label.Position(1) = 4;