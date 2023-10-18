%% prueba calor
% 1D t(i+1) -2Ti+T(i-1) = 0
% 2D -T(i+1,j)-T(i-1,j)-T(i,j+1)-T(i,j-1)+4*T(i,j) = 0
clear; clc; close all
%% Ejercicio 1
Tprueba = 25;
Tamb = 25; To = 0; %°C
L = 1; A = 0.1; %m | m2
K = 1; h = 11; %W/mK |W/m2K
q = 25; %W/m3
r = sqrt(A/pi);%m

nVolumes = 10000; %refinado

Qt = sparse(nVolumes,nVolumes); %temperatures equation matrix
dx = L/nVolumes; 
At = 2*pi*r*dx; %area transversal
V = pi*r^2*dx; %volumen de cada volumen finito


for iVol = 1:nVolumes
    
   if iVol == 1
       Qt(iVol,iVol) = -3*K*A-h*At;
       Qt(iVol,iVol+1) = K*A;
   elseif iVol == nVolumes
       Qt(iVol,iVol) = -K*A-h*At; 
       Qt(iVol,iVol-1) = K*A;  
   else
       Qt(iVol,iVol) = -2*K*A-h*At;
       Qt(iVol,iVol+1) = K*A;
       Qt(iVol,iVol-1) = K*A;
   end    
end

%boundary conditions
B = sparse(nVolumes,1);
B(1) = -2*K*A*To; %cond de borde temperatura fija
B(end) = 0; %cond de borde flujo nulo(aislado)
B(:) = B(:)-q*V*dx; %Generación de energía interna
B(:) = B(:)-h*At*Tamb; %Convección de calor

T = Qt\B;
fprintf('Ejercicio 1.\n')
fprintf('nivel de refinamiento: %d\n', nVolumes)
fprintf('Temperaturas: \n')
T(1)
T(end)


