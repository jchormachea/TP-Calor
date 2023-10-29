%% Ejercicio1 TP HT Hormachea 61439 - Nieto 61459
% Este código resuelve el ejercicio 1 del tp de transferencia de calor.
% Este código fue realizado con matlab R2020a, el uso de otra puede generar que no corra. 
%% incializar
clear; clc; close all

%% preprocesado

%eleccion de refinado
nVolumes = 10; % el minimo es 3

%declaración de variables
Tamb = 25; To = 0; %°C
L = 1; A = 0.1; %[m] | [m2]
K = 1; h = 11; %[W/mK] |[W/m2K]
q = 25; %[W/m3]
r = sqrt(A/pi);%[m]


dx = L/nVolumes; %[m]
As = 2*pi*r*dx; %area superficial[m2]
V = pi*r^2*dx; %[m^3]volumen de cada volumen finito

% armado de matriz de ecuaciones
Qt = sparse(nVolumes,nVolumes); %temperatures equation matrix

for iVol = 1:nVolumes
    
   if iVol == 1
       Qt(iVol,iVol) = -3*K*A-h*As*dx;
       Qt(iVol,iVol+1) = K*A;
   elseif iVol == nVolumes
       Qt(iVol,iVol) = -K*A-h*As*dx; 
       Qt(iVol,iVol-1) = K*A;  
   else
       Qt(iVol,iVol) = -2*K*A-h*As*dx;
       Qt(iVol,iVol+1) = K*A;
       Qt(iVol,iVol-1) = K*A;
   end    
end

%boundary conditions
B = sparse(nVolumes,1);
B(1) = -2*K*A*To; %cond de borde temperatura fija
B(end) = 0; %cond de borde flujo nulo(aislado)
B(:) = B(:)-q*V*dx; %Generación de energía interna
B(:) = B(:)-h*As*Tamb*dx; %Convección de calor

%% Solver
T = Qt\B;
T = full(T);
T = [To;T;T(end)]; %agrego las puntas de la barra

%% Post procesado
%solucion teórica
x = [0 0.5*dx:dx:(L-0.5*dx) L];
P = 2*pi*r; %perímetro
m = sqrt((h*P)/(K*A));
Tteo = (((To-Tamb-q/(K*m^2))/(1+exp(2*m*L)))*(exp(m*x)+exp(2*m*L)*exp(-m*x))+q/(K*m^2)+Tamb)';

% flujo de calor
qf(nVolumes) = 0; %flujo de calor [W/m^2]

for i = 2:size(T,1)
    
    if i == 1
        qf(i) = K/(0.5*dx)*(T(i)-T(i-1));
    else
        qf(i) = K/dx*(T(i)-T(i-1));
    end
end


%% Print de datos
fprintf('Ejercicio 1.\n')
fprintf('nivel de refinamiento: %d\n', nVolumes)
fprintf('Temperaturas: \n')

figure
plot(x,T,'b')
hold on; grid on
plot(x,Tteo,'r-.')
legend('FVM','Analitica')
title('Comparación de soluciones')
xlabel('posicion [m]')
ylabel('Temperatura [°C]')

figure
plot(x,qf,'r')
grid on
title('Flujo de calor')
xlabel('posicion [m]')
ylabel('Flujo [W/m^2]')