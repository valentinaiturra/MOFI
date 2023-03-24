clear all
close all
clc

load('Lavapie_272.mat');

% Definimos constantes
Rd =  287; % [J/K Kg]
g = 9.8;   % [m/s2]
z(1) =  10;  % altura de la estacion [m]
rkilo=r/1000; %Kg/Kg
T = T+273.15; % [K] pasamos la temperatura que nos dan en celsius a kelvin

%a) CAlcular la altura geopotencial para cada nivel de presion,
%considerando el aire seco y el aire humedo.
Pres = p*100;
Tv = T.*(1+0.61*rkilo); % [K] temperatura virtual

for i=2:1957
    z(i,1) = (Rd/g)*Tv(i)*log(p(i-1)/p(i))+z(i-1);
end

figure
plot(Pres,z,'-b')
xlabel('Presion')
ylabel('Altura geopotencial')


%b) Calcular el espesor 1000-500 hpa, usando la integracion del aire
%humedo. Comparar con el espesor obtenido a partir de la temperatura
%virtual promedio de la capa.

%Calculamos Tvm, la temperatura virtual promedio entre la capa de 1000hpa a
%500 hpa, buscando entre que posiciones la presion es 500 y 1000

press_500=find(p >= 500);
press_500 = press_500(end);

press_1000 = find(p >= 1000);

Tvm = (Tv(press_1000)+Tv(press_500))/2;

Espesor= ((Rd*Tvm)/g)*log(1000/500)


