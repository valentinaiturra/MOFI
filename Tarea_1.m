clc
clear all
close all

datos = readmatrix('2023_Tarea1.xlsx'); %Tambien sirve xlsread
%Presion (hPa), Temperatura (°C), HR(%)

%Pasamos los datos a matrices individuales
p = datos(:,1); % [hPa]
T = datos(:,2); % [°C]
HR = datos(:,3);

%1.a) Calcular la razon de mezcla. Entregar los valores en g/kg para los
%niveles estandar de presion

%Calculamos la presion de vapor saturada
for i=1:length(p) 
    e_s(i,:) = 6.112*exp((17.67*T(i))/(T(i)+243.5)); % [hPa]
end

%Calculamos la razon de mezcla de saturacion
for i=1:length(p) 
    r_s(i,:)=0.622*(e_s(i)/p(i)); % [g/Kg]
end

%Determinamos la razon de mezcla
for i = 1:length(p)
    r_v(i,:) = (HR(i)*r_s(i))/100; % [g/Kg]
end

%Nos piden los datos en los siguientes niveles:
presion=[1000,925,850,700,500,200];

%Buscamos las ubicaciones de esos niveles de presion
for i=1:length(presion) 
    fila(i) = find(p == (presion(i)));
end

%Finalmente los datos obtenidos son
for i = 1:length(presion)
    rv(i,:) = r_v(fila(i)); % [g/Kg]
end


%1.b) Calcular la altura geopotencial. Entregar los valores en m para los
%niveles estandar de presion
z= [];
z(1)=10; %[m]
Rd = 287; % [J/Kg*K]
g = 9.8;   % [m/s2]
T = T+273.15; %[K]
r_v = r_v/1000; %[kg/kg]
Tv = T.*(1+0.61*r_v); %[K]

for i=2:length(p)
    Tmv(i) = (Tv(i)+Tv(i-1))/2;
    z(i) = (Rd/g)*Tmv(i)*log(p(i-1)/p(i))+z(i-1);
end

plot(T,z)

clear all
%2. Calcular el coeficiente de expansión térmica en la atmósfera. 
% Dejar fijo el valor de razón de mezcla de vapor de agua en 5 g kg-1 , 
% y de presión en 1000 hPa. Mostrar, gráficamente, que el coeficiente de 
% expansión térmica en la atmósfera es proporcional a T-1.

%Notemos que nuestros datos son:
p = 200:10:1000; % [hPa]
r_v = 0:1:20; % [g/Kg]
T = -70:1:30; % [°C]

%Determinamos coeficientes y pasamos al sistema internacional
Rd = 287; % [J/Kg*K]
T = T+273.15; % [K]
p = p*100; % [Pa]
r_v = r_v/1000; % [Kg/Kg]

for h = 1:length(p)
    for i = 1:length(T)
        for j = 1:length(r_v)
            Tv = (T(i))*(1+0.61*r_v(j));
            rho(h,i,j) = p(h)/(Rd*Tv);
        end
    end
end

%Notemos que rho esta en funcion de rho(p,T,r) y que nos dicen que r y p
%tienen valores fijos, por tanto:
rv_fijo = find(r_v == 5/1000);
p_fijo = find(p == 1000*100);

%con esta informacion tenemos entonces que 
rhoT = rho(p_fijo,:,rv_fijo);

%Entonces determinamos la derivada
drho_dT = rhoT(2:end)-rhoT(1:end-1);


for i = 2:length(rhoT)
    ipshylon(i-1) =  - drho_dT(i-1) / rhoT(i);
end

figure(1)
plot(T(2:end), ipshylon)