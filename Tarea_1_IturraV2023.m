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


for i=1:length(p)
    %Calculamos la presion de vapor saturada
        e_s(i,:) = 6.112*exp((17.67*T(i))/(T(i)+243.5)); % [hPa]
    %Calculamos la razon de mezcla de saturacion
        r_s(i,:)=0.622*(e_s(i)/(p(i)-e_s(i))); % [Kg/Kg]
    %Determinamos la razon de mezcla
        r_v(i,:) = (HR(i)*r_s(i))/100; % [Kg/Kg]
end


%Nos piden los datos en los siguientes niveles:
presion=[1000,925,850,700,500,200];


for i=1:length(presion)
    %Buscamos las ubicaciones de esos niveles de presion
        fila(i) = find(p == (presion(i)));
    %Finalmente los datos obtenidos son
        hola(i,6) = r_v(fila(i))*1000; % [g/Kg] rv
        hola(i,4) = e_s(fila(i));
        hola(i,5) = r_s(fila(i))*1000;
        hola(i,2) = T(fila(i));
        hola(i,3) = HR(fila(i));
        hola(i,1) = p(fila(i));
end

%1.b) Calcular la altura geopotencial. Entregar los valores en m para los
%niveles estandar de presion
z= [];
z(1) = 0; % [m]
Rd = 287; % [J/Kg*K]
g = 9.8;   % [m/s2]
T = T+273.15; % [K]
Tv = T.*(1+0.61*r_v); %[K]

Tmv(1)=Tv(1);
for i=2:length(p)
    Tmv(i,:) = (Tv(i)+Tv(i-1))/2; % [K]
    z(i,:) = (Rd/g)*Tmv(i)*log(p(i-1)/p(i))+z(i-1); %[m]
end


%Finalmente los datos obtenidos son
for i = 1:length(presion)
    zv(i,:) = z(fila(i)); % [m]
    Tmvv(i,:) = Tmv(fila(i)); % [k]
    pv(i,:) = p(fila(i));
end


figure(1)
subplot (1,2,1)
    plot(Tmv,z,'-b','LineWidth',1)
    hold on
    plot(Tmvv,zv,'o','LineWidth',2,'Color', 'r')
    hold off
    legend('Todos los datos', 'Para las presiones dadas')
    xlabel('Temperatura [K]')
    ylabel('Altura geopotencial [m]')
    axis tight
    grid minor
subplot (1,2,2)
    plot(p,z,'-b','LineWidth',1)
    hold on
    plot(pv,zv,'o','LineWidth',2,'Color','g')
    hold off
    legend('Todos los datos', 'Para las presiones dadas')
    xlabel('Presión [hPa]')
    ylabel('Altura geopotencial [m]')
    axis tight
    grid minor

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

figure(2)
plot(T(2:end), ipshylon,'-r','LineWidth',2)
xlabel('Temperatura [K]')
ylabel('Coeficiente de expansión térmica atmosférica')
axis tight
grid minor