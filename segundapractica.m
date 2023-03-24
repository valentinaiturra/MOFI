clc
clear all
close all

p=[200:10:1000];
r=[0:1:20];
T=[-70:1:30];

Rd=287;

for h = 1:length(p)
    for i = 1:length(T)
        for j = 1:length(r)
            Tv = (T(i)+273.15)*(1+0.61*r(j)/1000);
            rho(h,i,j) = p(h)*100/(Rd*Tv);
        end
    end
end

% rho(p,r) con T fijo
Tfijo = find(T == 0); % posicion de T = 0 °C

figure(1)
contourf(p,r,squeeze(rho(:,Tfijo,:))')
xlabel('Presión [Pa]')
ylabel('Razón de mezcla [kg/kg]')
% que significa que las isopicnas (lineas de igual densidad) se vean verticales?
 %La densidad es invariable cr a la razon de mezcla, es mas importante la
 %presion


% rho(T,r) con p fijo
pfijo = find(p == 1000); % posicion de p = 1000 hPa

 figure(2)
contourf(T,r,squeeze(rho(pfijo,:,:))') 
xlabel('Temperatura [K]')
ylabel('Razón de mezcla [kg/kg]')

% que significa que las isopicnas (lineas de igual densidad) se vean verticales?
%razon de mezcla juega un papel mucho menos importante que la temperatura
%al determinar densidad



%% Parte 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Tfijo = find(T == 0); % posicion de T = 0 °C
rfijo = find(r == 5); % posicion de r = 5 g/kg
rhoP = rho(:,Tfijo,rfijo);

drho_dp = rhoP(2:end)-rhoP(1:end-1); % [kg m-3/ 10 hPa] 

% cambio a unidad 10'5 Pa equivalente a 1 bar (para comparar con el
% oceano)
drho_dp = drho_dp*100; % [kg m-3 / bar]

% calculo coeficiente de compresibilidad en la atmosfera... dividimos por
% un valor promedio de la densidad en superficie

gama_atm = mean(drho_dp)/1.29; % en bar-1 % dejar este valor anotado para la comparacion con el oceano

%1.29 densidad promedio en superficie
%gama_atm=0,98

%% Parte 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varia gama_atm con T? segun yo si

% rho(p,T,r)
rfijo = find(r == 5); % posicion de r = 5 g/kg
rhoPT = squeeze(rho(:,:,rfijo)); %Rho depends pressure nd temperature

figure
contourf(p,T,rhoPT')
xlabel('Presión [hPa]')
ylabel('Temperatura [°C]')
 % las isopicnas no son verticales... mayor dependencia de ambas variables,
 % tanto presion como temperatura, cerca de superficie son mas diagonales,
 % luego en superficie es mas dependiente de ambas variables

%La temperatura sin fijar
for i = 1:length(T)
    drho_dp(:,i) = rhoPT(2:end,i)-rhoPT(1:end-1,i); % [kg m-3/ 10 hPa] 
end

% cambio a unidad 10'5 Pa equivalente a 1 bar (para comparar con el
% oceano)
drho_dp = drho_dp*100; % [kg m-3 / bar]

figure()
plot(p(2:end),drho_dp) % es un valor constante de 1,27 kg m-3 /bar

% calculo coeficiente de compresibilidad en la atmosfera... dividimos por
% un valor promedio de la densidad en superficie

gama_atm = drho_dp/1.29; 

figure()
contourf(p(2:end),T,gama_atm'), colorbar % es dependiente de T?
xlabel('Presión [hPa]')
ylabel('Temperatura [°C]')
title('Coeficiente de compresibilidad de la atmósfera con T cambiante')

% es dependiente de T? al estar horizontales totalmente las lineas dependen
% unicamente de la Temperatura mas que de la presion

%% Parte 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcular el coeficiente de compresibilidad en el oceano y comparar con la
% atmosfera

% Ecuación de estado del agua de mar

% Variables de entrada:
% T, temperatura en C
% p, presión en dbar
% S, Salinidad en psu

clear all

p = 0:1:400; %entre 0 y 400 bar (que equivale a 4000 m de profundidad)
S = 30:1:40;
T = -2:1:30;

for h = 1:length(p)
    for i = 1:length(T)
        for j = 1:length(S)

            % calculo de densidad(S,T,0)
            A = 999.842594+6.793952E-2*T(i)-9.095290E-3*T(i)^2+1.001685E-4*T(i)^3-...
                1.120083E-6*T(i)^4+6.536332E-9*T(i)^5;
            B = 8.24493E-1-4.0899E-3*T(i)+7.6438E-5*T(i)^2-8.2467E-7*T(i)^3+...
                5.3875E-9*T(i)^4;
            C = -5.72466E-3+1.0227E-4*T(i)-1.6546E-6*T(i)^2;
            D = 4.8314E-4;

            rhoST(h,i,j) = A+B*S(j)+C*S(j)^(3/2)+D*S(j)^2;

            % cálculo del módulo de volumen secante K(S,T,p)
            E = 19652.21+148.4206*T(i)-2.327105*T(i)^2+1.360477E-2*T(i)^3-5.155288E-5*T(i)^4;
            F = 54.6746-0.603459*T(i)+1.09987E-2*T(i)^2-6.1670E-5*T(i)^3;
            G = 7.944E-2+1.6483E-2*T(i)-5.3009E-4*T(i)^2;
            H = 3.239908+1.43713E-3*T(i)+1.16092E-4*T(i)^2-5.77905E-7*T(i)^3;
            I = 2.2838E-3-1.0981E-5*T(i)-1.6078E-6*T(i)^2;
            J = 1.91075E-4;
            M = 8.50935E-5-6.12293E-6*T(i)+5.2787E-8*T(i)^2;
            N = -9.9348E-7+2.0816E-8*T(i)+9.1697E-10*T(i)^2;
            % Cálculo de la densidad(S,T,p)

            K = E+F*S(j)+G*S(j)^(3/2)+(H+I*S(j)+J*S(j)^(3/2))*(p(h)/10)+...
                (M+N*S(j))*(p(h)/10)^2;

            rhoW(h,i,j) = rhoST(h,i,j)/(1-(p(h)/10)/K);
        end
    end
end

% calculo de coeficientes que determinan la variacion de la densidad (rho)
% a los cambios de las variables T, p y r
% rho(p,T,r)

Tfijo = find(T == 0); % posicion de p = 0 °C
Sfijo = find(S == 35); % posicion de S = 30 psu o partes por millon

rhoP = squeeze(rhoW(:,Tfijo,Sfijo));

drho_dp = rhoP(2:end)-rhoP(1:end-1); % [kg m-3/ bar] 

figure()
plot(p(2:end),drho_dp) % valor casi constante de 0,0047 kg m-3 /bar

% calculo coeficiente de compresibilidad en el oceano... dividimos por
% un valor promedio de la densidad
gama_oce = mean(drho_dp)/1027; % en bar-1 % comparar con la atmosfera

%Densidad promedio del oceano es de 1027 kg /m3
%gama_oce=4,61x10^-6

%El de atm es mucho mayor que el de oceano, lo que quiere decir que atm es
%mucho mas compresibleque el oceano, ademas demuestra que el agua es un
%fluido incompresible


%% Parte 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo de la concentración salina y dependencia con T, a p = 0

pfijo = find(p == 0); % posicion de p = 0 bar
rhoTS = squeeze(rhoW(pfijo,:,:)); 


figure()
contourf(T,S,rhoTS'), colorbar % es dependiente de S?
xlabel('Temperatura [°C]')
ylabel('Salinidad [psu]')
%isolineas curvas y diagonales,lo que demuestra que la densidad varia por 
% salinidad y temperatura, igualmente puede haber otro factor que influya

for i = 1:length(T)
    drho_dS(i,:) = rhoTS(i,2:end)-rhoTS(i,1:end-1); % [kg m-3/ psu] 
end
contrac_salina = drho_dS/1027; %llamado beta

figure()
contourf(T,S(2:end),contrac_salina'), colorbar % es dependiente de S?
xlabel('Temperatura [°C]')
ylabel('Salinidad [psu]')
%Lineas verticales, se van haciendo horizontales a medida que aumenta la
%temperaatura, esto debido a que densidad depende mas de la propiedad del
%agua en temperaturas bajas, mientras que a altas temperaturas, ambas
%variables, temp y salinidad ponen de su parte

posS = find(S(2:end) == 35);
temp = [-2 0 5 10 15 20 25 30];
for i = 1:length(temp)
    posT(i) = find(T == temp(i));
end

contrac_salina(posT,posS); 
% son valores bastante aproximados a la Tabla en la pagina 10 de la clase 2
