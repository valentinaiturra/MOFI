clc
clear all
close all

M=[];
%Masas moleculares de los componentes
M(1) = 28.013; % N2
M(2) = 32;     % O2
M(3) = 39.95;  % Ar
M(4) = 44.01;  % CO2
M(5) = 20.18;  % Ne
M(6) = 4;      % He
M(7) = 16.04;  % CH4
M(8) = 83.8;   % Kr
M(9) = 2.02;   % H2

%masa de los componentes

m(1) = 75.52*10;   % N2
m(2) = 23.14*10;   % O2
m(3) = 1.29*10;    % Ar
m(4) = 0.063*10;   % CO2
m(5) = 0.0013*10;  % Ne
m(6) = 0.0007*10;  % He
m(7) = 0.001*10;   % CH4
m(8) = 0.00029*10;  % Kr
m(9) = 0.000003*10; % H2

% Sumamos las masas
sum_m = cumsum(m);
%Hacemos la suma sabiendo que moles n = m./M
sum_n = cumsum(m./M);
% Masa molecular del aire seco
Md = sum_m./sum_n; 

%A partir de usar cumsum, la suma acumulativa a hecho que podamos ver en
%cada caso si se consideran menos de nueve constituyentes de la atmosfera.

%Ahora determinamos el error, coinsideramos el caso de 9 constituyentes, es
%decir Md(9) como el caso mas correcto
deltaMd = (Md-Md(end))/Md(end)*100;

% Graficamos una figura
figure(1)
plot(1:9,deltaMd,'o')
xlabel({'N° de moléculas consideradas', 'en la composición de la atmósfera'})
ylabel('Error de Md [%]')
grid on

%CONCLUSIONES
%a partir del cuarto valor, es decir tomando solo cuatro componentes
%atmosfericos el error disminuye bastante


%Podemos determinar el valor de la densidad, es decir a partir de la
%ecuacion de estado de la atmosfera


% Calculamos R para el aire seco considerando los diversos valores de Md
Rd = 8.3145./Md; 

% Escribimos los valores de la temperatura virtual Tv [K] y la presión [Pa]
% usando valores fijos
p = 1000*100;    % [Pa]
Tv = 25+273.15;  % [K]
rho = p./(Tv.*Rd);

% Veamos el porcentaje de error en rho
deltarho = (rho-rho(end))/rho(end)*100; % porcentaje de error en calculo de rho

