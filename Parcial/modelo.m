%% PARCIAL - SISTEMAS INTELIGENTES 2025-1
% SOPHIA UPEGUI ROBLEDO

clc
close all
clear

%% Modelo del sistema

K = 1.25;
Tau = 25;
Theta = 0.5; 

Num = [K];
Den = [Tau 1];
Gs = tf(Num, Den, 'inputdelay', Theta); % lazo abierto

k =1;
Gs_lc = feedback(Gs,k)

%% Tiempo de muestreo

% Criterio de tau equivalente
w = sqrt(3) % del ultimo termino del polinomio denominador en lazo cerrado
xi = 1/(2*w)
tau_eq = 1/(xi*w); % si es de primer 
Ts1 = 0.2*tau_eq;

sysd1 = c2d(Gs,Ts1)
step(Gs, sysd1)


% Criterio de tiempo de establecimiento
ts = stepinfo(Gs_lc).SettlingTime; % tiemplo de establecimiento del step de la ft_lc
 % el valor elegido de 0.05*ts<Ts<0.15*ts
lb = 0.05*ts;
ub = 0.15*ts;
lim = [lb up]
Ts3 = 2.5;

Gz = c2d(Gs, Ts3, 'zoh')
step(Gs,Gz)
