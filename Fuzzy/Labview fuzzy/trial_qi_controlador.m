clc 
%clear
close all

%% Definición de la planta discreta
Tau = 47.5;
K = 3;
Theta = 2.1;
Ts = 1.05;

Gs = tf([K], [Tau 1], 'InputDelay', Theta);  % Planta continua
Gz = c2d(Gs, Ts, 'zoh');                     % Planta discreta

%% Controlador PID discreto
% Cambia estos valores para probar otros controladores

%error ss
% q0 = 0.8;
% q1 = -1.4;
% q2 = 0.6;

%overshoot 33%
% q0 = 1.5;
% q1 = -2;
% q2 = 0.6;
% 
% % q0 = 1.4;
% % q1 = -1.8;
% % q2 = 0.5;

% q0 = 14;
% q1 = -20;
% q2 = 7;

% Controlador C(z) = (q0 + q1*z^-1 + q2*z^-2) / (1 - z^-1)
NumC = [q0 q1 q2];
DenC = [1 -1 0];
C = tf(NumC, DenC, Ts);

% Sistema en lazo cerrado
GC = feedback(C * Gz, 1);

% Gráfica de la respuesta al escalón
figure
step(GC)
title('Respuesta del sistema con controlador discreto PID')
xlabel('Tiempo [s]')
ylabel('Salida')
grid on


Mp = stepinfo(GC).Overshoot