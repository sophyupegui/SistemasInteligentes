clc
clear
close all

%% Modelo de la planta 

Tau = 32.5;
K = 1.4951;
Theta = 0.4; 

Gs = tf([K], [Tau 1], 'InputDelay', Theta);  % Planta continua
Ts = 0.2;                                   % Tiempo de muestreo
Gz = c2d(Gs, Ts, 'zoh');                     % Planta discreta


%% Ziegler-Nichols - Curva de Reacción

% PID continuo
Kp1 = (1.2 * Tau) / (K * Theta);
Tau_i1 = 2 * Theta;
Tau_d1 = 0.5 * Theta;
Ki1 = Kp1 / Tau_i1;
Kd1 = Kp1 * Tau_d1;

% Controlador continuo PID
C1 = tf([Kd1 Kp1 Ki1], [1 0]);
C1z = c2d(C1, Ts, 'tustin');  % Discretización por Tustin

% Obtener coeficientes q0, q1, q2
[num1, ~] = tfdata(C1z, 'v');
q0_CR = num1(1);
q1_CR = num1(2);
q2_CR = num1(3);

% Mostrar resultados
disp('=== Ziegler-Nichols: Curva de Reacción ===');
fprintf('Kp = %.4f | Ki = %.4f | Kd = %.4f\n', Kp1, Ki1, Kd1);
fprintf('q0 = %.6f\n', q0_CR);
fprintf('q1 = %.6f\n', q1_CR);
fprintf('q2 = %.6f\n', q2_CR);

% Lazo cerrado
GCR = feedback(C1 * Gs, 1);

% Gráfica
figure(1)
step(GCR)
title('ZN Curva de Reacción - Respuesta al escalón')
xlabel('Tiempo'), ylabel('Salida')
grid on

%% Ziegler-Nichols - Ganancia Límite

% Obtener Ganancia límite y frecuencia límite
[Klim, ~, Wlim, ~] = margin(Gs);

if isinf(Klim)
    error('Klim = Inf: no se puede aplicar Ziegler-Nichols por ganancia límite');
end

Ku = Klim;
Tu = 2 * pi / Wlim;

% Ganancias PID
Kp2 = 0.6 * Ku;
Tau_i2 = 0.5 * Tu;
Tau_d2 = 0.125 * Tu;
Ki2 = Kp2 / Tau_i2;
Kd2 = Kp2 * Tau_d2;

% Coeficientes q0, q1, q2 (discretización directa)
q0_GL = Kp2 * (1 + Ts/(2*Tau_i2) + Tau_d2/Ts);
q1_GL = Kp2 * (Ts/(2*Tau_i2) - 2*Tau_d2/Ts - 1);
q2_GL = Kp2 * Tau_d2 / Ts;

% Mostrar resultados
disp('=== Ziegler-Nichols: Ganancia Límite ===');
fprintf('Ku = %.4f | Tu = %.4f\n', Ku, Tu);
fprintf('Kp = %.4f | Ki = %.4f | Kd = %.4f\n', Kp2, Ki2, Kd2);
fprintf('q0 = %.6f\n', q0_GL);
fprintf('q1 = %.6f\n', q1_GL);
fprintf('q2 = %.6f\n', q2_GL);

% Crear controlador discreto directo
C2z = tf([q0_GL q1_GL q2_GL], [1 -1 0], Ts);

% Lazo cerrado con planta discreta
GCGL = feedback(C2z * Gz, 1);

% Gráfica
figure(2)
step(GCGL)
title('ZN Ganancia Límite - Respuesta al escalón')
xlabel('Tiempo'), ylabel('Salida')
grid on
