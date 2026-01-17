clc
clear
close all

%% Datos

% Cargar datos desde un archivo Excel
loadedData = load('Datos6.txt');

% Definicion / asignacion de cada columna
tiempo = loadedData(:, 1)'; % 1st Columna (tiempo t) 
salida = loadedData(:, 2)'; % 2nd Columna (salida y(k)) 
entrada = loadedData(:, 3)'; % 3rd Columna (entrada u(k))

figure(1)
plot(tiempo, salida, 'b')
hold on
plot(tiempo, entrada, 'r', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
hold off

% Separar indices
i1 = 238;
i1end = 514;
i2 = i1end+1;
i2end = 877;
i3 = i2end+1;
i3end = 1207;
i4 = i3end+1;
i4end = 1589;

% Separo tiempo
tiempo1 = tiempo(i1:i1end);
tiempo2 = tiempo(i2:i2end);
tiempo3 = tiempo(i3:i3end);
tiempo4 = tiempo(i4:i4end);

% Separo escalones
escalon1 = entrada(i1:i1end);
escalon2 = entrada(i2:i2end);
escalon3 = entrada(i3:i3end);
escalon4 = entrada(i4:i4end);

% Salidas esperadas
salida1 = salida(i1:i1end);
salida2 = salida(i2:i2end);
salida3 = salida(i3:i3end);
salida4 = salida(i4:i4end);

%% Modelo: escalon 1

entrada_escalon = 10;
K = (salida1(end) - salida1(1)) / entrada_escalon;

% Estimar Tiempo muerto theta (L)
idx_inicio = find(salida1 > salida1(1) + 0.02*(salida1(end)-salida1(1)), 1);
L = tiempo1(idx_inicio) - tiempo1(1);

% Estimar Constante de Tiempo T
% Aproximación: tiempo donde la salida alcanza el 63.2% de su valor final
Y63 = salida1(1) + 0.632*(salida1(end) - salida1(1));
idx_T = find(salida1 >= Y63, 1);
T = tiempo1(idx_T) - tiempo1(idx_inicio);

sys = tf(K, [T 1], 'InputDelay', L);

entrada_1 = entrada_escalon * ones(size(tiempo1));
[respuesta_modelo, ~] = lsim(sys, entrada_1, tiempo1);

figure(2)
plot(tiempo1, escalon1, 'r')
hold on
plot(tiempo1, salida1, 'b','LineWidth', 1.6)
plot(tiempo1, respuesta_modelo, 'g', 'LineWidth', 1.2);
hold off
