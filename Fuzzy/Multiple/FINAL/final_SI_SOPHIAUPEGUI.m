clc
% clear
close all

%% Datos

% Cargar datos desde un archivo Excel
filename = 'datos (1).xlsx';
loadedData = xlsread(filename);

% Definicion / asignacion de cada columna
tiempo = loadedData(:, 1)'; 
salida = loadedData(:, 3)'; 
entrada = loadedData(:, 2)'; 

figure(1)
plot(tiempo, salida, 'b')
hold on
plot(tiempo, entrada, 'r', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
hold off

% Separar indices
i1 = 301; 
i1end = 600;
i2 = i1end+1;
i2end = 900;
i3 = i2end+1;
i3end = 1200;
i4 = i3end+1;
i4end = 1500;

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

% Modelo segun system identification
Tau1 = 48.009;
K1 = 1.32;
Theta1 = 0.929; 

% Modelo
Num1 = [K1];
Den1 = [Tau1 1];
Gs1 = tf(Num1, Den1, 'inputdelay', Theta1) 

[respuesta_modelo1, ~] = lsim(Gs1, escalon1, tiempo1);

figure(2)
plot(tiempo1, escalon1, 'r','LineWidth', 1.2)
hold on
plot(tiempo1, salida1, 'b','LineWidth', 1.3)
plot(tiempo1, respuesta_modelo1, 'g', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
y_inf = salida1(1);
y_sup = salida1(end) + 3;
ylim([y_inf y_sup])
xlim([tiempo1(1) tiempo1(end)])
hold off

% Criterio de Tau equivalente
Tau_eq = Tau1 + Theta1/2;
Gs_sint = tf(Num1, Den1);
LC = feedback(Gs_sint, 1);
[numLC, denLC] = tfdata(LC, 'v');
if denLC(2) ~= 1
    Teq = denLC(1)/denLC(2);
else
    Teq = denLC(1);
end

Ts1_r = [0.2*(Teq+Theta1) 0.6*(Teq+Theta1)]; % rango
Ts1_p = (0.2*(Teq+Theta1) + 0.6*(Teq+Theta1))/2; % ponderado


%% Modelo: escalon 2

% Modelo segun system identification
Tau2 = 40.8427;
K2 = 1.05;
Theta2 = 0; 

% Modelo
Num2 = [K2];
Den2 = [Tau2 1];
Gs2 = tf(Num2, Den2, 'inputdelay', Theta2) 

u = (escalon2(end)-escalon1(end)) * ones(size(tiempo2));
[respuesta_modelo2, ~] = lsim(Gs2, u, tiempo2);
respuesta_modelo2 = respuesta_modelo2 + salida2(1);

figure(3)
plot(tiempo2, escalon2, 'r','LineWidth', 1.2)
hold on
plot(tiempo2, salida2, 'b','LineWidth', 1.3)
plot(tiempo2, respuesta_modelo2,'g', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
y_inf = salida2(1) - 3;
y_sup = salida2(end) + 3;
ylim([y_inf y_sup])
xlim([tiempo2(1) tiempo2(end)])
hold off

% Criterio de Tau equivalente
Tau_eq = Tau2 + Theta2/2;
Gs_sint = tf(Num2, Den2);
LC = feedback(Gs_sint, 1);
[numLC, denLC] = tfdata(LC, 'v');
if denLC(2) ~= 1
    Teq = denLC(1)/denLC(2);
else
    Teq = denLC(1);
end

Ts2_r = [0.2*(Teq+Theta2) 0.6*(Teq+Theta2)]; % rango
Ts2_p = (0.2*(Teq+Theta2) + 0.6*(Teq+Theta2))/2; % ponderado


%% Modelo: escalon 3

% Modelo segun system identification
Tau3 = 50.6613;
K3 = 1.44;
Theta3 = 0;

% Modelo
Num3 = [K3];
Den3 = [Tau3 1];
Gs3 = tf(Num3, Den3, 'inputdelay', Theta3) 

u = (escalon3(end) - escalon2(end)) * ones(size(tiempo3));
[respuesta_modelo3, ~] = lsim(Gs3, u, tiempo3);
respuesta_modelo3 = respuesta_modelo3 + salida3(1);

figure(4)
plot(tiempo3, escalon3, 'r','LineWidth', 1.2)
hold on
plot(tiempo3, salida3, 'b','LineWidth', 1.3)
plot(tiempo3, respuesta_modelo3,'g', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
y_inf = escalon3(1) - 3;
y_sup = salida3(end) + 3;
ylim([y_inf y_sup])
xlim([tiempo3(1) tiempo3(end)])
hold off

% Criterio de Tau equivalente
Tau_eq = Tau3 + Theta3/2;
Gs_sint = tf(Num3, Den3);
LC = feedback(Gs_sint, 1);
[numLC, denLC] = tfdata(LC, 'v');
if denLC(2) ~= 1
    Teq = denLC(1)/denLC(2);
else
    Teq = denLC(1);
end

Ts3_r = [0.2*(Teq+Theta3) 0.6*(Teq+Theta3)]; % rango
Ts3_p = (0.2*(Teq+Theta3) + 0.6*(Teq+Theta3))/2; % ponderado


%% Modelo: escalon 4

% Modelo segun system identification
Tau4 = 49.2614;
K4 = 1.01;
Theta4 = 0;

% Modelo
Num4 = [K4];
Den4 = [Tau4 1];
Gs4 = tf(Num4, Den4, 'inputdelay', Theta4) 

u = (escalon4(end) - escalon3(end)) * ones(size(tiempo4));
[respuesta_modelo4, ~] = lsim(Gs4, u, tiempo4);
respuesta_modelo4 = respuesta_modelo4 + salida4(1);

figure(4)
plot(tiempo4, escalon4, 'r','LineWidth', 1.2)
hold on
plot(tiempo4, salida4, 'b','LineWidth', 1.3)
plot(tiempo4, respuesta_modelo4,'g', 'LineWidth', 1.2)
xlabel('Tiempo (s)')
ylabel('Salida')
y_inf = escalon4(1) - 3;
y_sup = salida4(end) + 3;
ylim([y_inf y_sup])
xlim([tiempo4(1) tiempo4(end)])
hold off

% Criterio de Tau equivalente
Tau_eq = Tau4 + Theta4/2;
Gs_sint = tf(Num4, Den4);
LC = feedback(Gs_sint, 1);
[numLC, denLC] = tfdata(LC, 'v');
if denLC(2) ~= 1
    Teq = denLC(1)/denLC(2);
else
    Teq = denLC(1);
end

Ts4_r = [0.2*(Teq+Theta4) 0.6*(Teq+Theta4)]; % rango
Ts4_p = (0.2*(Teq+Theta4) + 0.6*(Teq+Theta4))/2; % ponderado


%% Discretizacion

out = [Ts1_r; Ts2_r; Ts3_r; Ts4_r]
valor_min = max(out(:, 1)); % el valor más grande de los mínimos
valor_max = min(out(:, 2)); % el valor más pequeño de los máximos

% Verifica si hay intersección
if valor_min <= valor_max
    fprintf('El rango de intersección es: [%.4f, %.4f]\n', valor_min, valor_max);
else
    disp('No hay intersección entre los rangos.');
end

Ts = round((valor_min + valor_max)/2 * 10) / 10

% Discretizar modelo
Gz1 = c2d(Gs1, Ts, 'zoh')
figure()
step(Gs1, Gz1)

% Discretizar modelo
Gz2 = c2d(Gs2, Ts, 'zoh')
figure()
step(Gs2, Gz2)

% Discretizar modelo
Gz3 = c2d(Gs3, Ts, 'zoh')
figure()
step(Gs3, Gz3)

% Discretizar modelo
Gz4 = c2d(Gs4, Ts, 'zoh')
figure()
step(Gs4, Gz4)

