% clc
% clear
% close all

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

% Modelo segun system identification
Tau1 = 32.5;
K1 = 1.4951;
Theta1 = 0.4; 

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
feedback
Ts1 = Tau_eq / 10;

% Criterio de tiempo de establecimiento
Gs_lc = feedback(Gs1,1)
ts = stepinfo(Gs_lc).SettlingTime; % Step de lazo cerrado 
rango = [0.05*ts 0.15*ts];
Ts2 = 0.05*ts; % el valor elegido de 0.05*ts<Ts<0.15*ts

% Definir tiempo de muestreo
Ts_ponderado = min([Ts1 Ts2]);

% Comparación con retardo
if Ts_ponderado < Theta1/2
    Tse1 = Ts_ponderado;
else
    Tse1 = Theta1/2;
end


%% Modelo: escalon 2

% Modelo segun system identification
Tau2 = 48.051;
K2 = 1.7502;
Theta2 = 3.0605; 

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
Ts1 = Tau_eq / 10;

% Criterio de tiempo de establecimiento
Gs_lc = feedback(Gs2,1)
ts = stepinfo(Gs_lc).SettlingTime; % Step de lazo cerrado 
rango = [0.05*ts 0.15*ts];
Ts2 = 0.05*ts; % el valor elegido de 0.05*ts<Ts<0.15*ts

% Definir tiempo de muestreo
Ts_ponderado = min([Ts1 Ts2]);

% Comparación con retardo
if Ts_ponderado < Theta2/2
    Tse2 = Ts_ponderado;
else
    Tse2 = Theta2/2;
end


%% Modelo: escalon 3

% Modelo segun system identification
Tau3 = 45.5505;
K3 = 0.98;
Theta3 = 5.1795;

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
Ts1 = Tau_eq / 10;

% Criterio de tiempo de establecimiento
Gs_lc = feedback(Gs3,1)
ts = stepinfo(Gs_lc).SettlingTime; % Step de lazo cerrado 
rango = [0.05*ts 0.15*ts];
Ts2 = 0.05*ts; % el valor elegido de 0.05*ts<Ts<0.15*ts

% Definir tiempo de muestreo
Ts_ponderado = min([Ts1 Ts2]);

% Comparación con retardo
if Ts_ponderado < Theta3/2
    Tse3 = Ts_ponderado;
else
    Tse3 = Theta3/2;
end



%% Modelo: escalon 4

% Modelo segun system identification
Tau4 = 50;
K4 = 0.795;
Theta4 = 7;

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
Ts1 = Tau_eq / 10;

% Criterio de tiempo de establecimiento
Gs_lc = feedback(Gs4,1)
ts = stepinfo(Gs_lc).SettlingTime; % Step de lazo cerrado 
rango = [0.05*ts 0.15*ts];
Ts2 = 0.05*ts; % el valor elegido de 0.05*ts<Ts<0.15*ts

% Definir tiempo de muestreo
Ts_ponderado = min([Ts1 Ts2]);

% Comparación con retardo
if Ts_ponderado < Theta4/2
    Tse4 = Ts_ponderado;
else
    Tse4 = Theta4/2;
end

%% Discretizacion

%Ts = min([Tse1 Tse2 Tse3 Tse4])
Ts = 1

% Discretizar modelo
Gz1 = c2d(Gs1, Ts, 'zoh')
step(Gs1, Gz1)

% Discretizar modelo
Gz2 = c2d(Gs2, Ts, 'zoh')
step(Gs2, Gz2)

% Discretizar modelo
Gz3 = c2d(Gs3, Ts, 'zoh')
step(Gs3, Gz3)

% Discretizar modelo
Gz4 = c2d(Gs4, Ts, 'zoh')
step(Gs4, Gz4)

