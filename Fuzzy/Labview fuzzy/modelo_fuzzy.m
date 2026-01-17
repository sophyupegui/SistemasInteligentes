%clear 
close all
clc

%% Cargar base de datos 

loadedData = load('Datos6.txt');

% Definicion / asignacion de cada columna
tiempo = loadedData(:, 1)'; % 1st Columna (tiempo t) 
salida = loadedData(:, 2)'; % 2nd Columna (salida y(k)) 
entrada = loadedData(:, 3)'; % 3rd Columna (entrada u(k))

%% Graficar base de datos

figure(1)
hold on
plot(tiempo, entrada, Color='r', LineWidth=2)
grid on
plot(tiempo, salida, Color='b', LineWidth=1.3)
legend('Entrada','Salida','Location','southeast')
title(['Base de datos'])
xlabel("Tiempo")
ylabel("Temperatura")
hold off

% Step seleccionado
t_inicio = 540;
t_final = 890;

entrada_sel = entrada(t_inicio:t_final);
salida_sel = salida(t_inicio:t_final);
tiempo_sel = tiempo(t_inicio:t_final);

figure(2)
hold on
plot(tiempo_sel, entrada_sel, Color='r', LineWidth=2)
grid on
plot(tiempo_sel, salida_sel, Color='b', LineWidth=1.3)
legend('Entrada','Salida','Location','southeast')
title(['Step seleccionado'])
xlabel("Tiempo")
ylabel("Temperatura")
hold off


%% Planta

% Parametros del modelo
Tau = 47.5;
K = 3;
Theta = 2.1; 

% Modelo
Num = [K];
Den = [Tau 1];
Gs = tf(Num, Den, 'inputdelay', Theta) 

% Criterio de Tau equivalente
Tau_eq = Tau + Theta/2;
Ts1 = Tau_eq / 10;

% Criterio de tiempo de establecimiento
Gs_lc = feedback(Gs,1)
ts = stepinfo(Gs_lc).SettlingTime; % Step de lazo cerrado 
rango = [0.05*ts 0.15*ts]
Ts2 = 0.05*ts; % el valor elegido de 0.05*ts<Ts<0.15*ts

% Definir tiempo de muestreo
Ts_ponderado = min([Ts1 Ts2]);

% Comparación con retardo
if Ts_ponderado < Theta/2
    Ts = Ts_ponderado;
else
    Ts = Theta/2;
end

% Discretizar modelo
Gz = c2d(Gs, Ts, 'zoh')

% Graficar
step(Gs,Gz)
legend('Modelo continuo','Modelo discreto',location='southeast')


[out_Gs, ty_Gs] = step(Gs);
[out_Gz, ty_Gz] = step(Gz);
plot(ty_Gs, squeeze(out_Gs),LineWidth=2,Color='b')
hold on
plot(ty_Gz, squeeze(out_Gz),LineWidth=1.3,Color='r')
legend('Modelo continuo','Modelo discreto',location='southeast')
title('Respuesta al escalón de Gs vs Gz')
grid on
hold off

step(Gz)


%% Comparacion

U = entrada_sel(200) - entrada_sel(1);
offset = salida_sel(1)

Y = Gz*U + offset;
[out, ty] = step(Y);

figure(3)
hold on
plot(tiempo_sel, entrada_sel, Color='r', LineWidth=2)
plot(tiempo_sel, salida_sel, Color='b', LineWidth=1.3)
plot(ty, squeeze(out), LineWidth=2,Color='g')
grid on
legend('Entrada', 'Datos', 'Modelo', 'Location','southeast')
title('Validación de la aproximación POR del nivel')
xlabel("Tiempo")
ylabel("Temperatura")
xlim([200 700])
hold off


%% Comparacion
% Calcular el cambio en la entrada
deltaU = entrada_sel(200) - entrada_sel(1); 
offset = salida_sel(1);

% Crear vector de tiempo para simulación
t_sim = tiempo_sel - tiempo_sel(1);

% Simular la respuesta del modelo
[y_model, t_model] = step(Gs * deltaU, t_sim(end));

% Ajustar el tiempo del modelo para que coincida con los datos
t_model = t_model + tiempo_sel(1);

% Añadir el offset a la respuesta del modelo
y_model = y_model + offset;

% Graficar comparación
figure(3)
hold on
plot(tiempo_sel, entrada_sel, 'r', 'LineWidth', 2)
plot(tiempo_sel, salida_sel, 'b', 'LineWidth', 1.3)
plot(t_model, y_model, 'g', 'LineWidth', 2)
grid on
legend('Entrada', 'Datos', 'Modelo', 'Location','southeast')
title('Validación de la aproximación PC de la temperatura')
xlabel("Tiempo")
ylabel("Temperatura")
hold off

