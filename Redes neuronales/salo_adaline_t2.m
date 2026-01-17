clc;
clear;
close all;

%% Cargar la base de datos
data = load('Datos_2025.txt'); % Cargar datos
time = data(:,1); % Tiempo
u_o = data(:,3);  % Entrada
y_o = data(:,2);  % Salida (Temperatura)

% Regresores
regre = 2;
num_samples = length(u_o) - regre;  
x = zeros(num_samples, 4); 
T_seleccionados = zeros(num_samples, 1);

% Construcción de regresores
for k = 1:num_samples
    x_o(k, :) = [y_o(k+1), y_o(k), u_o(k+2), u_o(k+1)];
    T_seleccionados(k) = y_o(k+2);
end

% Mezclar aleatoriamente los datos
total_data = size(x_o, 1);
indices = randperm(total_data); % Permutación aleatoria de índices

% Aplicar la mezcla
x = x_o(indices, :);
T_seleccionados = T_seleccionados(indices);

% División de los datos
train_size = round(0.7 * total_data);
val_size = round(0.15 * total_data);
test_size = total_data - train_size - val_size; % 15% restante para prueba

% Asignar subconjuntos
x_train = x(1:train_size, :);
T_train = T_seleccionados(1:train_size);

x_val = x(train_size+1:train_size+val_size, :);
T_val = T_seleccionados(train_size+1:train_size+val_size);

x_test = x(train_size+val_size+1:end, :);
T_test = T_seleccionados(train_size+val_size+1:end);

% Inicialización de pesos y bias
w = rand(1,4);
b = rand(1);

% Tasa de aprendizaje
n = 0.0000001;

%% Entrenamiento

% Error permitido con base al error cuadrático medio
e_perm = 0.001;

% Inicializar error cuadrático medio
EMC = 100; % Un valor inicial grande
EMC_G = [];

% Inicializo el número de épocas
iter = 0;

while ((EMC > e_perm) && (iter < 1000)) % Sale por error permitido o # de iteraciones
    iter = iter + 1; % Aumentar iteración

    % Entrenamiento
    for k = 1:train_size  % Iterar sobre todas las muestras
        % Evaluar salida Yt
        YT(k) = w * x_train(k, :)' + b; % w (1x4) * x(k,:); (4x1) = escalar

        % Calcular el error e
        e(k) = T_train(k)' - YT(k);

        % Actualizar pesos y bias
        w = w + n * (e(k) * x_train(k, :)) ;
        b = b + n * e(k);
    end 

    % Calcular el error cuadrático medio
    EMC = (1/train_size) * sum(e.^2);
    EMC_G(iter) = EMC;
    
    % Validación
    for k = 1:val_size  

        YT_v(k) = w * x_val(k, :)' + b; 

        % Calcular el error e
        e_v(k) = T_val(k)' - YT_v(k);
    end 

    % Calcular el error cuadrático medio de validación
    EMC_V = (1/val_size) * sum(e_v.^2);
    EMC_G_V(iter) = EMC_V;
end

% Prueba
for k = 1:total_data  

    YT_prueba(k) = w * x_o(k, :)' + b; 
end 

% Mostrar resultados
disp('Iteraciones:')
disp(iter)
disp('Pesos finales:')
disp(w)
disp('Bias')
disp(b)
figure;
semilogx(EMC_G)
title('Error cuadrático medio entrenamiento')
xlabel('Iteraciones')
ylabel('EMC')
grid on;

figure;
semilogx(EMC_G_V)
title('Error cuadrático medio validación')
xlabel('Iteraciones')
ylabel('EMC')
grid on;

figure;
plot(y_o, 'b');
hold on;
plot(YT_prueba, 'r');
hold on;
plot(u_o, 'k')
legend('Real', 'Predicción', 'Entrada');
