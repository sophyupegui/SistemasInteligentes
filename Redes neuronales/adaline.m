%% Adaline

clc
clear
close all

%% Defino parametros

% Entradas
x1 = [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];
x2 = [0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
x3 = [0 1 1 0 0 1 1 0 0 1 1 0 0 1 1];
x4 = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
x = [x1; x2; x3; x4];

% Target outputs
T = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];

%% Inicialización de parametros

% Inicializo pesos
w = rand(1, 4); %(1, cantidad de entradas)

% Inicializo sesgo
b = rand(1);

% Taza de aprendizaje
n = 0.001;

% Error permitido con base al error cuadratico medio
e_perm = 0.001;

% Error cuadratico medio
EMC = 100; %un valor cualquiera con tal de que sea muy grande
EMC_G = [];

%% Seleccionar muestra de datos

% Número total de datos
num_total_datos = length(x1);

% 70% del total de datos
num_datos_seleccionados = round(0.7 * num_total_datos);

% Seleccion de datos
x_seleccionados = x(:,1:num_datos_seleccionados); % Patron de entrenamiento
T_seleccionados = T(1:num_datos_seleccionados); % Asegúrate de seleccionar también las salidas correspondientes

% Tamaño final
[f, c] = size(x_seleccionados);

%% Entrenamiento

% Inicializo # iteraciones (# iter = cantidad de epocas pasadas)
iter = 0;

while ((EMC > e_perm) && (iter<2000)) % sale por error permitido y # de iteraciones

    for k = 1:c

        % Evaluar salida Yt (por ser activador de recta Yp y Yt son lineales)
        YT(k) = w*x_seleccionados(:,k) + b; 
        
        % Calcular el error e
        e(k) = T_seleccionados(k) - YT(k); 

        % Acualizamos pesos y bias
        Delta_w = w + n*(e(k)*x_seleccionados(:,k)');
        Delta_b = b + n*e(k);
        w = Delta_w;
        b = Delta_b;

    end 
      
    iter = iter+1; % aumento iteracion

    % Error cuadratico medio
    EMC = (1/c)*sum(e.^2); %pongo el . (tras el e) para que lo haga elemento a elemento
    EMC_G(iter) = EMC; % i think this is the same as: EMC_epoca = [iter, EMC]  

end

w
YT
semilogx(EMC_G)

