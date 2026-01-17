%% Adaline u(k-1) con data base

clc
clear
close all

% Enunciado: 

% Buenos días estudiantes, comparto la base de datos que vamos a usar para modelar un sistema dinámico (planta de temperatura). 
% Por ahora deben:
% Cargar la base de datos
% Graficar la base datos 
% Identificar las columnas (salida (Temperatura), entrada (u))
% Preparar 70% de los datos para entrenamiento
% Preparar 15% de los datos para Validación  
% Preparar 15% de los datos para Prueba 
% Programa la estrategia para generar los regresores (y(-1),y(k-2).....u(k),u(k-1)....).


%% Cargar base de datos 
loadedData = load('Datos_2025.txt');

% Definicion / asignacion de cada columna
tiempo = loadedData(:, 1)'; % 1st Columna (tiempo t) 
salida = loadedData(:, 2)'; % 2nd Columna (entrada u(k)) 
entrada = loadedData(:, 3)'; % 3rd Columna (tiempo t)

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

%% Cantidad de regresores para estimación

% Determino cantidad de regresores que voy a usar
cant_regresores = 2;

% Regresores de la salida (y)
y_k_n = zeros(cant_regresores, length(entrada)); %y_k-1 y y_k-2

% Regresores de la entrada (u)
u_k_n = zeros(cant_regresores, length(entrada)); %u_k y u_k-1

% Construcción de la matriz de regresores
for i = cant_regresores+1:length(entrada)
    
    % Regresores de salida
    for k = 1:cant_regresores+1
        y_k_n(k,i) = salida(i-(k-1));
    end

    % Regresores de entrada % u_k (entrada actual)
    for k = 1:cant_regresores 
        u_k_n(k, i) = entrada(i-(k-1)); 
    end

end

y_selec = [y_k_n(3,:)];
x_reg = [y_k_n(1,:)' y_k_n(2,:)' u_k_n']

%% Aleatorio dato

% Total
num_total_datos = length(x_reg);

% 70% para el entrenamiento
num_datos_entrenamiento = round(0.7 * num_total_datos);

% 15% validación
num_datos_validacion = round(0.15 * num_total_datos);

% Generar índices aleatorios para seleccionar datos 
indices_entren = randperm(num_total_datos, num_datos_entrenamiento);
indices_val = randperm(num_total_datos, num_datos_validacion);
indices_pruebas = randperm(num_total_datos);

% Seleccionar los datos aleatorios entrenamiento
x_entren = x_reg(indices_entren,:); % Patron de entrenamiento
T_entren = y_selec(indices_entren);

% Seleccionar los datos aleatorios validacion
x_val = x_reg(indices_val,:);
T_val = y_selec(indices_val);

%% Entrenamiento & validacion

% Inicialización de pesos
w = rand(1,4);

% Inicialización del sesgo
b = rand(1);

% Taza de aprendizaje
n = 0.000001;

% Error permitido con base al error cuadratico medio
e_perm = 0.001;

% Error cuadratico medio
EMC_entren = 100; %un valor cualquiera con tal de que sea muy grande
EMC_G_entren = [];

EMC_val = 100; %un valor cualquiera con tal de que sea muy grande
EMC_G_val = [];

% Inicializo # iteraciones (# iter = cantidad de epocas pasadas)
iter = 0;

while ((EMC_entren > e_perm) && (iter<2000)) % sale por error permitido y # de iteraciones

    for k = 1:length(x_entren)

        YT_entren(k) = w*x_entren(k,:)' + b; % w (1x4) * x(k,:); (4x1) = escalar
        
        % Calcular el error e
        e_entren(k) = T_entren(k) - YT_entren(k); 

        % Acualizamos pesos y bias
        w = w + n*(e_entren(k)*x_entren(k,:));
        b = b + n*e_entren(k);

    end 

    iter = iter+1; % aumento iteracion

    % Error cuadratico medio de entrenamiento
    EMC_entren = (1/length(x_entren))*sum(e_entren.^2); %pongo el . (tras el e) para que lo haga elemento a elemento
    EMC_G_entren(iter) = EMC_entren; % i think this is the same as: EMC_epoca = [iter, EMC]  


    % Validacion
    for k = 1:length(x_val)

        YT_val(k) = w*x_val(k,:)' + b; 
        
        % Calcular el error e
        e_val(k) = T_val(k) - YT_val(k); 

    end
    
    % Error cuadratico medio de validacion
    EMC_val = (1/length(x_val))*sum(e_val.^2); %pongo el . (tras el e) para que lo haga elemento a elemento
    EMC_G_val(iter) = EMC_val; % i think this is the same as: EMC_epoca = [iter, EMC]  


end

disp('Iteraciones:')
disp(iter)

disp('Pesos finales:')
disp(w)

disp('Bias')
disp(b)

disp('Deseados vs Entrenamiento:')
prueba_entren = [T_entren' YT_entren'];
disp(prueba_entren)

disp('Deseados vs Validacion:')
prueba_val = [T_val' YT_val'];
disp(prueba_val)


figure(2)
hold on
semilogx(EMC_G_entren)
grid on
title("EMC Entrenamiento")
xlabel("Iteraciones")
ylabel("EMC")
hold off

figure(3)
hold on
semilogx(EMC_G_val)
grid on
title("EMC Validación")
xlabel("Iteraciones")
ylabel("EMC")
hold off


%% Prueba (100% de datos)

for k = 1:num_total_datos  
    YT_prueba(k) = w * x_reg(k, :)' + b; 
end 


figure(4)
hold on
plot(salida, 'b')
plot(entrada, 'r')
plot(YT_prueba, 'g')
legend('Salida', 'Entrada', 'Prediccion');
hold off

