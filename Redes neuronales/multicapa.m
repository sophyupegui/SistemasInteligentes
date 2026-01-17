%% Red neuronal multicapa

clc
close all
clear

%% Cargar base de datos 

loadedData = load('Datos_2025.txt');

% Definicion / asignacion de cada columna
tiempo = loadedData(:, 1)'; % 1st Columna (tiempo t) 
salida = loadedData(:, 2)'; % 2nd Columna (entrada u(k)) 
entrada = loadedData(:, 3)'; % 3rd Columna (tiempo t)

%% Normalizacion (0-1)

u_min = min(entrada);
u_max = max(entrada);
u_norm = (entrada - u_min) / (u_max - u_min);

y_min = min(salida);
y_max = max(salida);
y_norm = (salida - y_min) / (y_max - y_min);


%% Regresores para estimación

% Determino cantidad de regresores que voy a usar
cant_regresores = 2;

% Inicialización de matrices de regresores
y_k_n = zeros(cant_regresores, length(entrada)); % y_k-1 y y_k-2
u_k_n = zeros(cant_regresores, length(entrada)); % u_k y u_k-1

% Construcción de la matriz de regresores
for i = cant_regresores+1:length(entrada)
    
    % Regresores de salida
    for k = 1:cant_regresores
        y_k_n(k, i) = y_norm(i-k); % Se corrige el índice
    end

    % Regresores de entrada
    for k = 1:cant_regresores 
        u_k_n(k, i) = u_norm(i-k+1); % Se corrige el índice
    end

end

% Selección de la salida deseada
y_selec = y_norm(cant_regresores+1:end); % Se selecciona la salida real en el tiempo actual

% Construcción de la matriz de características (regresores)
x_reg = [y_k_n(:,cant_regresores+1:end)' u_k_n(:,cant_regresores+1:end)']; 


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

%% Parametros

% Dimensiones
[f c] = size(x_reg);

% Cantidad de neuronas por cada
cant_neu_c1 = c+2;
cant_neu_c2 = round(cant_neu_c1/2); % tener en cuenta que mi /2 afecta si mi numero de capas iniciales es entero o impar
cant_neu_c3 = 1;

% Inicialización de pesos:
% Primera capa oculta:
% 6 neuronas, osea que necesito 6 pesos sinapticos por entrada
w_c1 = rand(cant_neu_c1, c); % rand(c,cant_neu_c1) o rand(cant_neu_c1,c) dependiendo si quiero que mis filas o mis columnas sean mis entradas
[f_c1 c_c1] = size(w_c1); % me da matriz 6 x 4 osea: (6 neuronas) x (4 entradas)

% Segunda capa oculta:
w_c2 = rand(cant_neu_c2,f_c1); % considerar que las salidas de la primera capa son mis entradas de la segunda
[f_c2 c_c2] = size(w_c2); % me da matriz 3 x 4 osea: (3 neuronas) x (6 entradas)

% Capa final:
w_cout = rand(cant_neu_c3,f_c2);
[f_cout c_cout] = size(w_cout); % me da matriz 1 x 4 osea: (1 peso final) x (3 entradas)

% Inicialización del sesgo
b_c1 = rand(f_c1,1);
b_c2 = rand(f_c2,1);
b_cout = rand(f_cout,1);

% Inicializo iteraciones (# iter = cantidad de epocas pasadas)
iter = 0;
iter_max = 1000;

% Taza de aprendizaje
n = 0.00001;
emc=[];emc_eval=[];

%% Entrenamiento y validacion

while iter < iter_max
    
    % Entrenamiento
    for i = 1:length(x_entren)
       
        % Capa oculta 1
        Yp_1(:,i) = w_c1*x_entren(i,:)' + b_c1; 
        Yt_1(:,i) = sigmoide(Yp_1(:,i));
       
        % Capa oculta 2
        Yp_2(:,i) = w_c2*Yt_1(:,i) + b_c2;
        Yt_2(:,i) = sigmoide(Yp_2(:,i));
    
        % Capa final
        Yp_out(:,i) = w_cout*Yt_2(:,i) + b_cout;  
        %Yt_out(:,i) = sigmoide(Yp_out(:,i));
        Yt_out(:,i) = Yp_out(:,i);
    
        % Error
        e(:,i) = T_entren(:,i) - Yt_out(:,i);
        
        % Retropropagación
        e_23(:,i) = w_cout' * e(:,i);
        e_21(:,i) = w_c2' * e_23(:,i);

        % Actualizo pesos
        w_cout = w_cout + n * (e(:,i) .* der_sigmoide(Yp_out(:,i))) * Yt_2(:,i)';
        %w_cout = w_cout + n * (e(:,i) .* Yt_2(:,i)');
        w_c2 = w_c2 + n * (e_23(:,i) .* der_sigmoide(Yp_2(:,i))) * Yt_1(:,i)';
        w_c1 = w_c1 + n * (e_21(:,i) .* der_sigmoide(Yp_1(:,i))) * x_entren(i,:);
        
        % Actualizo sesgos
        b_cout = b_cout + n * (e(:,i) .* der_sigmoide(Yp_out(:,i)));
        b_c2 = b_c2 + n * (e_23(:,i) .* der_sigmoide(Yp_2(:,i)));
        b_c1 = b_c1 + n * (e_21(:,i) .* der_sigmoide(Yp_1(:,i)));
    
    end

    % Validación
    for i = 1:length(x_val)
       
        % Capa oculta 1
        Yp_1_val(:,i) = w_c1*x_val(i,:)' + b_c1; 
        Yt_1_val(:,i) = sigmoide(Yp_1_val(:,i));
       
        % Capa oculta 2
        Yp_2_val(:,i) = w_c2*Yt_1_val(:,i) + b_c2;
        Yt_2_val(:,i) = sigmoide(Yp_2_val(:,i));
    
        % Capa final
        Yp_out_val(:,i) = w_cout*Yt_2_val(:,i) + b_cout;  
        %Yt_out_val(:,i) = sigmoide(Yp_out_val(:,i));
        Yt_out_val(:,i) = Yp_out_val(:,i);

        
        % Error
        e_val(:,i) = T_val(:,i) - Yt_out_val(:,i);

    end
    
    iter = iter + 1;

    % Error cuadratico medio
    emc(iter)=(1/f)*sum(e.^2);
    emc_eval(iter)=(1/f)*sum(e_val.^2);

end

% Grafico error cuadratico medio
figure(1)
plot(emc)
title("EMC Entrenamiento")

figure(2)
plot(emc_eval)
title("EMC Validacion")

disp('Iteraciones:')
disp(iter)

disp('Pesos finales capa salida:')
disp(w_cout)
disp('Pesos finales capa oculta 2:')
disp(w_c2)
disp('Pesos finales capa oculta 1:')
disp(w_c1)

disp('Bias capa salida:')
disp(b_cout)
disp('Bias capa oculta 2:')
disp(b_c2)
disp('Bias capa oculta 1:')
disp(b_c1)

% disp('Deseados vs Entrenamiento:')
% prueba_entren = [T_entren' Yt_out'];
% disp(prueba_entren)
% 
% disp('Deseados vs Validacion:')
% prueba_val = [T_val' Yt_out_val'];
% disp(prueba_val)

%% Prueba

for k = 1:num_total_datos
    
    % Capa oculta 1
    Yp_1_prueba(:,k) = w_c1*x_reg(k,:)' + b_c1; 
    Yt_1_prueba(:,k) = sigmoide(Yp_1_prueba(:,k));
       
    % Capa oculta 2
    Yp_2_prueba(:,k) = w_c2*Yt_1_prueba(:,k) + b_c2;
    Yt_2_prueba(:,k) = sigmoide(Yp_2_prueba(:,k));
    
    % Capa final
    Yp_out_prueba(:,k) = w_cout*Yt_2_prueba(:,k) + b_cout;  
    Yt_out_prueba(:,k) = (Yp_out_prueba(:,k));

end

% Desnormalizar
Yt_prueba_desnormalizado = Yt_out_prueba * (y_max - y_min) + y_min;
y_selec_desnormalizado = y_selec * (y_max - y_min) + y_min;

% Comparar real vs prediccion
figure(3)
plot(Yt_prueba_desnormalizado,'b')
hold on
plot(salida, 'r', LineWidth=1.5)
plot(y_selec_desnormalizado, 'm')

prueb = [y_selec_desnormalizado', Yt_prueba_desnormalizado']


%% Sigmoides

function yT = sigmoide(yP)
    % Aplica la función sigmoide a cada elemento de la matriz
    yT = 1 ./ (1 + exp(-yP));
end

function yT = der_sigmoide(yP)
    % Aplica la derivada de la función sigmoide a cada elemento de la matriz
    sigmoidMatrix = 1 ./ (1 + exp(-yP));
    yT = sigmoidMatrix .* (1 - sigmoidMatrix);
end

