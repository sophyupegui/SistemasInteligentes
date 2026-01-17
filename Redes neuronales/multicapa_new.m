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
T_selec = y_norm(cant_regresores+1:end)'; % Se selecciona la salida real en el tiempo actual

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
x_train = x_reg(indices_entren,:); % Patron de entrenamiento
T_train = T_selec(indices_entren);

% Seleccionar los datos aleatorios validacion
x_val = x_reg(indices_val,:);
T_val = T_selec(indices_val);


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
w_c3 = rand(cant_neu_c3,f_c2);
[f_c3 c_c3] = size(w_c3); % me da matriz 1 x 4 osea: (1 peso final) x (3 entradas)

% Inicialización del sesgo
b_c1 = rand(f_c1,1);
b_c2 = rand(f_c2,1);
b_c3 = rand(f_c3,1);

% Inicializo iteraciones (# iter = cantidad de epocas pasadas)
iter = 0;
iter_max = 1000;

% Taza de aprendizaje
n = 0.01;

% Error permitido con base al error cuadrático medio
e_perm = 0.000001;

% Inicializar error cuadrático medio
EMC = 100; % Un valor inicial grande
EMC_hist = [];
EMC_V = 100;
EMC_V_hist = [];


%% Entrenamiento

while ((EMC > e_perm) && (iter < iter_max))
    
    iter = iter + 1;

    % Entrenamiento 
    for k = 1:length(x_train)

        % Capa oculta 1
        Yp_1 = w_c1 * x_train(k,:)' + b_c1; 
        Yt_1 = sigmoide(Yp_1);
       
        % Capa oculta 2
        Yp_2 = w_c2 * Yt_1 + b_c2;
        Yt_2 = sigmoide(Yp_2);
    
        % Capa final
        Yp_out = w_c3 * Yt_2 + b_c3;  
        Yt_out(k) = Yp_out; % sin sigmoide porque se satura
    
        % Error
        e(k) = T_train(k) - Yt_out(k);

        % Retropropagación del error
        e3 = e(k); 
        e2 = (w_c3' * e3); 
        e1 = (w_c2' * e2);

        % Actualización de pesos y bias
        w_c3 = w_c3 + n * e3 .* der_sigmoide(Yp_out) * Yt_2';
        b_c3 = b_c3 + n * e3 .* der_sigmoide(Yp_out);

        w_c2 = w_c2 + n * e2 .* der_sigmoide(Yp_2) * Yt_1';
        b_c2 = b_c2 + n * e2 .* der_sigmoide(Yp_2);
        
        w_c1 = w_c1 + n * e1 .* der_sigmoide(Yp_1) * x_train(k, :);
        b_c1 = b_c1 + n * e1 .* der_sigmoide(Yp_1);

    end

    % Calcular el error cuadrático medio del entrenamiento
    EMC = (1/num_datos_entrenamiento) * sum(e.^2);
    EMC_hist(iter) = EMC;

    % Validación
    for k = 1:length(x_val)

        % Capa oculta 1
        Yp_1_val = w_c1*x_val(k,:)' + b_c1; 
        Yt_1_val = sigmoide(Yp_1_val);
       
        % Capa oculta 2
        Yp_2_val = w_c2*Yt_1_val + b_c2;
        Yt_2_val = sigmoide(Yp_2_val);
    
        % Capa final
        Yp_out_val = w_c3*Yt_2_val + b_c3;  
        Yt_out_val(k) = Yp_out_val;

        % Error
        e_val(k) = T_val(k) - Yt_out_val(k); 

    end

    % Calcular el error cuadrático medio de validación
    EMC_V = (1/num_datos_validacion) * sum(e_val.^2);
    EMC_V_hist(iter) = EMC_V;

end


% Grafico error cuadratico medio
figure(1)
semilogx(EMC_hist)
title("EMC Entrenamiento")

figure(2)
semilogx(EMC_V_hist)
title("EMC Validacion")

disp('Iteraciones:')
disp(iter)

disp('Pesos finales capa salida:')
disp(w_c3)
disp('Pesos finales capa oculta 2:')
disp(w_c2)
disp('Pesos finales capa oculta 1:')
disp(w_c1)

disp('Bias capa salida:')
disp(b_c3)
disp('Bias capa oculta 2:')
disp(b_c2)
disp('Bias capa oculta 1:')
disp(b_c1)

disp('Deseados vs Entrenamiento:')
dif_entren = (abs(T_train - Yt_out')./ T_train)*100;
prueba_entren = [T_train Yt_out' dif_entren];
disp(prueba_entren)
 
disp('Deseados vs Validacion:')
dif_val = (abs(T_val - Yt_out_val')./ T_val)*100;
prueba_val = [T_val Yt_out_val' dif_val];
disp(prueba_val)


%% Prueba

for k = 1:num_total_datos
    x = x_reg(k, :)';
    Yt_1_test = sigmoide(w_c1 * x + b_c1);
    Yt_2_test = sigmoide(w_c2 * Yt_1_test + b_c2);
    Y_pred(k) = w_c3 * Yt_2_test + b_c3;
end

% Desnormalizar
Y_pred_real = Y_pred * (y_max - y_min) + y_min;

figure(3)
plot(salida, 'b')
hold on
plot(Y_pred_real, 'r')
plot(entrada, 'k')
legend('Real', 'Predicción', 'Entrada')
xlabel('Muestras')
ylabel('Salida')
title('Comparación de salida esperada vs. predicción')
grid on
hold off


%% Funcion

function yT = sigmoide(yP)
    yT = 1 ./ (1 + exp(-yP));
end

function ds = der_sigmoide(yP) % asi si uso Yp (como se define por derivada de sigmoide)
    % f*(1-f)
    sigmoid = 1 ./ (1 + exp(-yP));
    ds = sigmoid .* (1 - sigmoid);
end

% function ds = der_sigmoide(yT) % asi si uso Yt (como se define proceso)
%     % f*(1-f)
%     ds = yT .* (1 - yT))
% end
