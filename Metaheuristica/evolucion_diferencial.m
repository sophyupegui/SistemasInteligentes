%% Evolución Diferencial (DE) para Sintonización de un Controlador PID

clc
close all
clear

%% Modelo del sistema

Tau = 125;
K = 0.8;
Theta = 2.1; 

Num = [K];
Den = [Tau 1];
Gs = tf(Num, Den, 'inputdelay', Theta);

%% Inicialización de parámetros de Evolución Diferencial

N = 20;         % Tamaño de la población
D = 3;          % Numero de parámetros (Kp, Ki, Kd)
CR = 0.9;       % Tasa de cruce (0 - 1)
F_min = 0.9;    % Minimo valor de factor de escala (min: 0)
F_max = 1.2;    % Maximo valor de factor de escala (max: 2)

% Límites de búsqueda para Kp, Ki, Kd
lb = [0.1, 0.01, 0.01];
ub = [10, 5, 5]; 

% Inicializar población aleatoria dentro de los límites
pob = lb + rand(N, D) .* (ub - lb);
fitness = zeros(N, 1);

%% Controlador inicial

% Evaluar función de costo para la población inicial
peso = 0.5; % a criterio mio
for i = 1:N
    fitness(i) = objective_function(pob(i, :), Gs, peso);
end

% Guardar el primer controlador generado y su costo
[first_fitness, first_idx] = min(fitness);
first_solution = pob(first_idx, :); % q0, q1, q2 iniciales

% Grafica
Cs_inicial = pid(first_solution(1), first_solution(2), first_solution(3));
Control_inicial = feedback(series(Cs_inicial, Gs), 1);
figure(1)
step(Control_inicial);
title('Controlador inicial');

%% Evolución Diferencial (Iteraciones)

% Inicializacion
MaxIt = 50; % Número máximo de iteraciones    

% Proceso iterativo ED
for iter = 1:MaxIt

    for i = 1:N
        
        % Seleccionar 3 individuos aleatorios distintos
        idx = randperm(N, D+1);
        idx(idx == i) = [];          % Asegurar que i no esté en la selección
        r1 = pob(idx(1), :);
        r2 = pob(idx(2), :);
        r3 = pob(idx(3), :);
        
        % Factor de escala (mutación)
        F = F_min + (F_max - F_min) * rand;

        % Mutación
        V = r1 + F * (r2 - r3);
        V = max(min(V, ub), lb);  % Restringir rango
        
        % Cruce
        j_rand = randi(D); % Asegurar que al menos un valor cambie

        for j = 1:D
            if rand <= CR || j == j_rand
                U(j) = V(j);
            else
                U(j) = pob(i, j);
            end
        end

        % Asegurar que U esté dentro de los límites
        U = max(min(U, ub), lb); % Restringuir limites
        
        % Evaluar el nuevo vector de prueba
        new_fitness = objective_function(U, Gs, peso);
        
        % Selección: Reemplazar si la nueva solución es mejor
        if new_fitness < fitness(i)
            pob(i, :) = U;
            fitness(i) = new_fitness;
        end
        
    end
    
    % Guardar la mejor solución de esta iteración
    [best_fitness, best_idx] = min(fitness);
    best_solution = pob(best_idx, :);

    % Guardar el mejor fitness en el historial
    fitness_history(iter) = best_fitness;
    
    % Mostrar progreso
    fprintf('Iteración %d: Mejor fitness = %.5f\n', iter, best_fitness);

end

%% Resultados Finales

fprintf('\nMejor solución encontrada: Kp = %.5f, Ki = %.5f, Kd = %.5f', best_solution);
fprintf('\nCosto del mejor controlador: %.5f\n', best_fitness);

% Comparación con el primer controlador generado
fprintf('\nPrimer controlador generado:\nKp = %.5f, Ki = %.5f, Kd = %.5f', first_solution);
fprintf('\nCosto del primer controlador: %.5f\n', first_fitness);

% Controlador final
Cs_final = pid(best_solution(1), best_solution(2), best_solution(3));
Control_final = feedback(series(Cs_final, Gs), 1);

% Grafica comparacion de controladores
figure(2)
step(Control_final)
hold on
title('Controlador inicial vs optimizado')
step(Control_inicial)
legend('Controlador final','Controlador inicial', Location='southeast')
hold off

% Gráfica de minimización
figure(3)
plot(1:MaxIt, fitness_history, 'LineWidth', 1.5)
xlabel('Iteraciones');
ylabel('Valor de la función de costo');
title('Evolución de la función de costo');
grid on;

%% Función Objetivo (Evaluación del Controlador)
function costo = objective_function(q, Gs, peso)
    Kp = q(1);
    Ki = q(2);
    Kd = q(3);
    
    % Crear y evaluar el controlador PID
    Cs = pid(Kp, Ki, Kd);
    Control = feedback(series(Cs, Gs), 1);
    
    % Obtener métricas de desempeño
    C_aux = stepinfo(Control);
    
    if isnan(C_aux.SettlingTime) || isnan(C_aux.Overshoot)
        costo = inf; % Sistema es inestable
    else
        costo = peso * C_aux.SettlingTime + C_aux.Overshoot;
    end
   
end 
