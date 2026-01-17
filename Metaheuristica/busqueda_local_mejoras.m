%% Busqueda local - n mejoras

clc
close all
clear

%% Modelo - Búsqueda local con variación aleatoria

% Función de transferencia (sistema continuo)
Num = 0.264;
Den = [112 1];
Gp_s = tf(Num, Den);

% Grafico
figure(1)
step(Gp_s)
hold on

% Sistema discreto
t_muestreo = 20; % discretizo con tiempo muestreo 20
HGp_z = c2d(Gp_s, t_muestreo)  
step(HGp_z) % grafico en comparacion del continuo
title("Modelo: Continuo vs Discreto")
legend("Continuo", "Discreto", Location="southeast")
hold off

%% Solucion inicial (controlador discreto inicial)

% Coeficientes iniciales
q0 = 1.8;  % positivo
q1 = -1.1; % negativo
q2 = 0.01; % |q0| > |q1| >> |q2|

% Creo controlador discreto inicial
% Cz = (q0*z^2 + q1*z + q2) / ((z*(z-1))
Cz_inicial = tf([q0 q1 q2], [1 -1 0], t_muestreo); %si aumento q1 se vuelve 1st grado, si disminuyo se vuelve 2nd
G_z_inicial = series(Cz_inicial, HGp_z); % modelo de planta con controlador
Y_z_inicial = feedback(G_z_inicial,1); % realimentacion

% Grafico
figure(2)
step(Y_z_inicial)
title("Rspta al Escalon del Controlador Inicial")

% Caracteristicas de respuesta temporal (del controlador inicial)
R_t = stepinfo(Y_z_inicial)

% Funcion objetivo
[tss Mp] = costo_j(R_t);
f_obj_inicial = 0.01*tss + Mp;

%% Genero vecindarios / controladores

% Cambio maximo de coeficientes 
beta_q0 = 0.05;
beta_q1 = 0.01;
beta_q2 = 0.002;

% Número de iteraciones
iter_max = 20;

% Mejoras a almacenar
num_mejoras = 5; % Variable para definir la cantidad de mejoras deseadas
mejoras = [];
cont = 0;

while cont < num_mejoras && length(mejoras) < iter_max
    
    cont = cont + 1;
    
    % Generar perturbaciones aleatorias
    delta_q0 = (rand - 0.5) * 2 * beta_q0; % Rango de aleatorios: [−beta_q0, beta_q0]
    delta_q1 = (rand - 0.5) * 2 * beta_q1;
    delta_q2 = (rand - 0.5) * 2 * beta_q2;
    
    % Crear vecindario con valores perturbados
    Controladores = [q0+delta_q0, q1, q2;
                     q0, q1+delta_q1, q2;
                     q0, q1, q2+delta_q2];
    
    % Recorrer controladores - Evaluar cada conjunto de coeficientes
    for i = 1:length(Controladores)
        
        Cz = tf(Controladores(i,:), [1 -1 0], t_muestreo); % [Controladores(i,1), Controladores(i,2), Controladores(i,3)]

        % Evaluar cada controlador
        G_z = series(Cz, HGp_z); % modelo de planta con controlador
        Y_z = feedback(G_z,1); % realimentacion
        R_t = stepinfo(Y_z); % respuesta temporal
        [tss Mp] = costo_j(R_t); % calculo tss y Mp
        f_obj(i) = 0.01*tss + Mp; % funcion objetivo

    end
    
    % Seleccionar el mejor
    [f_min, indice] = min(f_obj);
    
    % Si mejora, actualizar coeficientes
    if f_min < f_obj(1)
        
        % Solucion final
        q0 = Controladores(indice, 1);
        q1 = Controladores(indice, 2);
        q2 = Controladores(indice, 3);

    end
    
    mejoras = [mejoras f_min];

end

% Display
f_obj
f_min

% Controlador final
Cz_final = tf([q0, q1, q2], [1 -1 0], t_muestreo); % con coeficientes finales
G_z_final = series(Cz_final, HGp_z);
Y_z_final = feedback(G_z_final, 1);

% Grafico
figure(3)
step(Y_z_final)
hold on
step(Y_z_inicial)
title("Controlador: Final vs Inicial")
legend("Controlador Final", "Control Inicial", Location="southeast")
hold off

% Respuesta temporal
R_t = stepinfo(Y_z_final)

% Funcion objetivo
[tss Mp] = costo_j(R_t);
f_obj = 0.01*tss + Mp;

% Graficar mejoras
figure(4)
plot(mejoras)
title("Evolución de la función objetivo")
xlabel("Iteración")
ylabel("f_{obj}")

% Ley de control
u_z_final = series(Cz_final, (1 - Y_z_final));
u_z_inicial = series(Cz_inicial, (1 - Y_z_inicial));

figure(5)
step(u_z_final)
hold on
step(u_z_inicial)
step(Y_z_final)
step(Y_z_inicial)
title("Ley de control: Final vs Inicial")
legend("Ley de control Final", "Ley de control Inicial", "Controlador Final", "Control Inicial", Location="southeast")
hold off

%% Función de costo
function [tss, Mp] = costo_j(R_t)
    Mp = R_t.Overshoot;
    tss = R_t.SettlingTime;
end
