%% Recosido simulado

clc
close all
clear

%% Modelo

Num = 0.8;
Den = [85 1];
Gs = tf(Num, Den, 'inputdelay', 1)

%% Controlador inicial - Solución inicial (So)

% Al menos debe cumplir con parametros de Ziegler-Nichols
Kp = 1.5;
Ki = 0.2;
Kd = 0.2;

% Creo controlador inicial
Cs_inicial = pid(Kp, Ki, Kd);

% Cierro lazo
Cs_G_s_inicial = series(Cs_inicial, Gs);
Control_inicial = feedback(Cs_G_s_inicial, 1); % podria poner esto como feedback(Cs*Gs, 1) y me ahorro el series de arriba

% Respuesta temporal
step(Control_inicial)
C_aux_inicial = stepinfo(Control_inicial) 

% Evaluo primer controlador
peso = 0.5; % a criterio mio
costo_inicial = peso*C_aux_inicial.SettlingTime + C_aux_inicial.Overshoot;

%% Inicializacion de parametros del RS

%Temperatura 
T0 = 100; % Temp inicial - asumo ambiente en este caso
T = T0; % Temp actual
Tmin = 0.0001; % Minima temp a la que puede llegar
Enfria = 0.2; % Ratio de enfriamiento

% Iteraciones
iter_max = 50;

% Almacenar mejoras
mejor_s = [Kp, Ki, Kd]; % mejor solucion (actualmente solo tengo la inicial entonces es esta, ira mejorando cada iteracion)
mejor_costo = costo_inicial; % mejor solucion (actualmente sparo de olo tengo la inicial entonces es esta, ira mejorando cada iteracion)
f_costos = [];

while T > Tmin
    
    for k = 1:iter_max
        
        % Nuevo controlador
        coef = 0.6; % el coeficiente lo voy variando a mi criterio dependiendo de que tanto quiero afectar mis ganancias
        beta_K = coef*(rand(size(mejor_s))); % aleatorizacion 
        C_new = mejor_s + beta_K; 
        Cs = pid(C_new(1), C_new(2), C_new(3));
        Cs_G_s = series(Cs, Gs);
        Control = feedback(Cs_G_s, 1); 

        % Respuesta temporal
        C_aux = stepinfo(Control) 
        
        if (isnan(C_aux.SettlingTime) || isnan(C_aux.Overshoot))
            % restriccion: si mi controlador se vuelve inestable, no recibirlo
            costo = inf;
        else
            peso = 0.5; % a criterio mio
            costo = peso*C_aux.SettlingTime + C_aux.Overshoot;
        end

        delta = costo - mejor_costo;

        % Estrategia para escapar minimo local
        if (delta < 0) || (rand<exp(-delta/T)) %nuevo controlador vs anterior: (es mejor, lo acepto) o (si es peor, lo acepto con una probabilidad)
            mejor_sol = C_new; % Traigo nueva solucion
            mejor_costo = costo;
        end
    
        f_costos = [f_costos mejor_costo];

    end % for

    k
    T=T*Enfria;

end % while


% Graficar controlador optimizado
C_opt = pid(mejor_sol(1), mejor_sol(2), mejor_sol(3));
Gop = feedback(C_opt*Gs, 1);
step(Gop)
hold on

% Control optimizado vs controlador inicial
step(Control_inicial)
hold off

% Costos
figure()
plot(f_costos)



