%% Perceptron

clc
clear
close all

%% Defino parametros

% Entradas
x1 = [5.5 22 16.5 11 27.5 33 38.5 49.5 44 33 49.5 38.5 44 27.5 33];
x2 = [14 7 35 21 28 14 56 35 49 63 56 42 63 49 42];
x = [x1; x2]; 

% Target outputs
T = [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];

% Inicializacion de pesos rand(filas, columnas)
w = rand(1,2); %(1, cantidad de entradas)

% Bias (sesgo)
b = rand(1);

% Error
e = 1; 

%% Seleccionar muestra de datos

% Número total de datos
num_total_datos = length(x1);

% 70% del total de datos
num_datos_seleccionados = round(0.7 * num_total_datos);

% 15% validación
num_datos_validacion = round(0.15 * num_total_datos);

% Generar índices aleatorios para seleccionar datos
indices_aleatorios = randperm(num_total_datos, num_datos_seleccionados);

% Seleccionar los datos aleatorios
x_seleccionados = x(:,indices_aleatorios); % Patron de entrenamiento
x1_seleccionados = x_seleccionados(1,:);
x2_seleccionados = x_seleccionados(2,:);
T_seleccionados = T(indices_aleatorios); % Asegúrate de seleccionar también las salidas correspondientes

% datdat = [x1_seleccionados' x2_seleccionados' T_seleccionados'] % Verificar que si sean valores acorde

[f c] = size(x_seleccionados); %f es el número de filas y c el número de columnas de P

%% (como profe hizo en clase)

x_seleccionados = [x1(1:11); x2(1:11)];
T_seleccionados = T(1:11);
[f c] = size(x_seleccionados); %f es el número de filas y c el número de columnas de P

%osea lo que esta mal es ese rand 

%% Calcular salida parcial

% Iteracion 
iter = 0; % inicial
max_iter = 1000; % límite de iteraciones

% Mientras E>Ep haga: (Entrenamiento perceptrón) --> pseudo codigo
while (numel(e(e~=0))>=1) && (iter < max_iter)  % mientas que la cantidad de elementos de e que sean diferentes de cero sea mayor o igual a 1
    % se utiliza para determinar si hay algún error diferente de cero en el vector de errores e
    % Por ejemplo, si e = [0, 1, -1, 0, 2], entonces e ~= 0 produce [false, true, true, false, true]
    % e(e ~= 0): Utiliza la indexación lógica para seleccionar solo los elementos de e que son diferentes de cero.
    % numel(e(e ~= 0)): La función numel cuenta el número de elementos en la matriz o vector. Aquí, cuenta cuántos elementos de e son diferentes de cero.
    % Siguiendo con el ejemplo, numel([1, -1, 2]) produce 3.
   
    for k = 1:c % desde 1 hasta la cantidad de datos total (podria reeemplazarlo por length(x_seleccionados)??)
        % para K=1 hasta cantidad de datos (pseudo)
        % : significa "todos los elementos" a lo largo de la dimensión especificada.

        % Evaluar salida parcial Yp
        Yp(k) = w*x_seleccionados(:,k) + b; 

        % Evaluar YT  -> este 0 o 1 es por mi limitador duro
        if (Yp(k) >= 0)
            YT(k) = 1;
        else 
            YT(k) = 0;
        end
        
        % Calcular el error e
        e(k) = T_seleccionados(k) - YT(k); 

        % Acualizamos pesos y bias
        Delta_w = w + (e(k)*x_seleccionados(:,k)');
        Delta_b = b + e(k);
        w = Delta_w;
        b = Delta_b;

    end %for
  
    iter = iter+1; % aumento iteracion
   
end %while

YT
T_seleccionados

%% Graficacion entrenamiento (ARREGLAR PARA QUE CONCUERDEN MIS VARIABLES)

% Recta discriminante / frontera
x_1 = 0:0.1:100; %x1 != x_1 son DIFERENTES
m = -W(1)/W(2); %pendiente de la recta
br = -b/W(2);
x_2 = m*x_1+br;

figure(1)
%title("Entrenamiento") !!!!!!!!!!!!!
plot(x_1,x_2)
hold on
plot(x1(1:11),x2(1:11),'o')


%% Prueba

% indices_aleatorios_validacion = randperm(num_total_datos, num_datos_validacion);

% como mi modelo de perceptron ya fue entrenado en el while, ya puedo
% entregarle el resto de los datos para que los clasifique (recordar que
% hice el entrenamiento con solo los datos del 1-11)

% usualmente la prueba es con un 30% de los datos, pero si a mi entranamiento
% le fue bien, puedo tirar de una sola vez el 100% de mis datos

P = [x_1(1:end);x_2(1:end)]; % 1:end = 100% ; 11:end = 30%
[f c] = size(P);

for k=1:c
    yp=W*P(:,k)+b; %aqui ya uso mi W final y mi b final
    if (yp>=0)
        yT(k)=1;
    else
        yT(k)=0;
    end
end

yT
T %verifico que T = yT

%% Grafica prueba

figure(2)
hold on
for k=1:length(x_1)
    if yT(k)==1
        plot(P(1,k),P(2,k),'b*') %patron en primera y segunda posicion y en la posicion k, es como x_1 y x_2 pero con el patron
    else
        plot(P(1,k),P(2,k),'ro')
    end
end

plot(x1,x2,"LineWidth",1,"Color","m") %misma linea de frontera que ya habia creado antes
text(15,20,"Clase1") %poner texto en la posicion x,y 
text(50,60,"Clase2")
text(40,10,"Recta discriminante / Frontera de decision")







