%% Perceptron

clc
clear
close all

%%
x_1= [5.5 22 16.5 11 27.5 33 38.5 49.5 44 33 49.5 38.5 44 27.5 33];
x_2= [14 7 35 21 28 14 56 35 49 63 56 42 63 49 42];

T= [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1];

%% Patrones de entrenamiento

P=[x_1(1:11);x_2(1:11)] % uso los datos de 1-11 para entrenar
[f c]=size(P)

%%Inicialización de pesos
W=rand(1,2)

%%Inicialización del b
b=rand(1)
e=1;
iter=0;

%% Entrenamiento

while (numel(e(e~=0))>=1)
   iter=iter+1;
   
   for k=1:c
       yp=W*P(:,k)+b;
       if (yp>=0)
           yT(k)=1;
       else
           yT(k)=0;
       end

  % Calcular el error 
  e(k)=T(k)-yT(k);
  
  % Actualización de pesos y b 
  Delta_w=W+(e(k)*P(:,k)');
  Delta_b=b+e(k);
  W=Delta_w;
  b=Delta_b;

   end
end

% yT = T

%% Graficacion entrenamiento

% Recta discriminante / frontera
x1 = 0:0.1:100; %x1 != x_1 son DIFERENTES
m = -W(1)/W(2); %pendiente de la recta
br = -b/W(2);
x2 = m*x1+br;

figure(1)
%title("Entrenamiento") !!!!!!!!!!!!!
plot(x1,x2)
hold on
plot(x_1(1:11),x_2(1:11),'o')


%% Prueba

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


