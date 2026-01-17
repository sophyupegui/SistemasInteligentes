% Cargar datos desde un archivo Excel
filename = 'Datos_segundo_escalon.xlsx';  % Cambia esto por el nombre de tu archivo
data = xlsread(filename);      % Lee todos los datos numéricos del Excel

% Extraer columna 1 y 2
columna1 = data(:,1);
columna2 = data(:,2);
columna3 = data(:,3);

% Corregir la escala dividiendo por 1 millón
columna1 = columna1 / 1e6;
columna2 = (columna2 -columna2(1))/ 1e6;
columna3 = columna3 / 1e6;

% Datos de entrada y salida
tiempo = columna1;   % Eje X en segundos
salida = columna2;   % Eje Y respuesta al escalón

% Determinar Ganancia K
entrada_escalon = 10;  % Cambia aquí si tu escalón fue otro valor
K = (salida(end) - salida(1)) / entrada_escalon;

% Estimar Tiempo muerto L
% Aproximación: el tiempo cuando la salida empieza a subir
idx_inicio = find(salida > salida(1) + 0.02*(salida(end)-salida(1)), 1);
L = tiempo(idx_inicio) - tiempo(1);

% Estimar Constante de Tiempo T
% Aproximación: tiempo donde la salida alcanza el 63.2% de su valor final
Y63 = salida(1) + 0.632*(salida(end) - salida(1));
idx_T = find(salida >= Y63, 1);
T = tiempo(idx_T) - tiempo(idx_inicio);

% Mostrar resultados
fprintf('Modelo estimado:\n');
fprintf('K = %.4f\n', K);
fprintf('L = %.4f segundos\n', L);
fprintf('T = %.4f segundos\n', T);
fprintf('\nFunción de Transferencia:\n');
fprintf('G(s) = %.4f / (%.4f s + 1) * exp(-%.4f s)\n', 3, T, 1.75);

% Datos originales
tiempo = columna1;   % Eje X en segundos
salida_real = columna2;  % Eje Y: respuesta medida

% Parámetros estimados (usa los que calculaste antes)
entrada_escalon = 10;  % valor del escalón aplicado
K = 1.75;
idx_inicio = find(salida_real > salida_real(1) + 0.02*(salida_real(end)-salida_real(1)), 1);
L = 3;
Y63 = salida_real(1) + 0.632*(salida_real(end) - salida_real(1));
idx_T = find(salida_real >= Y63, 1);
T = tiempo(idx_T) - tiempo(idx_inicio);

% Crear sistema en MATLAB
sys = tf(K, [T 1], 'InputDelay', L);

% Simular la respuesta con el mismo tiempo y escalón
entrada = entrada_escalon * ones(size(tiempo));
[respuesta_modelo, ~] = lsim(sys, entrada, tiempo);

% Graficar comparación
figure;
plot(tiempo, salida_real, 'b-', 'LineWidth', 1.5); hold on;
plot(tiempo, respuesta_modelo, 'r--', 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Salida');

title('Comparación: Datos Reales vs Modelo de Primer Orden con Retardo');
grid on;

