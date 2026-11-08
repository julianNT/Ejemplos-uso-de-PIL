function comparar_pil_normal(tFinal)
%COMPARAR_PIL_NORMAL Compara las ejecuciones Normal y PIL de pid_aw_lazo.
%   Ejecuta el modelo con el bloque Model en modo Normal y luego en modo
%   PIL, imprime la diferencia maxima de u, v e y, y grafica las diferencias.
%   Al terminar, deja el bloque en el modo que tenia.
%
%   comparar_pil_normal        simula 5 s
%   comparar_pil_normal('2')   simula 2 s
%
%   La placa debe estar conectada y el modelo configurado con su puerto COM.

if nargin < 1
    tFinal = '5';   % Tiempo simulado [s]
end

modelo = 'pid_aw_lazo';
senales = {'u', 'v', 'y'};
modos = {'Normal', 'Processor-in-the-loop (PIL)'};

load_system(modelo);
bloque = find_system(modelo, 'SearchDepth', 1, 'BlockType', 'ModelReference');
bloque = bloque{1};

modoOriginal = get_param(bloque, 'SimulationMode');
restaurar = onCleanup(@() set_param(bloque, 'SimulationMode', modoOriginal));

% Ejecuta el modelo en cada modo
salidas = cell(1, numel(modos));
for i = 1:numel(modos)
    set_param(bloque, 'SimulationMode', modos{i});
    fprintf('Ejecutando en modo: %s\n', get_param(bloque, 'SimulationMode'));
    inicio = tic;
    salidas{i} = sim(modelo, 'StopTime', tFinal);
    fprintf('  Duracion: %.1f s\n', toc(inicio));
end

% Deja los resultados en el workspace, para no repetir la simulacion si algo falla
assignin('base', 'salidasComparacion', salidas);

% Compara cada senal
tiempos = cell(1, numel(senales));
diferencias = cell(1, numel(senales));
maximos = zeros(1, numel(senales));
fprintf('\nDiferencia maxima entre Normal y PIL:\n');
for i = 1:numel(senales)
    normal = obtenerSenal(salidas{1}.logsout, senales{i});
    pil = obtenerSenal(salidas{2}.logsout, senales{i});

    if numel(normal.Data) ~= numel(pil.Data)
        error('La senal %s tiene %d muestras en Normal y %d en PIL.', ...
            senales{i}, numel(normal.Data), numel(pil.Data));
    end

    tiempos{i} = normal.Time(:);
    diferencias{i} = normal.Data(:) - pil.Data(:);
    maximos(i) = max(abs(diferencias{i}));
    fprintf('  %s: %g\n', senales{i}, maximos(i));
end

% Grafica las diferencias
figure('Name', 'Diferencia Normal - PIL');
for i = 1:numel(senales)
    subplot(numel(senales), 1, i);
    plot(tiempos{i}, diferencias{i});
    grid on;
    ylabel(senales{i});
    title(sprintf('%s: diferencia maxima = %g', senales{i}, maximos(i)));
end
xlabel('Tiempo [s]');
end

function valores = obtenerSenal(registro, nombre)
% Busca la senal registrada e indica cuales existen si no la encuentra.
elemento = registro.getElement(nombre);
if isempty(elemento)
    error('No se encontro la senal %s. Senales registradas: %s', ...
        nombre, strjoin(registro.getElementNames, ', '));
end
valores = elemento.Values;
end
