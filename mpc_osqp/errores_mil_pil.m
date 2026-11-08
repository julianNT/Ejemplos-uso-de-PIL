function T = errores_mil_pil(idMIL, idPIL, senales)
% ERRORES_MIL_PIL  Errores absoluto, relativo y porcentual entre dos runs del SDI.
%   T = errores_mil_pil()                  usa los 2 ultimos runs (penultimo=MIL, ultimo=PIL)
%   T = errores_mil_pil(idMIL, idPIL)      runs elegidos por ID
%   T = errores_mil_pil(idMIL, idPIL, {'x(1)','x(2)','V','iter'})
%
% Senales vectoriales: 'x(2)' = canal 2 de la senal 'x'.
% Referencia = MIL (a). Comparado = PIL (b).
%   e_abs  = max |a-b|                       [unidad de la senal]
%   e_fs   = max|a-b| / max|a| * 100         [% del maximo de la senal]
%   e_rel  = max |a-b|/|a| * 100             [%] (solo donde |a| > 1% del maximo)
%   e_2    = ||a-b||_2 / ||a||_2 * 100       [%] (error global, norma L2)

ids = Simulink.sdi.getAllRunIDs;
if nargin < 2, idMIL = ids(end-1); idPIL = ids(end); end
if nargin < 3, senales = {'x(1)','x(2)','V','iter'}; end

rA = Simulink.sdi.getRun(idMIL);
rB = Simulink.sdi.getRun(idPIL);

n = numel(senales);
e_abs = nan(n,1); e_fs = nan(n,1); e_rel = nan(n,1); e_2 = nan(n,1);
figure('Name','MIL vs PIL');
for k = 1:n
    [a, t, okA] = leer(rA, senales{k});
    [b, tb, okB] = leer(rB, senales{k});
    if ~okA || ~okB
        warning('Senal %s no encontrada en ambos runs.', senales{k});
        disp('Senales en MIL:'); disp(nombres(rA));
        continue
    end
    b = interp1(tb, b, t, 'previous', 'extrap');   % mismos instantes (ZOH)

    d = abs(a - b);
    e_abs(k) = max(d);
    e_fs(k)  = max(d) / max(abs(a)) * 100;
    m = abs(a) > 0.01*max(abs(a));                 % evita dividir por ~0
    e_rel(k) = max(d(m) ./ abs(a(m))) * 100;
    e_2(k)   = norm(a-b) / norm(a) * 100;

    subplot(n,1,k); plot(t, d); grid on
    ylabel(['|MIL-PIL| ' senales{k}]);
end
xlabel('t [s]');

T = table(senales(:), e_abs, e_fs, e_rel, e_2, ...
    'VariableNames', {'Senal','e_abs','e_fs_pct','e_rel_pct','e_2_pct'});
disp(T)
writetable(T, 'errores_mil_pil.csv');
end

function [y, t, ok] = leer(run, nombre)
% Devuelve la senal como vector columna double. 'x(2)' -> canal 2 de 'x'.
y = []; t = []; ok = false;
canal = 1;
s = getSignalsByName(run, nombre);
if isempty(s)
    tok = regexp(nombre, '^(.*)\((\d+)\)$', 'tokens', 'once');
    if isempty(tok), return, end
    s = getSignalsByName(run, tok{1});
    canal = str2double(tok{2});
    if isempty(s), return, end
end
ts = s(1).Values;
t  = ts.Time(:);
d  = double(ts.Data);
if ndims(d) > 2, d = squeeze(d); end
if size(d,1) ~= numel(t), d = d.'; end            % deja tiempo en filas
if canal > size(d,2), return, end
y  = d(:, canal);
ok = true;
end

function c = nombres(run)
n = run.SignalCount;
c = cell(n,1);
for i = 1:n
    s = getSignalByIndex(run, i);
    c{i} = s.Name;
end
end
