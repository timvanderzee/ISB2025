function [vts, Fts, toc, idF, idFd] = design_length_input_vector(vmax, RT, V_rel, N)

%% Design velocity vector
vt = [0 -.2 0 .2 0 -.4 0 .4 0 -.6 0 .6 0 -.8 0 .8 0 V_rel -V_rel 0 V_rel] * vmax;

% for velocity part, we just need to evaluate until evaluation time
idv = 1:(length(vt)-5);

ts = .2 * ones(size(vt));
for i = 1:length(idv)
    ts(idv(i)) = .12 / max(abs(vt(idv(i))), .6);
end

idS = (length(vt)-4):length(vt);
ts(idS) = [.3 .1 .1 RT .1];
Ts = [0 cumsum(ts)];
toc = linspace(0,sum(ts),1000);

% velocity vector
vts = zeros(1,N);
for i = 1:(length(Ts)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

% indices
idF = nan(1, length(idv)+1);
idF(1) = find(toc < Ts(2), 1, 'last');
for i = 1:length(idv)
    idF(i+1) = find(toc < Ts(idv(i)+1), 1, 'last');
end

idFd = [find(toc > Ts(end-4),1); find(toc > Ts(end-1),1)];

%% get the target force
d = [ -0.3183   -8.1492   -0.3741    0.8856];
Fts = d(1) * log((d(2)*vts/vmax + d(3)) + sqrt((d(2)*vts/vmax + d(3)).^2 + 1)) + d(4);

end