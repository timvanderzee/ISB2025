function [vts, Fts, toc, idF, idFd] = design_length_input_vector(vmax, N)

%% Design velocity vector
% vmax = 5;

vt = [0 -2 2 -4 4 -6 0 5 -5 0 5]/10 * vmax;

% for velocity part, we just need to evaluate until evaluation time
idv = [2 3 4 5 6];

ts = .2 * ones(size(vt));
for i = 1:length(idv)
    ts(i+1) = .12 / abs(vt(idv(i)));
end

idS = 7:length(vt);
ts(idS) = [.3 .1 .1 .1 .1];
Ts = [0 cumsum(ts)];

toc = linspace(0,sum(ts),1000);
% N = length(toc);

% model constraints
Ns = floor(linspace(0, N, length(vt)+1));

vts = zeros(1,N);
for i = 1:(length(Ns)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

% close all
% figure(1)
% plot(toc, vts); hold on

%
idF = nan(1, length(idv)+1);

idF(1) = find(toc < Ts(2), 1, 'last');
for i = 1:length(idv)
    idF(i+1) = find(toc < Ts(i+1), 1, 'last');
end

% plot(toc(idF),vts(idF),'o')

idFd = [find(toc > Ts(end-4),1); find(toc > Ts(end-1),1)];
% plot(toc(idFd), vts(idFd),'x')

%% get the force
Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilda = linspace(0,2);
vH = vmax/e2*(sinh((FMvtilda-e4)/e1)-e3); % can be inverted = simpler
Fts = interp1(vH, FMvtilda, vts);

end