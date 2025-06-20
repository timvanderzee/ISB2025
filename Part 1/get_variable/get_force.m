function[y] = get_force(t, Ts, F, lmtc, parms)

% get active and passive force
F_pas = parms.Fpe_func(lmtc, parms);
F_act = F * parms.Fscale;
Ftot = F_act(:) + F_pas(:);

% subtract the first 3 intervals (center around 2nd stretch)
tnew = t - sum(Ts(1:4));

[~,id] = unique(tnew);

% store force
y.ti = min(tnew):.001:max(tnew);
y.Fi_act = interp1(tnew(id), F_act(id), y.ti);
y.lmtc = interp1(tnew(id), lmtc(id), y.ti);
y.Fi_pas = interp1(tnew(id), F_pas(id), y.ti);
y.Fi = interp1(tnew(id), Ftot(id), y.ti);

% average over the 1 s before the first stretch
id1 = tnew > (-sum(Ts(2:4)) - 1) & tnew < -sum(Ts(2:4));
y.F0 = mean(Ftot(id1));

% SRS
dt = .01;

Ts_rel = cumsum(Ts) - sum(Ts(1:4)); % relative to second stretch

tnew = t - sum(Ts(1:4));

id1 = isfinite(tnew) & tnew >= 0 & tnew < dt;

% fit a linear on the force-length
if sum(id1)>1
    p = polyfit(lmtc(id1), Ftot(id1), 1);
else
    p = [nan nan];
end

% baseline force
id2 = tnew < Ts_rel(1) & tnew > (Ts_rel(1)-5); 

y.Fpre =  mean(Ftot(id2),'omitnan');
y.SRS = p(1);


end