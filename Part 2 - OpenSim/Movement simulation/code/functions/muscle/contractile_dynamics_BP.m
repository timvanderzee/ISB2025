function[vM, Xdot] = contractile_dynamics_BP(a, FMltilda, x, vMT, kT, kP, cos_a, params, parms)

% v2: with dampening
lMo = params(2);
lTs = params(3);

%% define parameters
gamma = 100; % length scaling
delta = parms.Fscale;

%% simulate biophysical dynamics
eps = 1e-6;
FM = x(1);
Q0 = x(2);

Q0 = max(Q0, eps);
Q1 = FM / delta - parms.ps * Q0;
Q2 = x(3);
Non = x(4);
DRX = x(5);

% Q2 = max(Q2, 0);
 
% get geometry
% F = Q1 + parms.ps * Q0;
% dlse = parms.Lse_func(F, parms);
% kS =  parms.kse_func(dlse, parms);
kS = 100;
kT = (kT .* lMo(:)./lTs(:)) / (delta*gamma); % F0/lTs to F0/lMo to Q1max/h;
kP = kP / (delta*gamma); % F0/lTs to F0/lMo;

parms.vmtc = vMT / lMo * gamma; % l0/s to h/s
parms.act = a;
parms.Noverlap = FMltilda;
parms.no_tendon = 0;

% simulate biophysical dynamics
X = [Non Q0 Q1 Q2 DRX];
[Xd, Ld] = biophysical_dynamics(parms.Ca, X, kS, kP, kT, cos_a, parms);

vM = Ld / gamma .* lMo(:);

Xdot = [Xd(2); Xd(4); Xd(1); Xd(5)];

end