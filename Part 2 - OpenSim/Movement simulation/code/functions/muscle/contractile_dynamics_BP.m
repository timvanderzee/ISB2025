function[vM, Xdot] = contractile_dynamics_BP(a, FMltilda, x, vMT, kT, kP, cos_a, params)

% v2: with dampening
lMo = params(2);
lTs = params(3);

%% define parameters
gamma = 100; % length scaling
delta = 1.9207;

parms.w = 0.2;

parms.f = 1e3;
parms.k11 = 55.3664;
parms.k12 = 2;
parms.k21 =  451.0874;
parms.k22 =  0.2328;

parms.f = 1.1711e+03;
parms.k11 =  14.2557;
parms.k12 = 2;
parms.k21 =   613.9038;
parms.k22 =  0;

parms.vmtc = 0;
parms.Ca = 1;
parms.ps = 1;

parms.J1 = 89.2823;
parms.J2 = 4 * parms.J1;
parms.JF = 70.8;

parms.koop = 19.4;
parms.kon = 34.7;
parms.koff = 80;

parms.gaussian.IGef{1} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)));
parms.gaussian.IGef{2} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*(c(2)-c(3)*k(2)/2);
parms.gaussian.IGef{3} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*((c(2)-c(3)*k(2)/2).^2+c(3)/2);

%% simulate biophysical dynamics
eps = 1e-6;
FM = x(1);
Q0 = x(2);
Q1 = FM / delta - parms.ps * Q0;
Q2 = x(3);
Non = x(4);
DRX = x(5);

Q0 = max(Q0, eps);
Q2 = max(Q2, eps);
 
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