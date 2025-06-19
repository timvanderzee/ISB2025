function[Qd, vM, Nond, DRXd] = contractile_dynamics_BP_v3(a, FMltilda, x, vMT, kT, kP, cos_alpha, params)

% v2: with dampening
lMo = params(2);
lTs = params(3);

%% spatial integrals
gamma = 100; % length scaling

delta =   2.018;
parms.w = 0.2;
parms.k11 = 106.5;
parms.k12 = 2;
parms.k21 = 480.6;
parms.k22 = .23;
parms.f = 1000;
parms.vmtc = 0;
parms.d = 1e-05;
parms.Ca = 1;

parms.J1 = 0.96;
parms.J2 = 4 * parms.J1;
parms.JF = 70.8;

parms.koop = 19.4;
parms.kon = 34.7;
parms.koff = 80;

parms.gaussian.IGef{1} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)));
parms.gaussian.IGef{2} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*(c(2)-c(3)*k(2)/2);
parms.gaussian.IGef{3} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*((c(2)-c(3)*k(2)/2).^2+c(3)/2);

%% spatial integrals
eps = 1e-6;
Q(1,1) = x(2);
Q(1) = max(Q(1),eps);
Q(2,1) = x(1) / delta - Q(1);
Q(3,1) = max(x(3), eps);

X = [x(5) Q(1) Q(2) Q(3) [] x(6) lMT];

parms.kT = (kT .* lMo(:)./lTs(:)) / (delta*gamma); % F0/lTs to F0/lMo to Q1max/h;
parms.kP = kP / (delta*gamma); % F0/lTs to F0/lMo;
parms.vMT = vMT / lMo * gamma; % l0/s to h/s

parms.act = a;
parms.FMltilda =FMltilda;
parms.cos_alpha = cos_alpha;

[Xd, Ftot, Ld] = ripping_model_func_exp_v2(t, X, parms, Ca);

vM = Ld / gamma .* lMo(:);

Nond = Xd(1);
Qd = Xd(2:4);
DRXd = Xd(6);

end