function[xdot, Ftot, Ld] = fiber_dynamics(t, x, parms, Ca)

% retrieve states
Non     = x(1);
Q       = x(2:end-3);
DRX     = x(end-2);
lmtc    = x(end-1);
L       = x(end-0);

% get force
[Q0, Q1, parms.xi] = force_from_distribution(Q, L, parms);

eps = 1e-6;
Q0 = max(Q0, eps);
FXB = Q0 + Q1;

% make sure values are within range
Non = max(Non, eps);

% get geometry
dLS = parms.Lse_func(FXB, parms);
kS = parms.kse_func(dLS, parms);
kT = 1000;
kP = parms.kpe;
cos_a = 1;

% simulate biophysical dynamics
X = [Non; Q(:); DRX];
 
[Xd, Ld] = biophysical_dynamics(Ca, X, kS, kP, kT, cos_a, parms);

% total force
F_pas = parms.Fpe_func(lmtc, parms);
F_act = FXB * parms.Fscale;
Ftot = F_act(:) + F_pas(:);

xdot = [Xd(:); parms.vmtc; Ld];

%     xi = parms.xi;
%     n = n_func(Q, parms.xi);

end