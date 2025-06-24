function[xdot, Ftot, Ld] = fiber_dynamics_old(t, x, parms, Ca)

% retrieve states
Non     = x(1);
Q0      = x(2);
Q1      = x(3);
Q2      = x(4);
DRX     = x(5);
lmtc    = x(6);

% make sure values are within range
eps = 1e-6;
Non = max(Non, eps);
Q0 = max(Q0,eps);
% Q1 = max(Q1, -Q0);
Q2 = max(Q2, eps.^2);
FXB = max(Q0 * parms.ps + Q1, 0);

% get geometry
dLS = parms.Lse_func(FXB, parms);
kS = parms.kse_func(dLS, parms);
kT = 1000;
kP = parms.kpe;
cos_a = 1;

% simulate biophysical dynamics
X = [Non Q0 Q1 Q2 DRX];
[Xd, Ld] = biophysical_dynamics_old(Ca, X, kS, kP, kT, cos_a, parms);

% total force
F_pas = parms.Fpe_func(lmtc, parms);
F_act = FXB * parms.Fscale;
Ftot = F_act(:) + F_pas(:);

xdot = [Xd(:); parms.vmtc; Ld];

%     xi = parms.xi;
%     n = n_func(Q, parms.xi);

end