function[xdot, Ftot, Ld] = fiber_dynamics(t, x, parms, Ca)

% retrieve states
Non     = x(1); % thin filment activation
Q       = x(2:end-3); % cross-bridge states (moments or number of XBs in each bin)
DRX     = x(end-2); % DRX state
lmtc    = x(end-1); % CE + SE length
L       = x(end-0); % CE length

% get moments and strain vector
[Q0, Q1, parms.xi] = force_from_distribution(Q, L, parms);

% Q0: number of cross-bridges
% Q1: (delta) force
% Q2: (delta) elastic strain energy
% because xi is centered around 0, Q1 and Q2 are actually delta terms wrt the powerstroke

eps = 1e-6;
Q0 = max(Q0, eps); % prevent division by zero
FXB = parms.ps * Q0 + Q1; % assuming normalized power stroke = 1

% make sure values are within range
Non = max(Non, eps); % prevent division by 0 

% get geometry
dLS = parms.Lse_func(FXB, parms); % SE length-force
kS = parms.kse_func(dLS, parms); % in-series with XBs
kT = 1000; % high stiffness, because no tendon
kP = parms.kpe; % parallel stiffness
cos_a = 1;

% force-length (need to turn off for force-velocity and/or part 1?)
half_s_len_norm = parms.s/2/parms.h;
curv_w = 0.5;
overlap_func = @(L, parms) max(0.00001, ...
    -(L-half_s_len_norm*curv_w)*(L+half_s_len_norm*curv_w)...
    /(half_s_len_norm*curv_w)^2);

parms.Noverlap = overlap_func(L, parms);

% simulate biophysical dynamics
X = [Non; Q(:); DRX];
 
% compute cross-bridge state derivative and cross-bridge velocity
[Xd, Ld] = biophysical_dynamics(Ca, X, kS, kP, kT, cos_a, parms);

% total force
F_pas = parms.Fpe_func(lmtc, parms);
F_act = FXB * parms.Fscale;
Ftot = F_act(:) + F_pas(:);

xdot = [Xd(:); parms.vmtc; Ld];

%     xi = parms.xi;
%     n = n_func(Q, parms.xi);

end