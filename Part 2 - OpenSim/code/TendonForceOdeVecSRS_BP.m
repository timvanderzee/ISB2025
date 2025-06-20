function [dx, FT] = TendonForceOdeVecSRS_BP(t,x,t_input,A,LMT,VMT,mparams, Fvparam, Fpparam, Faparam, lMtilda_isom, ksrs, kT)

% Input
a = interp1(t_input, A, t);
lMT = interp1(t_input, LMT, t); % [m]
vMT = interp1(t_input, VMT, t); % [m/s]

% Length and overlap
fse = x(1);
[~, lMtilda, cos_alpha] = get_lM_from_fse(fse, lMT, mparams, kT);
FMltilda = get_overlap(lMtilda, Faparam);

% Force from short-range stiffness
dlM = (0.5*tanh(1000*(lMtilda - lMtilda_isom))+0.5) .*(lMtilda - lMtilda_isom);
dlM = (0.5*tanh(1000*(-dlM+5.7*10^(-3)))+ 0.5).*dlM + (0.5*tanh(1000*(dlM-5.7*10^(-3)))+0.5)*5.7*10^(-3);
Fsrs = ksrs * a .* FMltilda.*dlM;

% Parallel force and stiffness
[Fpe, kP] = get_parallel_force(lMtilda, Fpparam);

% Contractile element force
FMce = max(fse./cos_alpha-Fpe-Fsrs, 0);

% Contractile dynamics
[Qd, vM, Nond, DRXd] = contractile_dynamics_BP_v3(a, FMltilda, [FMce x(2) x(3) x(4) x(5) x(6)], vMT, kT, kP, cos_alpha, mparams);
Ldd = 0;

% Tendon velocity and force
FMo = mparams(1);
vT = vMT-vM./cos_alpha;
FT = fse .* FMo;

% Force-rate
lTs = mparams(3);
dfse = kT.*vT./lTs;

% State derivative
dx(1,1) = dfse;
dx(2,1) = Qd(1);
dx(3,1) = Qd(3);
dx(4,1) = Ldd;
dx(5,1) = Nond;
dx(6,1) = DRXd;

return