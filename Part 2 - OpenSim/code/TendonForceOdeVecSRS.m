function [dx, FT] = TendonForceOdeVecSRS(t,x,t_input,A,LMT,VMT,mparams, Fvparam, Fpparam, Faparam, lMtilda_isom, ksrs, kT)

% Input
a = interp1(t_input, A, t);
lMT = interp1(t_input, LMT, t);
vMT = interp1(t_input, VMT, t);

% Length and overlap
fse = x(1);
[~, lMtilda, cos_alpha] = get_lM_from_fse(fse, lMT, mparams, kT);
FMltilda = get_overlap(lMtilda, Faparam);

% Force from short-range stiffness
dlM = (0.5*tanh(1000*(lMtilda - lMtilda_isom))+0.5) .*(lMtilda - lMtilda_isom);
dlM = (0.5*tanh(1000*(-dlM+5.7*10^(-3)))+ 0.5).*dlM + (0.5*tanh(1000*(dlM-5.7*10^(-3)))+0.5)*5.7*10^(-3);
Fsrs = ksrs * a .* FMltilda.*dlM;

% Parallel force
Fpe = get_parallel_force(lMtilda, Fpparam);

% Contractile element force
FMce = max(fse./cos_alpha-Fpe-Fsrs, 0);

% Contractile dynamics
vM = contractile_dynamics(a, FMltilda, FMce, Fvparam, mparams);

% Tendon velocity and force
FMo = mparams(1);
vT = vMT-vM./cos_alpha;
FT = fse .* FMo;

% Force-rate
lTs = mparams(3);
dfse = kT.*vT./lTs;

% State derivative
dx(1,1) = dfse;



return