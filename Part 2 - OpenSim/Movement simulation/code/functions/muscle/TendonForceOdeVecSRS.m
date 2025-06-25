function [dx, FT] = TendonForceOdeVecSRS(t,x,t_input,A,LMT,VMT,mparams, Fvparam, Fpparam, Faparam, kT, type)

% Input
a = interp1(t_input, A, t);
lMT = interp1(t_input, LMT, t);
vMT = interp1(t_input, VMT, t);

% Length and overlap
fse = x(1);
[~, lMtilda, cos_alpha] = get_lM_from_fse(fse, lMT, mparams, kT);
FMltilda = get_overlap(lMtilda, Faparam);

% Parallel force
[Fpe, kP] = get_parallel_force(lMtilda, Fpparam);

% Contractile element force
FMce = max(fse./cos_alpha-Fpe, 0);

% Contractile dynamics
if strcmp(type, 'Hill')
    [vM, Xd] = contractile_dynamics(a, FMltilda, FMce, Fvparam, mparams);
elseif strcmp(type, 'Biophysical')
    [vM, Xd] = contractile_dynamics_BP(a, FMltilda, [FMce x(2) x(3) x(4) x(5)], vMT, kT, kP, cos_alpha, mparams);
end

% Tendon velocity and force
FMo = mparams(1);
vT = vMT-vM./cos_alpha;
FT = fse .* FMo;

% Force-rate
lTs = mparams(3);
dfse = kT.*vT./lTs;

% State derivative
dx = [dfse(:); Xd(:)];

return