function[vM] = contractile_dynamics(a, FMltilda, FMce, Fvparam, params)

vMmax = params(5); % in m/s

FMvtilda = FMce./(a.*FMltilda);

% Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

vMtilda = 1/e2*(sinh((FMvtilda-e4)/e1)-e3); % 0-1

vM = vMtilda.*vMmax;

end