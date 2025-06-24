function[Fpe, kP] = get_parallel_force(lMtilda, Fpparam)

% Calculate passive force length
% Fpe = musclepassiveforcelength(lMtilda);
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilda - 0.10e1) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);
kP = kpe / (e0*Fpparam(2)) * exp(kpe * (lMtilda - 1) / e0);

end