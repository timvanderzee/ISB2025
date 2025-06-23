function[auxdata] = get_muscle_parms()

% PARAMETERS
load('Fvparam.mat','Fvparam')
Fvparam(1) = 1.475*Fvparam(1);
Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3) = Fvparam(3) + 0.75;
Fvparam(4) = Fvparam(4) - 0.027;
auxdata.Fvparam = Fvparam;

load('Faparam.mat','Faparam')
auxdata.Faparam = Faparam;

e0 = 0.6;
kpe = 4;
t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1);
t7 = exp(kpe);
pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];
auxdata.Fpparam = Fpparam;

end