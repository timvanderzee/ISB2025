function[dydt] = OdeFun_FV_exp(t,y, parms)

Q0 = y(1);
Q1 = y(2);
Q2 = y(3);

if length(y)>3
    Non = y(4);
    [dNondt] = ThinEquilibrium_exp(parms.Ca, Q0, Non, parms.kon, parms.koff, parms.koop);
else
    dNondt = [];
    Non = 1;
end

[dQ0dt, dQ1dt, dQ2dt] = MuscleEquilibrium_alt_exp(Q0, Q1, Q2, parms.f, parms.k11, parms.k12, parms.k21, parms.k22, Non, parms.vMtilda);

dydt = [dQ0dt; dQ1dt; dQ2dt; dNondt];

end