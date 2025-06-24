function[error] = OdeFun_FV(t,y,yp, parms)

Q0 = y(1);
Q1 = y(2);
Q2 = y(3);

dQ0dt = yp(1);
dQ1dt = yp(2);
dQ2dt = yp(3);

vMtilda = interp1(parms.ti, parms.vts, t);

if length(y) == 3
    error_thin = [];
    error_thick = [];
    Non = parms.a;
    DRX = 1;
else
    Non = y(4);
    dNondt = yp(4);
    
    DRX = y(5);
    dDRXdt = yp(5);
    
    [error_thin, ~] = ThinEquilibrium(parms.Ca, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.act * parms.Noverlap);
    
    F = Q1 + Q0 + parms.d * parms.vMtilda;
    [error_thick, ~] = ThickEquilibrium(Q0, F, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.act * parms.Noverlap);
end

% Compute p and q
Q00 = max(Q0, 1e-6);
p = Q1./Q00; 
q = Q2./Q00 - p.^2;  

Act = Non;
% DRX = DRX;
% [error_fv] = MuscleEquilibrium_alt(Q0, Q1, Q2, dQ0dt, dQ1dt, dQ2dt, parms.f, parms.k11, parms.k12, parms.k21, parms.k22, Act, vMtilda, DRX);
[error_fv] = MuscleEquilibrium(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, parms.f, parms.k11, parms.k12, parms.k21, parms.k22, Act, vMtilda, DRX);

error = [error_thin; error_thick; error_fv];

% dX = [dX_muscle(:); dX_thin(:); dX_thick(:)];

end