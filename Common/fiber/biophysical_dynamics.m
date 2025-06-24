function[Xd, Ld] = biophysical_dynamics(Ca, x, kS, kP, kT, cos_a, parms)

    % retrieve states
    Non     = x(1);
    Q       = x(2:end-1);
    DRX     = x(end);
    
    %% cross-bridge dynamics
    if length(Q) == 3 % assume approximation
        [Xdot, Ld, Q0, Q1] = DM_dynamics(Q, parms.f, parms.w, [parms.k11 parms.k12], [parms.k21 -parms.k22], parms.gaussian.IGef, Non, DRX, kS, kP, kT, cos_a, parms);
    else % original formulation
        [Xdot, Ld, Q0, Q1] = FULL_dynamics(Q, parms.f, parms.w, [parms.k11 parms.k12], [parms.k21 -parms.k22], parms.xi, Non, DRX, kS, kP, kT, cos_a, parms);
    end
    
    %% thin filament activation
    Ntot = max(parms.act * parms.Noverlap, 1e-6);
    [Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, parms.kon, parms.koff, parms.koop, Ntot);
    Nond    = Jon - Joff; % Eq (7)

    %% thick filament dynamics 
    FXB = max(Q1 + parms.ps * Q0, 0);
    [J1, J2] = ThickFilament_Dynamics(FXB, DRX, parms.J1, parms.J2, parms.JF, Ntot);
    Dd = J1 - J2; 
    
    %% combined state derivative vector
    Xd = [Nond; Xdot(:); Dd];

end