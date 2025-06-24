function[Xd, Ld] = biophysical_dynamics(Ca, x, kS, kP, kT, cos_a, parms)

    % retrieve states
    Non     = x(1);
    Q0      = x(2);
    Q1      = x(3);
    Q2      = x(4);
    DRX     = x(5);
    
    %% cross-bridge dynamics
    % get mean and standard deviation
    eps = 1e-8;
    p = Q1/max(Q0, eps); % mean of the distribution
    q = max(Q2/max(Q0, eps) - p^2, eps); % standard deviation of the distribution

    % Cross-bridge dynamics
    [Q0dot_isom, Q1dot_isom, Q2dot_isom] = CrossBridge_Dynamics(Q0, p, q, parms.f, parms.w, [parms.k11 parms.k12], [parms.k21 -parms.k22], parms.gaussian.IGef, Non, DRX);
    
    % isometric change in force
    FXBdot_isom  = Q1dot_isom + parms.ps * Q0dot_isom;

    % from solving velocity constraint
    Ld  = (cos_a * parms.vmtc * kS * kT - cos_a^2 * (FXBdot_isom*(kS+kP)) - kT * FXBdot_isom) / (kS*kT + cos_a^2 * (kS * (Q0+kP) + kP*Q0) + kT * Q0);

    if parms.no_tendon
        Ld = parms.c * parms.vmtc;
    end
    
    % change in distribution    
    Q0dot = Q0dot_isom;
    Q1dot = Q1dot_isom + 1 * Ld * Q0;
    Q2dot = Q2dot_isom + 2 * Ld * Q1;
    
    %% thin filament activation
    Ntot = max(parms.act * parms.Noverlap, 1e-6);
    [Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, parms.kon, parms.koff, parms.koop, Ntot);
    Nond    = Jon - Joff; % Eq (7)

    %% thick filament dynamics 
    FXB = max(Q1 + parms.ps * Q0, 0);
    [J1, J2] = ThickFilament_Dynamics(FXB, DRX, parms.J1, parms.J2, parms.JF, Ntot);
    Dd = J1 - J2; 
    
    %% combined state derivative vector
    Xd = [Nond; Q0dot; Q1dot; Q2dot; Dd; parms.vmtc; Ld];

    % remove length if needed
    Xd = Xd(1:length(x));

end