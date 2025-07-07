function[Xdot, Ld, Q0, Q1] = DM_dynamics(Q, f, w, k1, k2, IGef, Non, DRX, kS, kP, kT, cos_a, parms)
    
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);

    % get mean and standard deviation
    eps = 1e-6;
    p = Q1/max(Q0, eps); % mean of the distribution
    q = max(Q2/max(Q0, eps) - p^2, eps); % standard deviation of the distribution

    % Cross-bridge dynamics
    [Q0dot_isom, Q1dot_isom, Q2dot_isom] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX);
    
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
    
    Xdot = [Q0dot; Q1dot; Q2dot];
    
end