function[Xd, Ld] = biophysical_dynamics(Ca, x, kS, kP, kT, cos_a, parms)

    % retrieve states
    Non     = x(1);
    Q0      = x(2);
    Q1      = x(3);
    Q2      = x(4);
    DRX     = x(5);
    
    FXB = max(Q1 + parms.ps * Q0, 0);
    
    %% cross-bridge dynamics
    % get mean and standard deviation
    eps = 1e-8;
    p = Q1/max(Q0, eps); % mean of the distribution
    q = max(Q2/max(Q0, eps) - p^2, eps); % standard deviation of the distribution

    % spatial integrals
    beta = [1 0 parms.w^2] * parms.f;

    phi(1) = -parms.gaussian.IGef{1}([Q0 p 2*q],[parms.k11 parms.k12]) -parms.gaussian.IGef{1}([Q0 p 2*q],[parms.k21 -parms.k22]); 
    phi(2) = -parms.gaussian.IGef{2}([Q0 p 2*q],[parms.k11 parms.k12]) -parms.gaussian.IGef{2}([Q0 p 2*q],[parms.k21 -parms.k22]); 
    phi(3) = -parms.gaussian.IGef{3}([Q0 p 2*q],[parms.k11 parms.k12]) -parms.gaussian.IGef{3}([Q0 p 2*q],[parms.k21 -parms.k22]); 
  
    % isometric changes in distribution (i.e. velocity-independent)
    Q0dot_isom = DRX .* (beta(1) .* (Non-Q0)) + phi(1);
    Q1dot_isom = DRX .* (beta(2) .* (Non-Q0)) + phi(2);
    Q2dot_isom = DRX .* (beta(3) .* (Non-Q0)) + phi(3);
    
    % isometric change in force
    FXBdot_isom  = Q1dot_isom + parms.ps * Q0dot_isom;

    % from solving velocity constraint
    Ld  = (cos_a * parms.vmtc * kS * kT - cos_a^2 * (FXBdot_isom*(kS+kP)) - kT * FXBdot_isom) / (kS*kT + cos_a^2 * (kS * (Q0+kP) + kP*Q0) + kT * Q0);
%             Ld  = (parms.vmtc * kT - FXBdot_isom.*cos_a) / (Q0 + kT./cos_a + kP);    

    if parms.no_tendon
        Ld = parms.c * parms.vmtc;
    end
    
    % change in distribution    
    Q0dot = Q0dot_isom;
    Q1dot = Q1dot_isom + 1 * Ld * Q0;
    Q2dot = Q2dot_isom + 2 * Ld * Q1;
    
    %% thin filament activation    
    Ntot = max(parms.act * parms.Noverlap, eps);
    Noff = max(Ntot - Non, 0);  
    
    Jon     = Ca * parms.kon * Noff * (1 + parms.koop * (Non/Ntot)); % Eq (1)
    Joff    = parms.koff * (Non - Q0) *     (1 + parms.koop * Noff/Ntot); % Eq (2)
    Nond    = Jon - Joff; % Eq (7)

    %% thick filament dynamics
    SRX = 1 - DRX;
    J1 = parms.J1 * (1 + parms.JF * max(FXB,0)/Ntot) .* SRX;
    J2 = parms.J2 .* DRX;
    Dd = J1 - J2; 
    
    %% combined state derivative vector
    Xd = [Nond; Q0dot; Q1dot; Q2dot; Dd; parms.vmtc; Ld];

    % remove length if needed
    Xd = Xd(1:length(x));

end