function[nd, Ld, Q0, Q1] = FULL_dynamics(ns, f, w, k1, k2, xi, Non, DRX, kS, kP, kT, cos_a, parms)
    
    % safety
    ns(ns<0) = 0;
          
    % compute moments
    Q = trapz(xi(:), [ns xi(:).*ns]);
    Q0 = Q(1);
    Q1 = Q(2);

    %% spatial integrals
    % attachment and detachment at each strain
    beta = parms.f_func(xi, f, w);
    phi = -(parms.g_func(xi, k1(1), -k1(2)) + parms.g_func(xi, k2(1), -k2(2))) .* ns';   
   
    % change in cross-bridge attachment
    nd = DRX * (beta * (Non - Q0)) + phi;
    
    % first determine contraction velocity
    Qdot = trapz(xi(:), [nd(:) xi(:).*nd(:)]);
    FXBdot_isom  = Qdot(2) + parms.ps * Qdot(1);

    % from solving velocity constraint
    Ld = Length_dynamics(FXBdot_isom, Q0, kS, kP, kT, cos_a, parms);
    
end