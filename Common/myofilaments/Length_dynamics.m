function[Ld] = Length_dynamics(FXBdot_isom, Q0, kS, kP, kT, cos_a, parms)

    % from solving velocity constraint
    Ld  = (cos_a * parms.vmtc * kS * kT - cos_a^2 * (FXBdot_isom*(kS+kP)) - kT * FXBdot_isom) / (kS*kT + cos_a^2 * (kS * (Q0+kP) + kP*Q0) + kT * Q0);

%     if cos_a < 1e-2
%         keyboard
%     end

%     Ld  = (parms.vmtc/cos_a * kS * kT - FXBdot_isom*(kS+kP+kT)) ...
%         / (Q0 * (kS+kP+kT) + kS * (kT + kP));

    if parms.no_tendon
        Ld = parms.c * parms.vmtc;
    end 
    
end