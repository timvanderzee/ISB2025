function[Ld] = Length_dynamics(FXBdot_isom, Q0, kS, kP, kT, cos_a, parms)

    % from solving velocity constraint
%     Ld  = (cos_a * parms.vmtc * kS * kT - cos_a^2 * (FXBdot_isom*(kS+kP)) - kT * FXBdot_isom) / (kS*kT + cos_a^2 * (kS * (Q0+kP) + kP*Q0) + kT * Q0);

%     if cos_a < 1e-2
%         keyboard
%     end

dL = parms.dL;
% 
% a = 1 - (kP+Q0+kP*Q0)/kT - ((kP+Q0+kP*Q0) * -dL)/(kT*FXBdot_isom/kS*cos_a) - Q0/kS;
% b = parms.vmtc / cos_a - ((kP+Q0+kP*Q0) * -dL)/(kT*(1+Q0/kS)*cos_a) ...
%     -FXBdot_isom * (1 + kP)/kT ...
%     - ((1+kP)*-dL)/(kT/kS * cos_a) - FXBdot_isom/kS;
% c = -FXBdot_isom * (1 + kP)/kT * -dL/((1+Q0/kS)*cos_a);
% 
% Ld = -b + sqrt(b.^2 - 4 * a * c) ./ (2 * a); 

    Ld  = (parms.vmtc/cos_a * kS * kT - FXBdot_isom*(kS+kP+kT)) ...
        / (Q0 * (kS+kP+kT) + kS * (kT + kP));

    if parms.no_tendon
        Ld = parms.c * parms.vmtc;
    end 
    
end