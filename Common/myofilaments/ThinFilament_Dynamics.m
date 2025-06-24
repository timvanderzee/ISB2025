function[Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, kon, koff, koop, Ntot)

    % quantities 
%     Noff = max(Noverlap - Non, 0);
    Noff = Ntot - Non;
    Jon     = Ca * kon * Noff  .* (1 + koop * (Non/Ntot)); % Eq (1)
    Joff    = koff * (Non - Q0) .*     (1 + koop * Noff/Ntot); % Eq (2)
    
    %     Noff(Noff<0) = 0; % could overshoot due to fast dynamics
    
    % thin filament dynamics    
%     Jon     = Act * kon * Noff  .* (1 + max(koop * (Non/Noverlap),0)); % Eq (1)
%     Joff    = koff * (Non - Q0) .*     (1 + max(koop * Noff/Noverlap, 0)); % Eq (2)
    

    
end