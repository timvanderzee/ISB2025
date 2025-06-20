function [dNondt] = ThinEquilibrium_exp(Act, Q0, Non, kon, koff, koop)

    Noverlap = 1;

    % quantities 
    Noff = Noverlap - Non;
    Noff(Noff<0) = 0; % could overshoot due to fast dynamics
    
    % thin filament dynamics    
    Jon     = Act * kon * Noff * (1 + koop * (Non/Noverlap)); % Eq (1)
    Joff    = koff * (Non - Q0) *     (1 + koop * Noff/Noverlap); % Eq (2)
    
    dNondt = (Jon - Joff); % Eq (7)

end