function [error, dX] = ThinEquilibrium(Act, Q0, Non, dNondt, kon, koff, koop, Noverlap)

%     Noverlap = 1;

    % quantities 
%     Noff = max(Noverlap - Non, 0);
    Noff = Noverlap - Non;
%     Noff(Noff<0) = 0; % could overshoot due to fast dynamics
    
    % thin filament dynamics    
%     Jon     = Act * kon * Noff  .* (1 + max(koop * (Non/Noverlap),0)); % Eq (1)
%     Joff    = koff * (Non - Q0) .*     (1 + max(koop * Noff/Noverlap, 0)); % Eq (2)
    
    Jon     = Act * kon * Noff  .* (1 + koop * (Non/Noverlap)); % Eq (1)
    Joff    = koff * (Non - Q0) .*     (1 + koop * Noff/Noverlap); % Eq (2)
    
%     dX = max(Jon,0) - max(Joff,0);
    dX = Jon - Joff;
    
    error   = dNondt - dX; % Eq (7)

end