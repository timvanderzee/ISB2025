function [error, dX] = ThinEquilibrium(Act, Q0, Non, dNondt, kon, koff, koop, Noverlap)

    [Jon, Joff] = ThinFilament_Dynamics(Act, Q0, Non, kon, koff, koop, Noverlap);
    
%     dX = max(Jon,0) - max(Joff,0);
    dX = Jon - Joff;
    
    error   = dNondt - dX; % Eq (7)

end