function [error, dX] = ThinEquilibrium(Ca, Q0, Non, dNondt, kon, koff, koop, Ntot)

    [Jon, Joff] = ThinFilament_Dynamics(Ca, Q0, Non, kon, koff, koop, Ntot);
    
%     dX = max(Jon,0) - max(Joff,0);
    dX = Jon - Joff;
    
    error   = dNondt - dX; % Eq (7)

end