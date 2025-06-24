function [error_force] = ForceEquilibrium(Fce, Fpe, lMT, lT, lM, fse)
%Force Equilibrium

FM_ext  = Fce(1,:) + Fpe(1,:);
FM_flex = Fce(2,:) + Fpe(2,:);

% Force equilibrium extensor
cos_alpha_ext   = (lMT(1,:)-lT(1,:))./lM(1,:);
error_force_ext = FM_ext.*cos_alpha_ext - fse(1,:); 

% Force equilibrium flexor
cos_alpha_flex    = (lMT(2,:)-lT(2,:))./lM(2,:);
error_force_flex  = FM_flex.*cos_alpha_flex - fse(2,:);

% Total error force
error_force = [error_force_ext error_force_flex]; 

end

