function [FT, lM, lT, fse, w, kse] = CalculateTendonForce(lMtilda, lM_projected, params_OS, lMT, shift)

%Calculate Tendon Force,muscle fiber length and tendon length
lMtilda_ext       = lMtilda(1,:);
lMtilda_flex      = lMtilda(2,:);
lM_projected_ext  = lM_projected(1,:);
lM_projected_flex = lM_projected(2,:); 

% w 
lMo    = params_OS.MT(2,:); 
alphao = params_OS.MT(4,:); 
w_ext  = lMo(1).* sin(alphao(1));
w_flex = lMo(2).* sin(alphao(2));
w      = [w_ext; w_flex];  


% lM (muscle fiber length) 
lM_ext  = lMtilda_ext.*  lMo(1); % voor de extensor  
lM_flex = lMtilda_flex.* lMo(2); % voor de flexor 
lM      = [lM_ext; lM_flex]; 

% lT (tendon length)
lMT_ext  = lMT(1,:);
lMT_flex = lMT(2,:);
lT_ext   = lMT_ext - lM_projected_ext;
lT_flex  = lMT_flex - lM_projected_flex;
lT       = [lT_ext; lT_flex];  

% lTs (Tendon slack length)
lTs    = params_OS.MT(3,:); 

% lTtilda
lTtilda_ext  = lT_ext./lTs(1); 
lTtilda_flex = lT_flex./lTs(2); 

% Fse
fse_ext      = (exp(35*(lTtilda_ext - 0.995)))/5-0.25 + shift;
fse_flex     = (exp(35*(lTtilda_flex - 0.995)))/5-0.25 + shift;
fse          = [fse_ext; fse_flex]; 

% kse
kse_ext = 35/5 * (exp(35*(lTtilda_ext - 0.995)));
kse_flex = 35/5 * (exp(35*(lTtilda_flex - 0.995)));

lo = params_OS.MT(2,:); 
Ts = params_OS.MT(3,:); 

% express kse wrt to lMopt
kse = [kse_ext; kse_flex] .* lo(:)./Ts(:);
%%

% FMo
FMo          = params_OS.MT(1,:);

% Compute tendon force
FT_ext       = FMo(1).* fse_ext; 
FT_flex      = FMo(2).* fse_flex; 
FT           = [FT_ext; FT_flex]; 
end

