function [error] = CalculateMusculoSkeletalDynamics_F2_ForceFeedback(q,qd,qdd,lMtilda, lM_projected, kFpe, vMtilda, a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, Fsrs, dFsrsdt, Fsrs_del, dFsrs_deldt, kR,Fce_del, dFce_deldt, I_opt)
% Calculate Muscle tendon lengths and moment arms 
[lMT, MA] = CalculateMuscleTendonLengthAndMomentArms(q, data_exp, coeff_LMT_ma); 

% Calculate Tendon Force, muscle length and tendon length
[FT, lM, lT, fse, w ] = CalculateTendonForce(lMtilda, lM_projected, params_OS, lMT, shift);

% Get ForceLengthVelocity Relationships
[Fpe, FMltilda, FMvtilda] = getForceLengthVelocityRelation(lMtilda, kFpe, params_OS, vMtilda);

% Fsrs 
dFsrsdt_cal = -Fsrs/0.05;  % Calculated value of derivative of Fsrs (exponential decay)
error_Fsrs  = dFsrsdt - dFsrsdt_cal; 

% Reflexes
tau = info.tau;

tres   = a_ext*FMltilda(1,1)*FMvtilda(1,1) + Fsrs(1); 
a_refl = kR * (Fce_del-tres); 
a      = a_ext + (0.5*tanh(1000*a_refl)+0.5).*a_refl;

% FMce 
Fce_ext  = a.*FMltilda(1,:).* FMvtilda(1,:) + Fsrs;
%(a_ext + kR * Fsrs_del + kY * dFsrs_deldt - tres).* FMltilda(1,:).* FMvtilda(1,:) + Fsrs;      % FMce = fse.* lM ./(lMT-lT) - Fpe;  
Fce_flex = a_flex.* FMltilda(2,:).* FMvtilda(2,:); 
Fce      = [Fce_ext; Fce_flex];

dFsrs_deldt_cal = (Fsrs - Fsrs_del)./tau; 
error_Fsrs_del  = dFsrs_deldt - dFsrs_deldt_cal; 

dFce_deldt_cal  = ((a_ext.* (FMltilda(1,:).* FMvtilda(1,:))  - Fce_del))./tau;
error_Fce_del   = dFce_deldt - dFce_deldt_cal; 

% Force equilibrium 
[error_force] = ForceEquilibrium(Fce, Fpe, lMT, lT, lM, fse); 

% Error lm and lM projected
error_musclegeometry_ext  = lM(1,:).^2 - w(1,:).^2 - lM_projected(1,:).^2; 
error_musclegeometry_flex = lM(2,:).^2 - w(2,:).^2 - lM_projected(2,:).^2;
error_musclegeometry      = [error_musclegeometry_ext error_musclegeometry_flex]; 

%Implicit skelet dynamics 
I = I_opt; 
m = params_OS.inert.mass_OS; 
l = params_OS.inert.lc_OS; 
g = 9.81; 
ma_ext = MA(1,:); 
ma_flex = MA(2,:); 

error_skelet = qdd*I + m*g*l*cos(q) - FT(1,:).*ma_ext - FT(2,:).*ma_flex + B*qd ;  

% Total error
error = [error_force error_musclegeometry error_skelet error_Fsrs error_Fsrs_del error_Fce_del ];
end

