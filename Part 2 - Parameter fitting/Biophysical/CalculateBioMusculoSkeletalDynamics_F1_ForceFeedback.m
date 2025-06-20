function [error] = CalculateBioMusculoSkeletalDynamics_F1_ForceFeedback(q,qd,qdd,lMtilda, lM_projected, kFpe, vMtilda, a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, I_opt, Q0, Q2, dQ0dt, dQ2dt)

% Calculate Muscle tendon lengths and moment arms 
[lMT, MA] = CalculateMuscleTendonLengthAndMomentArms(q, data_exp, coeff_LMT_ma); 

% calc vMT
lMo    = params_OS.MT(2,:); 
vMTtilda_ext = (MA(1,:) .* qd)./lMo(1);
vMTtilda_flex = (MA(2,:) .* qd)./lMo(2);
vMTtilda = [vMTtilda_ext; vMTtilda_flex];

% Calculate Tendon Force, muscle length and tendon length
[FT, lM, lT, fse, w, kse] = CalculateTendonForce(lMtilda, lM_projected, params_OS, lMT, shift);

% Get ForceLengthVelocity Relationships
[Fpe, ~, ~, kpe] = getForceLengthVelocityRelation(lMtilda, kFpe, params_OS, vMtilda);

cos_alpha_ext   = (lMT(1,:)-lT(1,:))./lM(1,:);
cos_alpha_flex   = (lMT(2,:)-lT(2,:))./lM(2,:);

Fce_ext = fse(1,:) .* cos_alpha_ext - Fpe(1,:);
Fce_flex = fse(2,:) .* cos_alpha_flex - Fpe(2,:);

% Force equilibrium 
FXB = [Fce_ext; Fce_flex];
cos_alpha = [cos_alpha_ext; cos_alpha_flex];

[error_forcerate, error_Q0, error_Q2] = ForceRateEquilibrium(a_ext, a_flex, FXB, cos_alpha, vMTtilda, kse, kpe, vMtilda, Q0, Q2, dQ0dt, dQ2dt); 

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
error = [error_forcerate error_Q0 error_Q2 error_musclegeometry error_skelet]';
end

