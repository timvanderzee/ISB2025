function [error_forcerate, error_Q0, error_Q2] = ForceRateEquilibrium(a_ext, a_flex, FXB, cos_alpha, vMTtilda, kse, kpe, vMtilda, Q0, Q2, dQ0dt, dQ2dt)

% Extensor
[error_forcerate_ext, error_Q0_ext, error_Q2_ext]       = MuscleEquilibrium(a_ext,  FXB(1,:), cos_alpha(1,:), vMTtilda(1,:), kse(1,:), kpe(1,:), vMtilda(1,:), Q0(1,:), Q2(1,:), dQ0dt(1,:), dQ2dt(1,:));

% Flexor
[error_forcerate_flex, error_Q0_flex, error_Q2_flex]    = MuscleEquilibrium(a_flex, FXB(2,:), cos_alpha(2,:), vMTtilda(2,:), kse(2,:), kpe(2,:), vMtilda(2,:), Q0(2,:), Q2(2,:), dQ0dt(2,:), dQ2dt(2,:));

error_forcerate = [error_forcerate_ext error_forcerate_flex]; 
error_Q0 = [error_Q0_ext error_Q0_flex]; 
error_Q2 = [error_Q2_ext error_Q2_flex]; 

end

