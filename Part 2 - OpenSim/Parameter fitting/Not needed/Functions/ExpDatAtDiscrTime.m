function [q_spline, qd_spline, N_spline, tvect_spline] = ExpDatAtDiscrTime(t_span, t_exp, q_exp, dt_spline)
%Interpolate data at discretized time using spline
%   1. Get time vector for simulation (tvect)
%   2. Interpolate experimental data

tvect_spline = t_span(1):dt_spline:t_span(end);     % Time vector for simulation
N_spline     = length(tvect_spline) ;                      % Number of steps in simulation

qspline              = spline(t_exp,q_exp);            % Spline fit
[q_spline,qd_spline] = SplineEval_ppuval(qspline,tvect_spline,1); % Get angles and velocities

q_spline             = q_spline';
qd_spline            = qd_spline';
tvect_spline         = tvect_spline'; 
end

