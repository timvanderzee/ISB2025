function J = DefineObjectiveFunction_time_Lars(x,xd,data_exp, info, vMtilda, a_ext, a_flex, kR, dt1, N_1, dt2, N)
%Define objective function that should be minimized 
w_q  = info.wq;
w_qd = info.wqd;

% We need to interpolate the experimental data to the time vector of the
% optimization
% Original data 
exp_time = data_exp.tspline-data_exp.tspline(1); 
exp_data = data_exp.qspline; 
exp_vel  = data_exp.qdspline; 

% Fase 1 
Qref_new  = interp1(linspace(1,N_1,length(exp_data(1:N_1))),exp_data(1:N_1),(1:1:N_1),'spline','extrap');
Qdref_new = interp1(linspace(1,N_1,length(exp_vel(1:N_1))),exp_vel(1:N_1),(1:1:N_1),'spline','extrap');

q_error_f1 = sumsqr(x(1:N_1) - Qref_new); 
qd_error_f1= sumsqr(xd(1:N_1) - Qdref_new); 

% Fase 2
time_exp_f2 = exp_time(N_1+1:end)-exp_time(N_1+1); 
time_opt_f2 = 0:0.002:time_exp_f2(end); 

interpolated_angles_f2 = interp1(time_exp_f2, exp_data(N_1+1:N),time_opt_f2);
interpolated_vel_f2    = interp1(time_exp_f2, exp_vel(N_1+1:N), time_opt_f2);

q_error_f2  = sumsqr(x(N_1+1:N) - interpolated_angles_f2);
qd_error_f2  = sumsqr(xd(N_1+1:N) - interpolated_vel_f2);

J = w_q * (q_error_f1+q_error_f2) + w_qd * (qd_error_f1+qd_error_f2) + 0.001 * (sumsqr(vMtilda)) + 0.1*(a_ext^2 + a_flex^2+ kR^2); 



end

