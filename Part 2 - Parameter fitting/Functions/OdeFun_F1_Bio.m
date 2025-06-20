function [error_f1] = OdeFun_F1_Bio(t,y,yp, params)
%Odefunction - forward simulation of F1 of pendulum simulaties
%   t = time
%   y = State
%   yp = State derivative

% States
%  q (y1), qd (y2), lMtilda_ext (y3), lMtilda_flex (y4), Fsrs_del (y5), Fce_del (y6)

% State derivatives
% qd (yp1), qdd (yp2), vMtilda_ext (yp3), vMtilda_flex (yp4), dFsrs_deldt (yp5),
% dFce_deldt (yp6)

% Params
% lMtilda_init, lM_projected, kFpe, a_ext, a_flex, data_exp, coeff_LMT_ma,
% params_OS, shift, B, info, kR, I_opt

% lM projected

% w 
lMo    = params.params_OS.MT(2,:); 
alphao = params.params_OS.MT(4,:); 
w_ext  = lMo(1).* sin(alphao(1));
w_flex = lMo(2).* sin(alphao(2));
w      = [w_ext; w_flex];  

lM_ext  = y(3).* lMo(1); % voor de extensor  
lM_flex = y(4).* lMo(2); % voor de flexor 
lM      = [lM_ext; lM_flex]; 

lM_projected = sqrt(lM.^2 - w.^2);

vMtilda = yp(3:4);

error_all = CalculateBioMusculoSkeletalDynamics_F1_ForceFeedback(y(1),y(2),yp(2), y(3:4), lM_projected,params.kFpe,vMtilda, ...
    params.a_ext, params.a_flex, params.data_exp, params.coeff_LMT_ma, params.params_OS, params.shift, params.B, params.info, params.I_opt, y(5:6), y(7:8), yp(5:6), yp(7:8));

error_vel = y(2) - yp(1);
error_f1 = [error_vel; error_all(9); error_all(1:6)];

end