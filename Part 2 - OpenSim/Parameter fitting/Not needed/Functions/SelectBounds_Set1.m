function [Ub, Lb] = SelectBounds_Set1(data_exp)
%Function to make bounds
% Change bounds here

% States
Lb.x    =  -pi;     Ub.x = 10*pi/180;
Lb.a    = 0.001;    Ub.a = 1;
Lb.kFpe = 0;        Ub.kFpe = 0.2;
Lb.B    = 0.001;   Ub.B = 0.15;
Lb.kR   = 1e-4;     Ub.kR = 5;
Lb.kY   = 1e-4;     Ub.kY = 10; 
Lb.tres = 1e-4;     Ub.tres = 1;

Lb.vMtilda= -10;    Ub.vMtilda = 10;
Lb.lMtilda = 0.2;   Ub.lMtilda = 1.5;
Lb.lM_projected = 1e-4;

Lb.dt1 = 0.001;
Ub.dt1 = 0.01;

bound_mean    = mean(data_exp.qspline(end-250:end));
bound_std     = std(data_exp.qspline(end-250:end));
Lb.x_end = bound_mean-2*bound_std;
Ub.x_end = bound_mean+2*bound_std;
end

