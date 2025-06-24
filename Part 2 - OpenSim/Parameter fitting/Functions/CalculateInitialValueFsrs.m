function [Fsrs_f1_cal] = CalculateInitialValueFsrs(lMtilda, a_ext, info, params_OS, data_exp)
% FMltilda
% Active muscle force-length characteristics
b11 = params_OS.Fa(1); b21 = params_OS.Fa(2); b31 = params_OS.Fa(3); b41 = params_OS.Fa(4);
b12 = params_OS.Fa(5); b22 = params_OS.Fa(6); b32 = params_OS.Fa(7); b42 = params_OS.Fa(8);
b13 = 0.1;             b23 = 1;               b33 = 0.5*sqrt(0.5);   b43 = 0;

num3     = lMtilda-b23;
den3     = b33+b43*lMtilda;
FMtilda3 = b13*exp(-0.5*num3.^2./den3.^2);

num1     = lMtilda-b21;
den1     = b31+b41*lMtilda;
FMtilda1 = b11*exp(-0.5*num1.^2./den1.^2);

num2     = lMtilda-b22;
den2     = b32+b42*lMtilda;
FMtilda2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilda = FMtilda1+FMtilda2+FMtilda3;

% kSRS & N_1
kSRS         = info.kSRS; 
N_1          = data_exp.N_1; 

% Stretch
stretch      = lMtilda(1,N_1)- lMtilda(1,1);

% Intial value 
tan_val_incr = (0.5* tanh(1000*(5.7e-3-stretch))+0.5).*FMltilda(1,N_1)*a_ext*kSRS.*stretch; 
tan_val_plat = (0.5* tanh(1000*(stretch-5.7e-3))+0.5).*FMltilda(1,N_1)*a_ext*kSRS*5.7e-3; 

Fsrs_f1_cal  = tan_val_incr + tan_val_plat; 
end

