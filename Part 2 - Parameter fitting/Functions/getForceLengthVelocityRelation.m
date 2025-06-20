function [Fpe, FMltilda, FMvtilda, kpe] = getForceLengthVelocityRelation(lMtilda, kFpe, params_OS, vMtilda)
%calculate Fpe, FMltilda, FMvtilda
%   Fpe      = Normalized passive muscle force 
%   FMltilda = Normalized force-length multiplier
%   FMvtilda = Normalized force-velocity multiplier

lMtilda_ext = lMtilda(1,:); 
lMtilda_flex= lMtilda(2,:); 
vMtilda_ext = vMtilda(1,:); 
vMtilda_flex= vMtilda(2,:);

% Fpe
% ext
e0          = 0.6;   
kp = 4;
t5_ext      = exp(kp * (lMtilda_ext - kFpe*10) / e0);        
Fpe_ext     = ((t5_ext - 0.10e1) - params_OS.Fp(1)) / params_OS.Fp(2); %Fpe = musclepassiveforcelength(lMtilda);
% Flexor
e0          = 0.6;   
kp = 4;
t5_flex     = exp(kp * (lMtilda_flex - kFpe*10) / e0);        
Fpe_flex    = ((t5_flex - 0.10e1) - params_OS.Fp(1)) / params_OS.Fp(2); %Fpe = musclepassiveforcelength(lMtilda);

Fpe         = [Fpe_ext; Fpe_flex];  

% kpe
kpe_ext = kp / (e0*params_OS.Fp(2)) * exp(kp * (lMtilda_ext - kFpe*10) / e0);
% kpe_ext(Fpe_ext<=0) = 0;

kpe_flex = kp / (e0*params_OS.Fp(2)) * exp(kp * (lMtilda_flex - kFpe*10) / e0);
% kpe_flex(Fpe_flex<=0) = 0;

kpe = [kpe_ext; kpe_flex];

% FMltilda
% Active muscle force-length characteristics
b11 = params_OS.Fa(1); b21 = params_OS.Fa(2); b31 = params_OS.Fa(3); b41 = params_OS.Fa(4);
b12 = params_OS.Fa(5); b22 = params_OS.Fa(6); b32 = params_OS.Fa(7); b42 = params_OS.Fa(8);
b13 = 0.1;             b23 = 1;               b33 = 0.5*sqrt(0.5);   b43 = 0;

% ext 
num3_ext     = lMtilda_ext-b23;
den3_ext     = b33+b43*lMtilda_ext;
FMtilda3_ext = b13*exp(-0.5*num3_ext.^2./den3_ext.^2);

num1_ext     = lMtilda_ext-b21;
den1_ext     = b31+b41*lMtilda_ext;
FMtilda1_ext = b11*exp(-0.5*num1_ext.^2./den1_ext.^2);

num2_ext     = lMtilda_ext-b22;
den2_ext     = b32+b42*lMtilda_ext;
FMtilda2_ext = b12*exp(-0.5*num2_ext.^2./den2_ext.^2);

FMltilda_ext = FMtilda1_ext+FMtilda2_ext+FMtilda3_ext;

% flex 
num3_flex     = lMtilda_flex-b23;
den3_flex     = b33+b43*lMtilda_flex;
FMtilda3_flex = b13*exp(-0.5*num3_flex.^2./den3_flex.^2);

num1_flex     = lMtilda_flex-b21;
den1_flex     = b31+b41*lMtilda_flex;
FMtilda1_flex = b11*exp(-0.5*num1_flex.^2./den1_flex.^2);

num2_flex     = lMtilda_flex-b22;
den2_flex     = b32+b42*lMtilda_flex;
FMtilda2_flex = b12*exp(-0.5*num2_flex.^2./den2_flex.^2);

FMltilda_flex = FMtilda1_flex+FMtilda2_flex+FMtilda3_flex;

FMltilda      = [FMltilda_ext; FMltilda_flex]; 

% FMvtilda
vMtildamax = params_OS.MT(5,:); 

% active force-velocity
e1 = params_OS.Fv(1); 
e2 = params_OS.Fv(2);
e3 = params_OS.Fv(3);
e4 = params_OS.Fv(4); 

FMvtilda_ext  = e1*log((e2*vMtilda_ext./vMtildamax(1)+e3)+sqrt((e2*vMtilda_ext./vMtildamax(1)+e3).^2+1))+e4; % extensor
FMvtilda_flex = e1*log((e2*vMtilda_flex./vMtildamax(2)+e3)+sqrt((e2*vMtilda_flex./vMtildamax(2)+e3).^2+1))+e4; % flexor, change vMtilda to flexor 
FMvtilda      = [FMvtilda_ext; FMvtilda_flex];  
end

