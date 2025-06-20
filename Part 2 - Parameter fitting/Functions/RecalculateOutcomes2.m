function [output] = RecalculateOutcomes2(R, info, N_1)
% Function to recalculate Forces and Lengths
% SRS is added (Compared to version 1)

%% Input
% States 
x            = R.x;           
xd           = R.xd;  
xdd          = R.xdd; 
lMtilda_ext  = R.lMtilda(1,:); 
lMtilda_flex = R.lMtilda(2,:); 

% Controls
vMtilda_ext      = R.vMtilda(1,:);     
vMtilda_flex     = R.vMtilda(2,:); 
lM_projected_ext = R.lMprojected(1,:); 
lM_projected_flex= R.lMprojected(2,:);

% States for Reflexes
Fsrs_del         = R.Fsrs_del;
dFsrs_deldt      = R.dFsrs_deldt; 
Fce_del          = R.Fce_del;
dFce_deldt       = R.dFce_deldt;

% Params
a_ext    = R.a_ext; 
a_flex   = R.a_flex;
B        = R.B;
kFpe     = R.kFpe; 
kR       = R.kR; 

% OS params
FMo      = R.OS.MT(1,:);
lMo      = R.OS.MT(2,:); 
lTs      = R.OS.MT(3,:); 
vMtildamax = R.OS.MT(5,:);

% Inertia 
m        = R.OS.inert.mass_OS; 
l        = R.OS.inert.lc_OS;
I        = R.OS.inert.I_OS; 

% Additional
offset   = R.exp.offset; 
m_offset = mean(offset);
shift    = R.shift; 
coeff_LMT_ma = R.coeff; 

% Cost function
output.J_kern = R.J -  0.001 * (sumsqr(R.vMtilda)) - 0.1*(a_ext + a_flex + kR); 

%% Results 
% Muscle Tendon Lengths
lMT_ext     = coeff_LMT_ma(1,1) + coeff_LMT_ma(2,1)*(x+m_offset) + coeff_LMT_ma(3,1)*(x+m_offset).^2 + coeff_LMT_ma(4,1)*(x+m_offset).^3;
output.lMT_ext = lMT_ext; 
lMT_flex     = coeff_LMT_ma(1,2) + coeff_LMT_ma(2,2)*(x+m_offset) + coeff_LMT_ma(3,2)*(x+m_offset).^2 + coeff_LMT_ma(4,2)*(x+m_offset).^3;
output.lMT_flex = lMT_flex; 

% Moment Arms
MA_ext  = -coeff_LMT_ma(2,1) + -coeff_LMT_ma(3,1)*(x+m_offset) + -coeff_LMT_ma(4,1)*(x+m_offset).^2;
output.MA_ext = MA_ext; 
MA_flex  = -coeff_LMT_ma(2,2) + -coeff_LMT_ma(3,2)*(x+m_offset) + -coeff_LMT_ma(4,2)*(x+m_offset).^2;
output.MA_flex = MA_flex; 

% lT 
lT_ext   = lMT_ext - lM_projected_ext; 
output.lT_ext = lT_ext; 
lT_flex   = lMT_flex - lM_projected_flex; 
output.lT_flex = lT_flex; 

% lTtilda
lTtilda_ext  = lT_ext./lTs(1); 
output.lTtilda_ext = lTtilda_ext; 
lTtilda_flex  = lT_flex./lTs(2); 
output.lTtilda_flex = lTtilda_flex; 

% lM 
lM_ext  = lMtilda_ext.* lMo(1);
output.lM_ext = lM_ext; 
lM_flex  = lMtilda_flex.* lMo(2);
output.lM_flex = lM_flex; 

% Fse
fse_ext      = (exp(35*(lTtilda_ext - 0.995)))/5-0.25 + shift;
output.fse_ext = fse_ext; 
fse_flex     = (exp(35*(lTtilda_flex - 0.995)))/5-0.25 + shift;
output.fse_flex = fse_flex; 

% FT
FT_ext       = FMo(1).* fse_ext; 
output.FT_ext = FT_ext; 
FT_flex      = FMo(2).* fse_flex; 
output.FT_flex = FT_flex; 

% Fpe
e0      = 0.6;   kpe = 4;
t5_ext       = exp(kpe * (lMtilda_ext - kFpe*10) / e0);         
Fpe_ext      = ((t5_ext - 0.10e1) - R.OS.Fp(1)) / R.OS.Fp(2);    
t5_flex      = exp(kpe * (lMtilda_flex - kFpe*10) / e0);   
Fpe_flex     = ((t5_flex - 0.10e1) - R.OS.Fp(1)) / R.OS.Fp(2);  
output.Fpe_ext  = Fpe_ext; 
output.Fpe_flex = Fpe_flex; 

% FMltilda
% Active muscle force-length characteristics
b11 = R.OS.Fa(1); b21 = R.OS.Fa(2); b31 = R.OS.Fa(3);       b41 = R.OS.Fa(4);
b12 = R.OS.Fa(5); b22 = R.OS.Fa(6); b32 = R.OS.Fa(7);       b42 = R.OS.Fa(8);
b13 = 0.1;        b23 = 1;          b33 = 0.5*sqrt(0.5);    b43 = 0;

% Ext 
num3_ext = lMtilda_ext-b23;
den3_ext = b33+b43*lMtilda_ext;
FMtilda3_ext = b13*exp(-0.5*num3_ext.^2./den3_ext.^2);

num1_ext = lMtilda_ext-b21;
den1_ext = b31+b41*lMtilda_ext;
FMtilda1_ext = b11*exp(-0.5*num1_ext.^2./den1_ext.^2);

num2_ext = lMtilda_ext-b22;
den2_ext = b32+b42*lMtilda_ext;
FMtilda2_ext = b12*exp(-0.5*num2_ext.^2./den2_ext.^2);

FMltilda_ext = FMtilda1_ext+FMtilda2_ext+FMtilda3_ext;
output.FMltilda_ext = FMltilda_ext;

% Flex
num3_flex = lMtilda_flex-b23;
den3_flex = b33+b43*lMtilda_flex;
FMtilda3_flex = b13*exp(-0.5*num3_flex.^2./den3_flex.^2);

num1_flex = lMtilda_flex-b21;
den1_flex = b31+b41*lMtilda_flex;
FMtilda1_flex = b11*exp(-0.5*num1_flex.^2./den1_flex.^2);

num2_flex = lMtilda_flex-b22;
den2_flex = b32+b42*lMtilda_flex;
FMtilda2_flex = b12*exp(-0.5*num2_flex.^2./den2_flex.^2);

FMltilda_flex = FMtilda1_flex+FMtilda2_flex+FMtilda3_flex;
output.FMltilda_flex = FMltilda_flex;

% FMvtilda
% active force-velocity
e1 = R.OS.Fv(1); 
e2 = R.OS.Fv(2);
e3 = R.OS.Fv(3);
e4 = R.OS.Fv(4); 

FMvtilda_ext  = e1*log((e2*vMtilda_ext./vMtildamax(1)+e3)+sqrt((e2*vMtilda_ext./vMtildamax(1)+e3).^2+1))+e4; % extensor
FMvtilda_flex = e1*log((e2*vMtilda_flex./vMtildamax(2)+e3)+sqrt((e2*vMtilda_flex./vMtildamax(2)+e3).^2+1))+e4; % flexor
output.FMvtilda_ext  = FMvtilda_ext; 
output.FMvtilda_flex = FMvtilda_flex; 

% Fsrs  
kSRS    = info.kSRS; 
stretch = lMtilda_ext(1:N_1)- lMtilda_ext(1); 

tan_val_incr = (0.5* tanh(1000*(5.7e-3-stretch))+0.5).*FMltilda_ext(1:N_1)*a_ext*kSRS.*stretch; 
tan_val_plat = (0.5* tanh(1000*(stretch-5.7e-3))+0.5).*FMltilda_ext(1:N_1)*a_ext*kSRS*5.7e-3; 
Fsrs_f1 = tan_val_incr + tan_val_plat; 

output.Fsrs_f1 = Fsrs_f1; 
output.Fsrs    = [Fsrs_f1 R.Fsrs]; 

% Reflexes 
tres   = a_ext*FMltilda_ext(1)*FMvtilda_ext(1) + output.Fsrs(1); 

a_refl = kR * (Fce_del-tres); 
a      = a_ext + (0.5*tanh(1000*a_refl)+0.5).*a_refl;

output.a_refl = a_refl;
output.a      = a; 

% Fce
Fce_ext = a_ext.* FMltilda_ext.* FMvtilda_ext + output.Fsrs;
%(a_ext + kR * Fsrs_del + kY * dFsrs_deldt - R.tres).* FMltilda_ext.* FMvtilda_ext + output.Fsrs; 
output.Fce_ext = Fce_ext; 

Fce_flex = a_flex.* FMltilda_flex.* FMvtilda_flex; 
output.Fce_flex = Fce_flex; 

% FM
FM_ext = Fce_ext + Fpe_ext; 
output.FM_ext = FM_ext; 

FM_flex = Fce_flex + Fpe_flex; 
output.FM_flex = FM_flex; 

% Muscle torque
T_muscle_ext = -FT_ext.*MA_ext;
output.T_muscle_ext = T_muscle_ext;

T_muscle_flex = -FT_flex.*MA_flex;
output.T_muscle_flex = T_muscle_flex;

% Gravity torque
g    = 9.81;
T_Fz = m*g*l*cos(x);
output.T_Fz      = T_Fz; 

% Inertie torque
T_inert = xdd*I ; 
output.T_inert = T_inert;

% Damping torque
T_damping = xd*B; 
output.T_damping = T_damping; 
end

