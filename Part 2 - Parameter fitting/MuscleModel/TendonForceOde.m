function [dfse, FT] = TendonForceOde(a,fse,lMT,vMT,params, Fvparam, Fpparam, Faparam)
% Try to vectorize it
FMo = ones(size(a,1),1)*params(1,:);
lMo = ones(size(a,1),1)*params(2,:);
lTs = ones(size(a,1),1)*params(3,:);
alphao = ones(size(a,1),1)*params(4,:);
vMmax = ones(size(a,1),1)*params(5,:);

FT = fse .* FMo;
FT(FT<0) = 0
lTtilda = log(5*(fse + 0.25))/35 + 0.995;

lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilda).^2)

cos_alpha = (lMT-lTs.*lTtilda)./lM;

lMtilda = lM./lMo;

% Calculate active force length
% FMltilda = muscleactiveforcelength(lMtilda);
% kShapeActive = 0.5;   
% 
% x=(lMtilda-1).*(lMtilda-1);
% FMltilda = exp(-x/kShapeActive);

% Gaussians
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilda-b23;
den3 = b33+b43*lMtilda;
FMtilda3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilda-b21;
den1 = b31+b41*lMtilda;
FMtilda1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilda-b22;
den2 = b32+b42*lMtilda;
FMtilda2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilda = FMtilda1+FMtilda2+FMtilda3

% Calculate passive force length
% Fpe = musclepassiveforcelength(lMtilda);
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilda - 0.10e1) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2)

FMce = fse./cos_alpha-Fpe;
FMce(FMce<0) = 0;
FMvtilda = FMce./(a.*FMltilda)

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

vMtilda = 1/e2*(sinh((FMvtilda-e4)/e1)-e3);

vM = vMtilda.*vMmax
vT = vMT-vM./cos_alpha

dfse = 7.*exp(35*(lTtilda - 0.995)).*vT./lTs;

pause;

return