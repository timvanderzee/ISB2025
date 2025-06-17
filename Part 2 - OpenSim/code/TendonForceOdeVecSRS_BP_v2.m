function [dx, FT] = TendonForceOdeVecSRS_BP_v2(t,x,t_input,A,LMT,VMT,params, m, Fvparam, Fpparam, Faparam, lMtilda_isom, ksrs, kT)

% length as a state
% Try to vectorize it
FMo = params(1,m);
lMo = params(2,m);
lTs = params(3,m);
alphao = params(4,m);
% vMmax = params(5,m); % in m/s

a = interp1(t_input, A, t);
lMT = interp1(t_input, LMT, t); % [m]
vMT = interp1(t_input, VMT, t); % [m/s]

lMtilda = x(1);
lM = lMtilda * lMo;
lMp = sqrt(lM.^2 - (lMo.*sin(alphao)).^2); % [m]

% calc fse
lT = lMT - lMp; % [m]
lTtilda = lT ./ lTs;
fse = kT * (lTtilda-1);

FT = fse .* FMo;
% lTtilda = log(5*(fse + 0.25))/35 + 0.995;
% lTtilda = fse/kT + 1;

cos_alpha = (lMT-lTs.*lTtilda)./lM;
% cos_alpha = 1;
Fpe = 0;
Fsrs = 0;
FMltilda = 1;

% lMtilda = lM./lMo; % []

% Calculate active force length
% FMltilda = muscleactiveforcelength(lMtilda);
% kShapeActive = 0.5;   
% 
% x=(lMtilda-1).*(lMtilda-1);
% FMltilda = exp(-x/kShapeActive);

% Gaussians
% b11 = Faparam(1);
% b21 = Faparam(2);
% b31 = Faparam(3);
% b41 = Faparam(4);
% b12 = Faparam(5);
% b22 = Faparam(6);
% b32 = Faparam(7);
% b42 = Faparam(8);
% 
% b13 = 0.1;
% b23 = 1;
% b33 = 0.5*sqrt(0.5);
% b43 = 0;
% num3 = lMtilda-b23;
% den3 = b33+b43*lMtilda;
% FMtilda3 = b13*exp(-0.5*num3.^2./den3.^2);
% 
% num1 = lMtilda-b21;
% den1 = b31+b41*lMtilda;
% FMtilda1 = b11*exp(-0.5*num1.^2./den1.^2);
% 
% num2 = lMtilda-b22;
% den2 = b32+b42*lMtilda;
% FMtilda2 = b12*exp(-0.5*num2.^2./den2.^2);
% 
% FMltilda = FMtilda1+FMtilda2+FMtilda3;
% FMltilda = 1;
% 
% % Force from short-range stiffness
% dlM = (0.5*tanh(1000*(lMtilda - lMtilda_isom))+0.5) .*(lMtilda - lMtilda_isom);
% dlM = (0.5*tanh(1000*(-dlM+5.7*10^(-3)))+ 0.5).*dlM + (0.5*tanh(1000*(dlM-5.7*10^(-3)))+0.5)*5.7*10^(-3);
% % dlM = min(max(0,lMtilda - lMtilda_isom),5.7*10^(-3));
% Fsrs = ksrs * a .* FMltilda.*dlM;
% Fsrs = 0;
% 
% % Calculate passive force length
% % Fpe = musclepassiveforcelength(lMtilda);
% e0 = 0.6;
% kpe = 4;
% t5 = exp(kpe * (lMtilda - 0.10e1) / e0);
% Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);
% kP = kpe / (e0*Fpparam(2)) * exp(kpe * (lMtilda - 1) / e0);

FMce = fse./cos_alpha-Fpe-Fsrs;

%% spatial integrals
delta = 2.4; % force scaling
gamma = 100; % length scaling

Q(1,1) = x(2);
Q(3,1) = x(3);
Q(2,1) = FMce / delta - Q(1);

eps = 1e-6;

Q(1) = max(Q(1),eps);
% Q(2) = max(Q(2), -Q(1));
% Q(3) = max(Q(3),0);

parms.w = 0.1667;
parms.k11 = 671;
parms.k12 = 2;
parms.k21 = 127;
parms.k22 = 2;
parms.f = 500;
parms.vmtc = 0;

parms.gaussian.IGef{1} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),1e3)));
parms.gaussian.IGef{2} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),1e3)))*(c(2)-c(3)*k(2)/2);
parms.gaussian.IGef{3} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),1e3)))*((c(2)-c(3)*k(2)/2).^2+c(3)/2);

% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

%% spatial integrals
% get gaussian coefficients
% mean of the distribution
p = Q(2)/Q(1); % Eq. 52

% the standard deviation of the distribution
q = max(Q(3)/Q(1) - p^2, 0);  % Eq. 52
% q = .1;
c1 = [Q(1) p 2*q];

% moments of the attachment function
C = [1 0 parms.w^2];

% points were integrals is evaluated
k1 = [parms.k11 parms.k12];
k2 = [parms.k21 -parms.k22];

% analytical
phi1 = nan(1,3);
phi2 = nan(1,3);
beta = nan(1,3);

for i = 1:3
    % breaking
    phi2(1,i) = -parms.gaussian.IGef{i}(c1,k1) -parms.gaussian.IGef{i}(c1,k2); 

    % attaching
    beta(1,i) = parms.f * C(i);
    phi1(1,i) = -parms.f * C(i) * Q(1);

end

%% cross-bridge dynamics
Q0dot = a * beta(1,1) + phi1(1) + phi2(1);
Q1dot = a * beta(1,2) + phi1(2) + phi2(2);

% velocity - independent derivative
Fdot  = Q1dot + Q0dot;

% scale stuff
% vMTtilda = vMT / lMo; % m/s to l0/s
vMT_CB = vMT * gamma; % l0/s to h/s

% kTc = kT .* lMo(:)./lTs(:); % F0/lTs to F0/lMo
kTc = kT .* lTs(:)./lMo(:); % F0/lTs to F0/lMo
kT_CB = kTc / (delta*gamma); % F0/lTs to F0/lMo
% kP_CB = kP / gamma; % F0/lTs to F0/lMo

% kTc = kT;

% velocity from force constraint
Ld  = (vMT_CB * kT_CB - Fdot.*cos_alpha) / (Q(1) + kT_CB./cos_alpha);
% Ld  = -Fdot / (Q(1) + kT_CB);

% change in distribution
Qr = [0; Q];
Qd = nan(size(Q));

for i = 1:3
    Qd(i,1)  = a * beta(i) + phi1(i) + phi2(i)...
             + (i-1) * Ld * Qr(i);
end


% vM = Ld / gamma * lMo; % h/s to m/s

%% combined state derivative vector
% vT = vMT - vM; % m/s
% vT = vMT - Ld .* lMo(:) / gamma;

% dfse = kT.*vT./lTs;
% dfse = kT.*vT;

dx(1,1) = Ld / gamma;
dx(2,1) = Qd(1);
dx(3,1) = Qd(3);

return