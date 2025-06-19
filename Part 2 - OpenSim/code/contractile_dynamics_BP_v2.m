function[Qd, vM, Ldd, Nond, DRXd] = contractile_dynamics_BP_v2(a, FMltilda, x, vMT, kT, kP, cos_alpha, params)

% v2: with dampening
lMo = params(2);
lTs = params(3);

%% spatial integrals
gamma = 100; % length scaling

delta =   2.018;
parms.w = 0.2;
parms.k11 = 106.5;
parms.k12 = 2;
parms.k21 = 480.6;
parms.k22 = .23;
parms.f = 1000;
parms.vmtc = 0;
parms.d = 1e-05;
parms.Ca = 1;

parms.J1 = 0.96;
parms.J2 = 4 * parms.J1;
parms.JF = 70.8;

parms.koop = 19.4;
parms.kon = 34.7;
parms.koff = 80;

parms.gaussian.IGef{1} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)));
parms.gaussian.IGef{2} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*(c(2)-c(3)*k(2)/2);
parms.gaussian.IGef{3} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*((c(2)-c(3)*k(2)/2).^2+c(3)/2);

%% spatial integrals
Ld = x(4);
Fd =  FMltilda * a * parms.d * Ld;

eps = 1e-6;
Q(1,1) = x(2);
Q(1) = max(Q(1),eps);
Q(2,1) = x(1) / delta - Q(1) - Fd;
Q(3,1) = max(x(3), eps);

% get gaussian coefficients
% mean of the distribution
p = Q(2)/Q(1); % Eq. 52

% the standard deviation of the distribution
if length(x)>2
    q = max(Q(3)/Q(1) - p^2, eps);  % Eq. 52
else
    q = .01;
end

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

%% thin filament
Ntot = max(a * FMltilda, eps);
Non = max(x(5), eps);
Noff = max(Ntot - Non, 0);
% Non = Ntot;

% Campbell 2018
Jon = parms.Ca *    parms.kon  * (Noff)  * (1 + max(parms.koop * Non/Ntot, 0));
Joff =              parms.koff * (Non-Q(1))    * (1 + max(parms.koop * (Noff)/Ntot,0));
Nond = max(Jon,0) - max(Joff,0);

%% thick filament
DRX = x(6);
F = Q(1) + Q(2);
J1 = (parms.J1 + parms.JF * max(F,0)/max(Ntot,1e-3)) .* (1-DRX);
J2 = parms.J2 .* DRX;
DRXd = J1 - J2; % this is a mistake!

%% cross-bridge dynamics
Q0dot = Non * beta(1,1) + phi1(1) + phi2(1);
Q1dot = Non * beta(1,2) + phi1(2) + phi2(2);

% velocity - independent derivative
Fdot = Q1dot + Q0dot;

% scale stuff
vMTtilda = vMT / lMo; % m/s to l0/s
vMT_CB = vMTtilda * gamma; % l0/s to h/s

kTc = kT .* lMo(:)./lTs(:); % F0/lTs to F0/lMo
% kTc = kT .* lTs(:)./lMo(:); % F0/lTs to F0/lMo
kT_CB = kTc / (delta*gamma); % F0/lTs to F0/lMo
kP_CB = kP / (delta*gamma); % F0/lTs to F0/lMo

% velocity from force constraint
% Ld  = (vMT_CB * kT_CB - Fdot.*cos_alpha) / (Q(1) + kT_CB./cos_alpha + kP_CB);
Ldd = (-Ld * (Q(1) + kT_CB./cos_alpha + kP_CB) + vMT_CB * kT_CB - Fdot .* cos_alpha) / (Non * parms.d);
 
% change in distribution
Qr = [0; Q];
Qd = nan(size(Q));

for i = 1:3
    Qd(i,1)  = Non * beta(i) + phi1(i) + phi2(i)...
             + (i-1) * Ld * Qr(i);
end

vM = Ld / gamma .* lMo(:);


end