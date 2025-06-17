function[Qd, vM] = contractile_dynamics_BP(a, FMltilda, x, vMT, kT, kP, cos_alpha, params)

% without dampening
lMo = params(2);
lTs = params(3);

%% spatial integrals
delta = 2.2897; % force scaling
gamma = 100; % length scaling

eps = 1e-6;
Q(1,1) = x(2);
Q(1) = max(Q(1),eps);
Q(2,1) = x(1) / delta - Q(1);
Q(3,1) = x(3);

parms.w = 0.2;
parms.k11 = 559.9780;
parms.k12 = 2;
parms.k21 = 120.8824;
parms.k22 = 2;
parms.f = 500;
parms.vmtc = 0;

parms.gaussian.IGef{1} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)));
parms.gaussian.IGef{2} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*(c(2)-c(3)*k(2)/2);
parms.gaussian.IGef{3} =  @(c,k)(c(1)*k(1)*exp(min(c(3)*k(2)^2/4-c(2)*k(2),2)))*((c(2)-c(3)*k(2)/2).^2+c(3)/2);

%% spatial integrals
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

%% cross-bridge dynamics
Q0dot = FMltilda * a * beta(1,1) + phi1(1) + phi2(1);
Q1dot = FMltilda * a * beta(1,2) + phi1(2) + phi2(2);

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
Ld  = (vMT_CB * kT_CB - Fdot.*cos_alpha) / (Q(1) + kT_CB./cos_alpha + kP_CB);

% change in distribution
Qr = [0; Q];
Qd = nan(size(Q));

for i = 1:3
    Qd(i,1)  = FMltilda * a * beta(i) + phi1(i) + phi2(i)...
             + (i-1) * Ld * Qr(i);
end

vM = Ld / gamma .* lMo(:);

end