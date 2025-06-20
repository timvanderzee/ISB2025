function [error_forcerate, error_Q0, error_Q2] = MuscleEquilibrium(a, FXB, cos_alpha, vMTtilda, kse, kpe, vMtilda, Q0, Q2, dQ0dt, dQ2dt)

j = 1;
% points where integrals is evaluated
k1 = [6 2];
k2 = [20 -0.5];
f = 50;
w = 0.2;

C = [1 0 w^2];

delta = 2; % force scaling
gamma = 100; % length scaling

% calc Q1
Q1 = (FXB(j,:) / delta) - Q0(j,:);

%% spatial integrals
% get gaussian coefficients
% mean of the distribution
Q00 = Q0(j,:);

p = Q1./Q00; % Eq. 52

% the standard deviation of the distribution
q = Q2(j,:)./Q00 - p.^2;  % Eq. 52

% moments of the attachment function
c1 = [Q0(j,:); p; 2*q];

% attaching
beta = f .* C;

IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

phi0 = -IGef{1}(c1,k1) -IGef{1}(c1,k2) - beta(1) * Q0(j,:);  
phi1 = -IGef{2}(c1,k1) -IGef{2}(c1,k2) - beta(2) * Q0(j,:);  
phi2 = -IGef{3}(c1,k1) -IGef{3}(c1,k2) - beta(3) * Q0(j,:);  

phi = [phi0; phi1; phi2];

%% cross-bridge dynamics
% first determine contraction velocity
Q0dot = (a * beta(1,1) + phi(1,:));
Q1dot = (a * beta(1,2) + phi(2,:));
Q2dot = (a * beta(1,3) + phi(3,:));

% velocity - independent derivative
FMdot  = (Q1dot + Q0dot) * delta;
km = Q0(j,:) * delta / gamma;

% velocity from force constraint
error_forcerate  = vMtilda(j,:) - (vMTtilda(j,:) .* kse(j,:) - FMdot.*cos_alpha(j,:)) ./ (km + kse(j,:)./cos_alpha(j,:) + kpe(j,:));
error_Q0 = dQ0dt(j,:) - Q0dot;
error_Q2 = dQ2dt(j,:) - (Q2dot + 2 * vMtilda(j) * gamma * Q1);

end