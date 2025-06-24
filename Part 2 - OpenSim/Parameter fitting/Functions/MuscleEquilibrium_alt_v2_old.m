function [error] = MuscleEquilibrium_alt_v2(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, f, k11, k12, k21, k22, a, vMtilda, DRX)

% no tendon
j = 1;

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];
w = 0.2;

C = [1 0 w^2];

gamma = 100; % length scaling

%% spatial integrals
% get gaussian coefficients
% mean of the distribution
% Q00 = max(Q0(j,:), 1e-6);
Q1 = p .* Q0;

% p = Q1./Q00; % Eq. 52

% the standard deviation of the distribution
% q = Q2(j,:)./Q00 - p.^2;  % Eq. 52

% moments of the attachment function
c1 = [Q0(j,:); p; 2*q];

% attaching
beta = f .* C;
% 
IGef{1} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(c(3,:)*k(2)^2/4-c(2,:)*k(2))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

phi0 = -IGef{1}(c1,k1) -IGef{1}(c1,k2);  
phi1 = -IGef{2}(c1,k1) -IGef{2}(c1,k2);  
phi2 = -IGef{3}(c1,k1) -IGef{3}(c1,k2);  

phi = [phi0; phi1; phi2];

%% cross-bridge dynamics
% first determine contraction velocity
Q0dot = DRX .* (beta(1,1) .* (a-Q0(j,:))) + phi(1,:);
Q1dot = DRX .* (beta(1,2) .* (a-Q0(j,:))) + phi(2,:);
Q2dot = DRX .* (beta(1,3) .* (a-Q0(j,:))) + phi(3,:);

error_Q0 = dQ0dt(j,:) - Q0dot;
error_Q1 = dQ1dt(j,:) - (Q1dot + 1 * vMtilda * gamma .* Q0);
error_Q2 = dQ2dt(j,:) - (Q2dot + 2 * vMtilda * gamma .* Q1);

% error_p = p(j,:) - (Q2dot + 2 * vMtilda(j) * gamma * Q1);

error = [error_Q0; error_Q1; error_Q2];

end