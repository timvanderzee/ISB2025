function [Q0dot, Q1dot, Q2dot] = MuscleEquilibrium_alt_exp(Q0, Q1, Q2, f, k11, k12, k21, k22, a, vMtilda)

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
Q00 = max(Q0(j,:), 1e-6);

p = Q1./Q00; % Eq. 52

% the standard deviation of the distribution
q = Q2(j,:)./Q00 - p.^2;  % Eq. 52

% moments of the attachment function
c1 = [Q0(j,:); p; 2*q];

% attaching
beta = f .* C;

% IGef{1} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3)));
% IGef{2} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*(c(2,:)-c(3,:)*k(2)/2);
% IGef{3} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

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
Q1dot = (a * beta(1,2) + phi(2,:)) + 1 * vMtilda(j) * gamma * Q0; 
Q2dot = (a * beta(1,3) + phi(3,:)) + 2 * vMtilda(j) * gamma * Q1;

end