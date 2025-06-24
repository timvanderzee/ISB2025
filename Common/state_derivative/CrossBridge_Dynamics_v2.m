function[Q0dot, Q1dot, Q2dot] = CrossBridge_Dynamics_v2(Q0, Q1, Q2, f, w, k1, k2, IGef, Non, DRX)

Q00 = max(Q0, 1e-6);

p = Q1./Q00; % Eq. 52

% the standard deviation of the distribution
q = Q2./Q00 - p.^2;  % Eq. 52

% attaching
beta = f .* [1 0 w^2];
c1 = [Q0; p; 2*q];

% detaching
phi0 = -IGef{1}(c1,k1) -IGef{1}(c1,k2);  
phi1 = -IGef{2}(c1,k1) -IGef{2}(c1,k2);  
phi2 = -IGef{3}(c1,k1) -IGef{3}(c1,k2);  
phi = [phi0; phi1; phi2];

%% cross-bridge dynamics
% first determine contraction velocity
Q0dot = DRX .* (beta(1,1) .* (Non-Q0)) + phi(1,:);
Q1dot = DRX .* (beta(1,2) .* (Non-Q0)) + phi(2,:);
Q2dot = DRX .* (beta(1,3) .* (Non-Q0)) + phi(3,:);

end