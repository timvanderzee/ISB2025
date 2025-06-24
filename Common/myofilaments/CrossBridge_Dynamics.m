function[Q0dot, Q1dot, Q2dot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, Non, DRX)

% attaching
beta = f .* [1 0 w^2];

% detaching
c1 = [Q0; p; 2*q];
phi0 = -IGef{1}(c1,k1) -IGef{1}(c1,k2);  
phi1 = -IGef{2}(c1,k1) -IGef{2}(c1,k2);  
phi2 = -IGef{3}(c1,k1) -IGef{3}(c1,k2);  
phi = [phi0; phi1; phi2];

% cross-bridge dynamics
Q0dot = DRX .* (beta(1,1) .* (Non-Q0)) + phi(1,:);
Q1dot = DRX .* (beta(1,2) .* (Non-Q0)) + phi(2,:);
Q2dot = DRX .* (beta(1,3) .* (Non-Q0)) + phi(3,:);

end