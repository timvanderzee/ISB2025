function [error] = MuscleEquilibrium_alt_v2(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, f, k11, k12, k21, k22, a, vMtilda, DRX)

% points where integrals is evaluated
k1 = [k11 k12];
k2 = [k21 -k22];
w = 0.2;

gamma = 100; % length scaling

IGef{1} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3)));
IGef{2} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*(c(2,:)-c(3,:)*k(2)/2);
IGef{3} = @(c,k)(c(1,:)*k(1).*exp(min(c(3,:)*k(2)^2/4-c(2,:)*k(2),3))).*((c(2,:)-c(3,:)*k(2)/2).^2+c(3,:)/2);

% Compute first-order moment
Q1 = p .* Q0;



% Compute Qdot
[Q0dot, Q1dot, Q2dot] = CrossBridge_Dynamics(Q0, p, q, f, w, k1, k2, IGef, a, DRX);

error_Q0 = dQ0dt - Q0dot;
error_Q1 = dQ1dt - (Q1dot + 1 * vMtilda * gamma .* Q0);
error_Q2 = dQ2dt - (Q2dot + 2 * vMtilda * gamma .* Q1);

error = [error_Q0; error_Q1; error_Q2];

end