function [Q0, Q1, xi, n] = force_from_distribution(x, lce, parms)

if length(x) == 3 % distribution-moment
    Q = x;
    Q0 = Q(1);
    Q1 = Q(2);
    xi = parms.xi;
    n = parms.n_func(xi, Q, 1e-6);
    
else % discretized
    
    n = x;
    
    % displacement from start
    xi = parms.xi + (lce - parms.lce0);
   
    % compute moments
    Q = trapz(xi(:), [n xi(:).*n]);
    Q0 = Q(1);
    Q1 = Q(2);
end
end