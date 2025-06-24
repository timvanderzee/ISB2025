function [Q0, Q1, xi] = force_from_distribution(Q, lce, parms)

if length(Q) < 4
    Q0 = Q(1);
    Q1 = Q(2);
    xi = parms.xi;
    
else
    
    n = Q;
    
    % displacement from start
    xi = parms.xi0 + (lce - parms.lce0);
%     iRel = ((xi(:) < 2) & (xi(:) > -1)) | (abs(n(:)) > 1e-8);
    iRel = 1:length(parms.xi0);
    
    % only select relevant portion
    parms.xi = xi(iRel);
    ns = n(iRel);
   
    % compute moments
    Q = trapz(parms.xi(:), [ns parms.xi(:).*ns]);
    Q0 = Q(1);
    Q1 = Q(2);
    xi = parms.xi0 + (lce - parms.lce0);

end




end