function[Xd, Ftot, n] = twostate_model_func_exp_v2(t, x, parms, Act)

    % retrieve states
    Q = x(1:3); 
    Ld = x(4);
    lmtc = x(5);
    
    eps = 1e-6;

    Q(1) = max(Q(1),eps);
    Q(2) = max(Q(2), -Q(1));
    Q(3) = max(Q(3),eps.^2);
    
    %% spatial integrals
    % get gaussian coefficients
    % mean of the distribution
    p = Q(2)/max(Q(1), eps); % Eq. 52

    % the standard deviation of the distribution
    q = sqrt(max(Q(3)/max(Q(1), eps) - (Q(2)/max(Q(1), eps))^2, eps));  % Eq. 52

    % G = @(x,c) c(1)*exp(-(x-c(2)).^2/c(3))
    c1 = [Q(1) p 2*q^2];

    % moments of the attachment function
    C = [1 0 parms.w^2];
    C2 = [1 0 parms.w^2];
    
    % points were integrals is evaluated
    k1 = [parms.k11 parms.k12];
    k2 = [parms.k21 -parms.k22];
%     k3 = [parms.k31 -parms.k32];

    % analytical
    phi1 = nan(1,3);
    phi2s = nan(2,3);
    beta = nan(1,3);
    
    for i = 1:3

        % breaking
        phi2s(1,i) = -parms.gaussian.IGef{i}(c1,k1) -parms.gaussian.IGef{i}(c1,k2) - parms.g1 * Q(i); 
       

        % ripping
        if parms.forcible_detachment
            phi2s(2,i) = -parms.k * (Q(i)/2 - parms.gaussian.IG{i}(parms.dLcrit, c1)) + parms.b * C2(i) * R(1);
        else
            phi2s(2,i) = 0;
        end
        
        % attaching
        beta(1,i) = parms.f * C(i);
        phi1(1,i) = -parms.f * C(i) * Q(1);
    
    end
    
    if sum(~isfinite(phi2s)) > 0
        keyboard
    end
    
    Rd = -phi2s(2,1);
    phi2 = sum(phi2s);

    %% cross-bridge dynamics
    % first determine contraction velocity
    F = max(Q(2) + parms.ps * Q(1), 0);
    
    Q0dot = Act * beta(1,1) + phi1(1) + phi2(1);
    Q1dot = Act * beta(1,2) + phi1(2) + phi2(2);

    if ~parms.no_tendon
        % length dependence
        dlse = parms.Lse_func(F, parms);

        % stiffness
        kse = parms.kse_func(dlse, parms);

        % velocity - independent derivative
        Fdot  = Q1dot + parms.ps * Q0dot;

        % velocity from force constraint
%         Ld  = (parms.vmtc * kse - Fdot) / (Q(1) + kse);
        
        Ldd = (-Ld * (Q(1) + kse) + parms.vmtc * kse - Fdot) / parms.d;
            
    else
        Ld = parms.c * parms.vmtc;
    end
    
    % change in distribution
    Qr = [0; Q];
    Qd = nan(size(Q));
    
    for i = 1:3
        Qd(i,1)  = Act * beta(i) + phi1(i) + phi2(i)...
                 + (i-1) * Ld * Qr(i);
    end
        
    % total force
    F_pas = parms.Fpe_func(lmtc, parms);
    F_act = F * parms.Fscale;
    Ftot = F_act(:) + F_pas(:);
    
    %% combined state derivative vector
    Xd = [Qd; Ldd; parms.vmtc];

    % remove length if needed
    Xd = Xd(1:length(x));
    
    if nargout > 2
        n = Q(1) ./ (sqrt(2*pi)*q) * exp(-((parms.xi-p).^2) / (2*q^2));  % Eq. 52, modified
    end

%     if sum(isnan(Xd)) > 0
%         keyboard
%     end
%     F = Q(2);

end