function[Xd, Ftot] = ripping_model_func_exp(t, x, parms, Ca)

    if nargin < 4
        Ca = parms.Ca; % expressed in uM
    end

    % activation from Ca
    Act = parms.actfunc(Ca, parms);

    % retrieve states
    Non = x(1);
    Q = x(2:4); 
    R = x(5);
    DRX = x(6);
    lmtc = x(7);
    lce = x(8);
    
    eps = 1e-6;

    Q(1) = max(Q(1),eps);
    Q(2) = max(Q(2), -Q(1));
    Q(3) = max(Q(3),eps.^2);
    
%     xi = parms.xi;
%     n = n_func(Q, parms.xi);

    %% thin filament activation
    if parms.max
        Non = 1;
    end

    % get filament overlap %
    parms.Noverlap = calc_filament_overlap(parms, lce);
    
    % quantities 
    Noff = parms.Noverlap - Non;
    Noff(Noff<0) = 0; % could overshoot due to fast dynamics
    
    % thin filament dynamics    
    Jon     = Act * parms.kon * Noff * (1 + parms.koop * (Non/parms.Noverlap)); % Eq (1)
    Joff    = parms.koff * (Non - Q(1)) *     (1 + parms.koop * Noff/parms.Noverlap); % Eq (2)
    Nond    = Jon - Joff; % Eq (7)
    
    if parms.max
        Nond = 0;
    end
    
    % in case we don't do thin filament cooperative activation
    if (parms.kon == 0) && (parms.koff == 0) && (parms.koop == 0)
        Non = Act;
    end

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
    
    Q0dot = DRX * (Non * beta(1,1) + phi1(1)) + phi2(1);
    Q1dot = DRX * (Non * beta(1,2) + phi1(2)) + phi2(2);

    if ~parms.no_tendon
        % length dependence
        dlse = parms.Lse_func(F, parms);

        % stiffness
        kse = parms.kse_func(dlse, parms);

        % velocity - independent derivative
        Fdot  = Q1dot + parms.ps * Q0dot;

        % velocity from force constraint
        Ld  = (parms.vmtc * kse - Fdot) / (Q(1) + kse);
            
    else
        Ld = parms.c * parms.vmtc;
    end
    
    % change in distribution
    Qr = [0; Q];
    Qd = nan(size(Q));
    
    for i = 1:3
        Qd(i,1)  = DRX * (Non * beta(i) + phi1(i)) + phi2(i)...
                 + (i-1) * Ld * Qr(i);
    end

    %% super-relaxed state dynamics
    SRX = min(1,(1 - DRX - Q(1) - R)); % super-relaxed = total - on - bound
    J1 = parms.J1 * (1 + parms.JF * F) * SRX; % Eq (3)
    J2 = parms.J2 * DRX; % Eq (4)

    % gains from super-relaxed, loses to super-relaxed, gains and loses to binding
    Dd = J1 - J2 - Q0dot;
        
    % total force
    F_pas = parms.Fpe_func(lmtc, parms);
    F_act = F * parms.Fscale;
    Ftot = F_act(:) + F_pas(:);
    
    %% combined state derivative vector
    Xd = [Nond; Qd; Rd; Dd; parms.vmtc; Ld];

    % remove length if needed
    Xd = Xd(1:length(x));

%     if sum(isnan(Xd)) > 0
%         keyboard
%     end
%     F = Q(2);

end