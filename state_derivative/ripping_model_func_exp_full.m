function[Xd, F,n,xi] = ripping_model_func_exp_full(t, x, parms, Ca)

    if nargin < 4
        Ca = parms.Ca; % expressed in uM
    end

    % activation from Ca
    Act = parms.actfunc(Ca, parms);

    % retrieve states
    Non = x(1);
    n = x(2:(parms.nbins+1)); 
    R = x(parms.nbins+2);
    DRX = x(parms.nbins+3);
    lce = x(parms.nbins+4);

    % safety
    n(n<0) = 0;
    
    % compute moments
    % displacement from start
    xi = parms.xi0 + (lce - parms.lce0);
%     iRel = ((xi(:) < 2) & (xi(:) > -1)) | (abs(n(:)) > 1e-16);
    
    iRel = 1:length(parms.xi0);
%     iRel = abs(n(:)) > 1e-16;
    
    % compute moments
    Q = nan(1,3);
    for i = 1:3
        Q(i) = trapz(xi, xi.^(i-1) .* n');
    end
    
    % only select relevant portion
    parms.xi = xi(iRel);
    ns = n(iRel);
    
    %% thin filament activation
    if parms.max
        Non = 1;
    end
    
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
    % attachment and detachment at each strain
    beta = parms.f_func(parms.xi, parms.f, parms.w);
    phi1 = -parms.f_func(parms.xi, parms.f, parms.w) .* Q(1);
    phi2_rd = -(parms.g_func(parms.xi, parms.k11, -parms.k12) + parms.g_func(parms.xi, parms.k21, parms.k22)) .* ns';

    if parms.forcible_detachment
        phi2_fd = -(parms.k_func(parms.xi, parms.k, parms.dLcrit)) .* ns' + parms.f_func(parms.xi, parms.b, parms.w) * R;
    else
        phi2_fd = zeros(size(parms.xi));
    end
    
    % add both types of detachment
    phi2 = phi2_rd + phi2_fd;
    
    % change in cross-bridge attachment
    ndot = DRX * (Non * beta + phi1) + phi2;

    % state derivates
    nd = zeros(size(parms.xi0'));
    nd(iRel,1) = ndot';
    
    Rd = trapz(parms.xi, -phi2_fd);

    %% cross-bridge dynamics
    % first determine contraction velocity
    F = max(Q(2) + parms.ps * Q(1), 0);
    Q0dot = trapz(parms.xi, ndot);
    Q1dot = trapz(parms.xi, ndot .* parms.xi);    

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
    
    %% super-relaxed state dynamics
    SRX = min(1,(1 - DRX - Q(1) - R)); % super-relaxed = total - on - bound
    J1 = parms.J1 * (1 + parms.JF * F) * SRX; % Eq (3)
    J2 = parms.J2 * DRX; % Eq (4)

    % gains from super-relaxed, loses to super-relaxed, gains and loses to binding
    Dd = J1 - J2 - Q0dot;
        
    %% combined state derivative vector
    Xd = [Nond; nd; Rd; Dd; Ld; parms.vmtc];

    % remove length if needed
    Xd = Xd(1:length(x));

end