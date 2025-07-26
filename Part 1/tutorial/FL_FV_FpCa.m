%% load parameters and set up simulations 
addpath(genpath([pwd,'/../..']))

warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

parms.forcible_detachment = 0;
parms.kse = 0;
parms.kpe = 0;
parms.no_tendon = 1;
odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

save_hill_properties = 0; 
% if set to be 1, it will prompt to save Hill properties at the end

%% Run F-l, F-v, F-pCa protocols and save the outcome as splines.
hill_properties = struct();
ref_pCa = 4.5; % << pCa during F-L, F-v simulation 

for cond_itr = 1:3 %  1: F-L, 2: F-v, 3: F-pCa
    figure;

    if(cond_itr==1) % F-L
        l0 = [-0.5:0.1:0.5] *half_s_len_norm;
    elseif(cond_itr==2) % F-v
        l0 = [-1:0.2:1] *half_s_len_norm / 20;
    else % F-pCa
        pCa_list = [4.5,5,5.5:0.1:6.5,7:0.5:9];
        l0 = zeros(size(pCa_list));
    end

    % subplot to monitor time-length 
    subplot(2,2,1, 'colorOrder', winter(length(l0))); hold on
    % subplot to monitor time-force 
    subplot(2,2,3, 'colorOrder', winter(length(l0))); hold on

    F0 = nan(size(l0));

    % run the simulation 
    for k=1:length(l0)

        % activation protocol 
        if(cond_itr==3)
            pCa = pCa_list(k);
        else
            pCa = ref_pCa;
        end
        Ca = 10^(-pCa+6);

        % initialize model 
        model = @fiber_dynamics;
        parms.xi0 = linspace(-15,15,nbins);
        parms.xi = parms.xi0;
        parms.nbins = nbins;
        parms.xss = zeros(1,parms.nbins + 4);
        parms.xss(end-2) = 0.0909;

        % set initial length according to the protocol 
        parms.xss(end-1) = l0(k);
        parms.xss(end) = l0(k);
        parms.lce0 = l0(k);

        % velocity protocol 
        if(cond_itr==1||cond_itr==3) % isometric for F-L & F-pCa
            us = [0,0];
            Ts = [2,1];
        else % isokinetic for F-v
            us = [0,-l0(k)*20,0];
            Ts = [3,0.05,3];
        end

        x0 = parms.xss;
        t = 0;
        x = x0;

        % solve ODE to get XB distributions over time 
        temp_idx = 1;
        for i = 1:length(us)
            parms.vmtc = us(i);
            T = Ts(i);
            [tnew,xnew] = ode15s(model, [0 T], x0, odeopt, parms, Ca);
            x0 = xnew(end,:);
            t = [t; tnew(2:end)+t(end)];
            x = [x; xnew(2:end,:)];
            if(i==2), temp_idx = height(x); end
        end

        % calculate force
        F = nan(1,height(x));
        for i = 1:height(x)
            [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
        end

        % plot time-length
        subplot(2,2,1)
        plot(t, x(:,end)/half_s_len_norm)
        hold on

        % plot time-force
        subplot(2,2,3)
        plot(t,F)
        hold on
        pause(0.1)
        F0(k) = F(temp_idx);
    end

    % plot outcome F-l, F-v, F-pCa curves obtained from simulations 
    subplot(2,2,[2,4])
    
    if(cond_itr==1) % F-L
        plot(l0/half_s_len_norm, F0);
        xlabel('\Delta length (l_{opt})')
        ylabel('F (F_0)')
        title('F-l')
        hill_properties.FL_spline = spline(l0/half_s_len_norm, F0);
    elseif(cond_itr==2) % F-v
        plot(l0/half_s_len_norm/0.05, F0);
        xlabel('velocity (l_{opt}/s)')
        ylabel('F (F_0)')
        title('F-v')
        hill_properties.FV_spline = spline(l0/half_s_len_norm/0.05, F0/F0(abs(l0)<0.0001));
    else % F-pCa
        plot(pCa_list, F0);
        xlabel('pCa')
        ylabel('F (F_0)')
        title('F-pCa')
        xlim([4.5 9])
        hill_properties.FPca_spline = spline(pCa_list, F0/F0(1));
    end

    % label subplots 
    subplot(2,2,1)
    ylabel('\Delta length (l_{opt})')
    subplot(2,2,3)
    ylabel('F (F_0)')
    xlabel('time (s)')
end

% if needed, save hill-properties (can be loaded into other codes)
if(save_hill_properties)
    uisave('hill_properties');
end

