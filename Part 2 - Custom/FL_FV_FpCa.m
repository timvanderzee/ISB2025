% clear; clc; close all

addpath(genpath([pwd,'/..']))

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

for cond_itr = 3
    figure;

    if(cond_itr==1) % F-L
        l0 = [-0.5:0.1:0.5] *half_s_len_norm;
    elseif(cond_itr==2) % F-v
        l0 = [-1:0.2:1] *half_s_len_norm / 10;
    else % F-pCa
        l0 = zeros(1,10);
        pCa_list = linspace(4.5, 9, length(l0));
    end

    subplot(2,2,1, 'colorOrder', winter(length(l0))); hold on
    subplot(2,2,3, 'colorOrder', winter(length(l0))); hold on

    F0 = nan(size(l0));

    for k=1:length(l0)
        if(cond_itr==3)
            pCa = pCa_list(k);
        else
            pCa = 4.5;
        end
        Ca = 10^(-pCa+6);
        model = @fiber_dynamics;

        parms.xi0 = linspace(-15,15,nbins);
        parms.nbins = nbins;
        parms.xss = zeros(1,parms.nbins + 4);
        %         parms.xss = zeros(1,7);

        parms.xss(end-2) = 0.0909;

        parms.xss(end-1) = l0(k);
        parms.xss(end) = l0(k);

        parms.lce0 = l0(k);

         if(cond_itr==1||cond_itr==3) % F-L & F-pCa
             us = [0,0];
             Ts = [2,1];
         else % F-v
             us = [0,-l0(k)*10,0];
             Ts = [3,0.1,3];
         end

        x0 = parms.xss;
        t = 0;
        x = x0;

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

        F = nan(1,height(x));
        for i = 1:height(x)
            [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
        end

        subplot(2,2,1)
        plot(t, x(:,end)/half_s_len_norm)
        hold on

        subplot(2,2,3)
        plot(t,F)
        hold on
        pause(0.1)
        F0(k) = F(temp_idx);
    end

    subplot(2,2,[2,4])
    if(cond_itr==1) % F-L
        plot(l0/half_s_len_norm, F0);
        xlabel('\Delta length (l_{opt})')
        ylabel('F (F_0)')
    elseif(cond_itr==2) % F-v
        plot(l0/half_s_len_norm*10, F0);
        xlabel('velocity (l_{opt}/s)')
        ylabel('F (F_0)')
    else % F-pCa
        plot(pCa_list, F0);
        xlabel('pCa')
        ylabel('F (F_0)')
    end
    subplot(2,2,1)
    ylabel('\Delta length (l_{opt})')
    subplot(2,2,3)
    ylabel('F (F_0)')
    xlabel('time (s)')
end


