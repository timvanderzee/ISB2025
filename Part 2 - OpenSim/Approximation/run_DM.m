clear all; close all; clc

% from Horslen paper: from 1 typical example fiber
cd('C:\Users\timvd\Documents\ISB2025\Part 1\parameters')
cd('C:\Users\timvd\Documents\ISB2025\Part 2 - OpenSim\Movement simulation\input\common')
load('parms.mat', 'parms')

% for a typical Horslen protocol
load('protocol.mat', 'XData')

parms.forcible_detachment = 0;
parms.Fscale = 1.5;

nbins = round(linspace(50,500,10));
% nbins = 500;

%% conditions
Ca = 10^(-XData.pCas+6);
[us, Ts] = get_usTs(XData.v(1,:), XData.AMPs(1,:), XData.tiso(1,:), XData.ISI(1,:), parms);

%% simulate
ls = {'-','--'};
close all

parms.act = 1; % active muscle volume
parms.cosa = 1; % cosine of pennation angle
parms.Noverlap = 1; % myofilament overlap
parms.n_func = @(xi, Q, eps) Q(1) ./ (sqrt(2*pi)*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))) * exp(-((xi-(Q(2)/max(Q(1), eps))).^2) / (2*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))^2)); 

model = @fiber_dynamics;
ti = linspace(9.9, 10.3, 100);
ni = nan(length(ti), 500, 2);
Fi = nan(length(ti),2);
li = nan(length(ti),2);
Li = nan(length(ti),2);

nbins = 500; % number of bins
    
for j = 1:2
    parms.xi = linspace(-15,15,nbins); % initial strain vector (power stroke)

    if j == 1 % DM
        parms.xss = zeros(1,7);

    else % discretized

        parms.nbins = length(parms.xi);
        parms.xss = zeros(1,parms.nbins + 4); % 4 non-cross-bridge states
        parms.xss(end-2) = 0.0909; % DRX state, given default parameters
    end
    
    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);

    L = x(:,end-1);
    lce = x(:,end);
    F = nan(1,length(x));
    n = nan(length(parms.xi), length(x));
    for i = 1:length(x)
        [~,F(i), ~, n(:,i)] = model(t(i), x(i,:)', parms, Ca);
    end
    
    if j == 1
        xi = parms.xi;
    else
        xi = parms.xi + repmat((lce(:) - parms.lce0), 1, length(parms.xi));
    end
    
    % interpolate n
    ni(:,:,j) = interp1(t, n', ti);
    li(:,j) = interp1(t, lce, ti);
    Li(:,j) = interp1(t, L, ti);
    Fi(:,j) = interp1(t, F, ti);
    
%     figure(1)
%     subplot(121)
%     plot(t,F, ls{j}, 'linewidth',2); hold on
%     xlim([9.8 10.6])
    
%     subplot(122)
%     plot(xi(end,:), n(:,end)', ls{j}); hold on
end

%% compare distributions
% interpolate n
close all
id = 89;

figure(1)
color = get(gca,'colororder');
subplot(222);
plot(ti-9.9, Li(:,2), '-', 'linewidth',2, 'color', color(1,:));
xlim([0 .4])
box off
% xlabel('Time (s)')
ylabel('\Delta Length (L_0)')
title('Length trajectory')
y1 = line('xdata', [ti(id) ti(id)]-9.9, 'ydata', [0 5], 'linestyle','--', 'color', [0 0 0]);
ylim([0 5])

subplot(224)
plot(ti-9.9, Fi(:,1), '--', 'linewidth',2); hold on
plot(ti-9.9, Fi(:,2), '-', 'linewidth',2, 'color', color(1,:));
xlim([0 .4])

y2 = line('xdata', [ti(id) ti(id)]-9.9, 'ydata', [0 1.8], 'linestyle','--', 'color', [0 0 0]);
ylim([0 1.8])
box off
xlabel('Time (s)')
ylabel('Force (F_0)')
title('Force trajectory')
% legend('Approximated','Original','location','best')
% legend boxoff

subplot(2,2,[1 3])
h = line('xdata',parms.xi * parms.h * 1e9, 'ydata', ni(id,:,1), 'linestyle','--', 'color', color(1,:),'linewidth',2);
g = line('xdata', (parms.xi + li(id,2)) * parms.h * 1e9, 'ydata', ni(id,:,2), 'linestyle','-', 'color', color(1,:),'linewidth',2);
axis([-24 24 0 2])
box off
xlabel('Cross-bridge strain (nm)')
ylabel('Fraction attached')
title('Cross-bridge distribution')
set(gcf,'units','normalized','position',[.2 .1 .3 .4])

%%
figure(1)
clear cframe;

for i = 1:size(ni,1)
    set(y1, 'xdata', [ti(i) ti(i)]-9.9);
    set(y2, 'xdata', [ti(i) ti(i)]-9.9);
    
    set(h, 'ydata', ni(i,:,1));
    set(g, 'xdata', (parms.xi + li(i,2)) * parms.h * 1e9, 'ydata', ni(i,:,2));
    drawnow
    
    cframe(i) = getframe(gcf);
end

%%
close all
writerObj=VideoWriter('distributions_realtime.mp4','MPEG-4');
writerObj.FrameRate = round(1/mean(diff(ti)));
open(writerObj);
writeVideo(writerObj,cframe);
close(writerObj);

%%
writerObj=VideoWriter('distributions_slowed_10x_with_approx.mp4','MPEG-4');
writerObj.FrameRate = round(1/mean(diff(ti))) / 10;
open(writerObj);
writeVideo(writerObj,cframe);
close(writerObj);



%% test effect of Ca
% simulate isometric at different calcium levels
pCas = [9 linspace(7, 4.5, 100)];
Cas = 10.^(-pCas+6);
parms.vmtc = 0;
F = nan(size(Cas));

for j = 1:length(Cas)
    
    Ca = Cas(j);
    [t,x] = ode15s(model, [0 10], parms.xss, [], parms, Ca);   
    
    [~,F(j)] = model(t(end), x(end,:)', parms, Ca);
        
end

figure(2)
plot(-pCas, F)

%% test overlap and activation effects
close all
c = linspace(0,1,5);
model = eval(['@',modelname]);
   
for k = 1:length(c)
%     parms.Noverlap = c(k);
    parms.act = c(k);

    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);
    F = nan(1,length(x));
        
    for i = 1:length(x)
        [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
    end
    
    figure(1)
    plot(t,F); hold on
    xlim([9.8 10.6])
end

    