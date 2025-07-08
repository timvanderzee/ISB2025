clear all; close all; clc

% from Horslen paper: from 1 typical example fiber
load('parms.mat', 'parms')

% for a typical Horslen protocol
load('protocol.mat', 'XData')

parms.forcible_detachment = 0;

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
model = @fiber_dynamics;
parms.n_func = @(xi, Q, eps) Q(1) ./ (sqrt(2*pi)*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))) * exp(-((xi-(Q(2)/max(Q(1), eps))).^2) / (2*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))^2)); 

ti = linspace(9.8, 10.6, 100);
ni = nan(length(ti), 500, 2);
Fi = nan(length(ti),2);
li = nan(length(ti),2);

for j = 1:2

    if j == 1 % DM
        parms.xss = zeros(1,7);

    else % discretized
        nbins = 500; % number of bins
        parms.xi = linspace(-15,15,nbins); % initial strain vector (power stroke)
        parms.nbins = length(parms.xi);
        parms.xss = zeros(1,parms.nbins + 4); % 4 non-cross-bridge states
        parms.xss(end-2) = 0.0909; % DRX state, given default parameters
    end
    
    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);

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
id = 57;

figure(1)
color = get(gca,'colororder');
subplot(121);

plot(ti, Fi(:,1), '-', 'linewidth',2); hold on
plot(ti, Fi(:,2), '--', 'linewidth',2);
xlim([9.8 10.6])

y = line('xdata', [ti(id) ti(id)], 'ydata', [0 1.8], 'linestyle','--', 'color', [0 0 0]);
ylim([0 1.8])
box off
xlabel('Time (s)')
ylabel('Force (-)')
title('Force trajectory')
legend('Approximated','Original','location','best')
legend boxoff

subplot(122)
h = line('xdata',parms.xi * parms.h * 1e9, 'ydata', ni(id,:,1), 'linestyle','-', 'color', color(1,:),'linewidth',2);
g = line('xdata', (parms.xi + li(id,2)) * parms.h * 1e9, 'ydata', ni(id,:,2), 'linestyle','--', 'color', color(2,:),'linewidth',2);
axis([-24 24 0 2])
box off
xlabel('XB strain (nm)')
ylabel('Fraction of attached XBs')
title('Cross-bridge (XB) distribution')

set(gcf,'units','normalized','position',[.2 .3 .5 .4])

%%
figure(1)
clear cframe;

for i = 1:length(ni)
    set(y, 'xdata', [ti(i) ti(i)]);
    
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


writerObj=VideoWriter('distributions_slowed_10x.mp4','MPEG-4');
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

    