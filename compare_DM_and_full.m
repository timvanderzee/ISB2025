clear all; clc; close all

load('parms.mat')
load('protocol.mat')

parms.forcible_detachment = 0;

nbins = round(linspace(50,500,10));
% nbins = 500;

%% conditions
Ca = 10^(-XData.pCas+6);
[us, Ts] = get_usTs(XData.v(1,:), XData.AMPs(1,:), XData.tiso(1,:), XData.ISI(1,:), parms);

%% simulate

for j = 1:(length(nbins)+1)
   disp(j)

    if j == 1
        modelname = 'ripping_model_func_exp';
        parms.xss = zeros(1,8);
        parms.xss(end-2) = 0.0909;
    else
        modelname = 'ripping_model_func_exp_full';
        parms.xi0 = linspace(-15,15,nbins(j-1));
        parms.nbins = length(parms.xi0);
        
        parms.xss = zeros(1,parms.nbins + 5);
        parms.xss(end-2) = 0.0909;
    end

    model = eval(['@',modelname]);
    
    tic
    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);
    pt(j,1) = toc;
    
    F = nan(1,length(x));
    pt2(j).t = nan(1, length(x));
    
    for i = 1:length(x)
        tic
        [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
        pt2(j).t(i) = toc;
    end
    
    pt(j,2) = mean(pt2(j).t);
    spt(j) = std(pt2(j).t);

    parms.xss(1) = 0;
    
    Data(j).Fmodel = F;
    Data(j).tmodel = t - sum(Ts(1:4));
end


%% plot
close all
figure(1)
color = get(gca,'colororder');
pcolors = parula(length(nbins));
colors = [color(2,:); pcolors];

Fi = nan(length(XData.texp), length(Data));
for i = 1:length(Data)
    Fi(:,i) = interp1(Data(i).tmodel, Data(i).Fmodel, XData.texp);
end

eps = nan(1,length(Data));
for i = 1:length(Data)
    eps(i) = mean(abs(Fi(:,i)-Fi(:,end)));
end

for j = length(nbins):-1:1
    
    figure(1)
    subplot(221)
    plot(Data(j+1).tmodel, Data(j+1).Fmodel,'linewidth',2,'color',pcolors(j,:)); hold on
    
    subplot(222)
    bar(nbins(j), pt(j+1,1),'facecolor',pcolors(j,:),'barwidth',mean(diff(nbins)),'edgecolor','none'); hold on
    
    subplot(223)
    bar(nbins(j), pt(j+1,2),'facecolor',pcolors(j,:),'barwidth',mean(diff(nbins)),'edgecolor','none'); hold on
    
    h=subplot(224);
    bar(nbins(j), eps(j+1),'facecolor',pcolors(j,:),'barwidth',mean(diff(nbins)),'edgecolor','none'); hold on
    set(h,'YScale','log');
end


subplot(221);
plot(Data(j).tmodel, Data(j).Fmodel,'--','linewidth',2,'color',color(2,:)); hold on
xlabel('Time (s)')
ylabel('Force (-)')
box off
xlim([-.5 .5])
title('Forces')

subplot(222)
yline(pt(1,1),'--','color',color(2,:),'linewidth',2)
box off
xlabel('# bins')
ylabel('Processing time (s)')
title('Computational cost')

subplot(223)
yline(pt(1,2),'--','color',color(2,:),'linewidth',2)
box off
xlabel('# bins')
ylabel('Processing time (s)')
title('Computational cost')

subplot(224)
yline(eps(1),'--','color',color(2,:),'linewidth',2)
box off
xlabel('# bins')
ylabel('Error (-)')
title('Error')

