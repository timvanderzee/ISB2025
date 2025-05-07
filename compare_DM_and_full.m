clear all; clc; close all

load('parms.mat')
load('protocol.mat')

nbins = linspace(10,500,10);

for j = 1:(length(nbins)+1)
    
    if j == 1
        model = eval('@ripping_model_func_exp');
        parms.xss = zeros(1,8);
        ls = '-';
    else
        parms.nbins = round(nbins(j-1));
        parms.xi0 = linspace(-30,30,parms.nbins);
        
        model = eval('@ripping_model_func_exp_full');
        parms.xss = zeros(1, parms.nbins + 5);
        ls = '--';
       
    end
    
    parms.xss(end-2) = 0.0909;
    parms.xss(1) = 0;

    % simulate multiple pCas
    tic
    [Fmodel,~,~,Y] = stretch_shorten_multiple(model, XData, parms);
    pt(j) = toc;
    
    Data(j).Fmodel = Fmodel;
    Data(j).tmodel = Y.texp;
end

%% plot
close all
figure(1)
color = get(gca,'colororder');
pcolors = parula(length(nbins));
colors = [color(2,:); pcolors];

eps = nan(1,length(Data));
for i = 1:length(Data)
    eps(i) = mean(abs(Data(i).Fmodel - Data(end).Fmodel));
end

for j = length(nbins):-1:1

figure(1)
subplot(131)
plot(Data(j+1).tmodel, Data(j+1).Fmodel,'linewidth',2,'color',pcolors(j,:)); hold on

subplot(132)
bar(nbins(j), pt(j+1),'facecolor',pcolors(j,:),'barwidth',mean(diff(nbins)),'edgecolor','none'); hold on

h=subplot(133);
bar(nbins(j), eps(j+1),'facecolor',pcolors(j,:),'barwidth',mean(diff(nbins)),'edgecolor','none'); hold on
set(h,'YScale','log');
end


subplot(131);
plot(Data(j).tmodel, Data(j).Fmodel,'--','linewidth',2,'color',color(2,:)); hold on
xlabel('Time (s)')
ylabel('Force (-)')
box off
xlim([-.5 .5])
title('Forces')

subplot(132)
yline(pt(1),'--','color',color(2,:),'linewidth',2)
box off
xlabel('# bins')
ylabel('Processing time (s)')
title('Computational cost')

subplot(133)
yline(eps(1),'--','color',color(2,:),'linewidth',2)
box off
xlabel('# bins')
ylabel('Error (-)')
title('Error')

