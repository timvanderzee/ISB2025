clear all; clc; close all

load('parms.mat')
load('protocol.mat')

ls = {'-','--'};

for j = 2:2
    
    if j == 1
        model = eval('@ripping_model_func_exp');
        parms.xss = zeros(1,8);
    else
        model = eval('@ripping_model_func_exp_full');
        parms.xss = zeros(1, parms.nbins + 5);
    end
    
    parms.xss(end-2) = 0.0909;
    parms.xss(1) = 0;

    % simulate multiple pCas
    tic
    [Fmodel,~,~,Y] = stretch_shorten_multiple(model, XData, parms);
    pt(j) = toc;
    
    figure(1)
    subplot(121)
    plot(Y.texp, Fmodel,ls{j},'linewidth',2); hold on
end


%% make nice
figure(1)

subplot(121);
xlabel('Time (s)')
ylabel('Force (-)')
box off
legend('DM','Discretized','location','best')
legend boxoff
xlim([-.5 .5])
title('Forces')

subplot(122)
bar(pt)
box off
ylabel('Processing time (s)')
title('Computational cost')