clear all; clc; close all

load('parms.mat')

parms.forcible_detachment = 0;

parms.d = .001;

vs = [-100 0 100];

figure(1)
color = get(gca,'colororder');
ls = {'-','--'};

pt = nan(2,3);
n = nan(2,3);

for i = 1:length(vs)
    parms.vmtc = vs(i);

for j = 1:2
    if j == 1
        x0 = zeros(1,5);
        modelname = 'twostate_model_func_exp_v2';
        
    else
        x0 = zeros(1,4);
        modelname = 'twostate_model_func_exp';
    end
    
    model = eval(['@',modelname]);

%     tic
    [t,x] = ode15s(model, [0 1], x0(:), [], parms, 1);

%     n(j,i) = length(t);
%     pt(j,i) = toc;
    
    figure(1)
    subplot(221)
    plot(t,x(:,1),ls{j},'color',color(i,:)); hold on

    subplot(222)
    plot(t,x(:,2),ls{j},'color',color(i,:));hold on

    subplot(223)
    plot(t,x(:,3),ls{j},'color',color(i,:));hold on

    if j == 1
        subplot(224)
        plot(t,x(:,4) * parms.d,ls{j},'color',color(i,:));hold on
    end
    
end
end
