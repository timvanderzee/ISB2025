clear all; close all; clc

load('parms.mat')
load('protocol.mat')

parms.forcible_detachment = 0;

nbins = round(linspace(50,500,10));
% nbins = 500;

%% conditions
Ca = 10^(-XData.pCas+6);
[us, Ts] = get_usTs(XData.v(1,:), XData.AMPs(1,:), XData.tiso(1,:), XData.ISI(1,:), parms);

%% simulate
ls = {'-','--'};
close all

for j = 1:2
    if j == 1
        modelname = 'ripping_model_func_exp';
    else
        modelname = 'ripping_model_func_exp_v2';
    end
    
    
    
    parms.xss = zeros(1,8);
    parms.xss(end-2) = 0.0909;


    model = eval(['@',modelname]);

    tic
    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);


    F = nan(1,length(x));


    for i = 1:length(x)

        [~,F(i)] = model(t(i), x(i,:)', parms, Ca);

    end

    figure(1)
    plot(t,F, ls{j}); hold on
end