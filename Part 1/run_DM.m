clear all; close all; clc

load('parms.mat', 'parms')
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

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

modelname = 'fiber_dynamics';

for j = 1

    if j == 1
        parms.xss = zeros(1,7);

    else
        nbins = 500;
        parms.xi0 = linspace(-15,15,nbins);
        parms.nbins = length(parms.xi0);
        parms.xss = zeros(1,parms.nbins + 4);
        parms.xss(end-2) = 0.0909;
    end
    
    model = eval(['@',modelname]);

    tic
    [t,x] = stretch_shorten(model, Ts, us, parms.xss, parms, Ca);

    F = nan(1,length(x));
    for i = 1:length(x)
        [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
    end

    figure(1)
    plot(t,F, ls{j}); hold on
    xlim([9.8 10.6])
end

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

    