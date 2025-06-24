% clear; clc; close all

addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;
pCa = 4.5;
Ca = 10^(-pCa+6);

parms.forcible_detachment = 0;
parms.kse = parms.kse*0.5;
parms.kpe = 0;
parms.no_tendon = 0;

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

%%

model = @fiber_dynamics;
% parms.xss = zeros(1,7);

parms.xi0 = linspace(-50,50,nbins);
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);

%%

parms.xss(end-2) = 0.0909;
x0 = parms.xss;
tic
[t_sim,x_sim] = ode15s(@(t,x)dAllStates(t,x,model,parms,Ca), [0:0.01:2], ...
    [0, 0, x0],odeopt);

F_sim = nan(height(x_sim),1);
for i = 1:height(x_sim)
    [~,F_sim(i)] = ...
        model(t_sim(i), x_sim(i,3:end)', parms, Ca);
end
toc
%%
% figure
% for i = 1:height(x_sim)
%     plot(x_sim (i, 4:nbins+4));
%     ylim([0 2])
%     pause(0.2)
% end

%%

sim_fig = figure;
hold on

Lce0 = 1; % * half_sarc_len
Lse0 = 0.4;

sim_fig.UserData.ce_x = [-pi:0.01:pi, pi:-0.01:-pi, -pi]/pi/2+0.5;
sim_fig.UserData.ce_y = [sech(-pi:0.01:pi) -sech(pi:-0.01:-pi) sech(-pi)]*0.1;
sim_fig.UserData.ce_handle = fill(sim_fig.UserData.ce_x*Lce0, ...
    sim_fig.UserData.ce_y, [1 .3 .5], 'LineStyle','none');

sim_fig.UserData.spring_x = [0 2 2.5:1:7.5 8 10]/10;
sim_fig.UserData.se_handle = plot(Lce0+sim_fig.UserData.spring_x*Lse0, ...
    [0 0 1 -1 1 -1 1 -1 0 0]*0.02, '-', 'lineWidth', 4, 'color', [.8 .7 .5]);

sim_fig.UserData.mass_handle = plot((Lce0+Lse0), 0, '.', ...
    'color', [.1 .5 .8], 'markerSize', 50);

sim_fig.UserData.Lce0 = Lce0;
sim_fig.UserData.Lse0 = Lse0;

text(1.3, 0.2, char(8594), 'fontSize', 40)
text(1.3, 0.25, 'const. F_{ext}')

xlim([0 Lce0+Lse0*1.5])
axis equal
pause(0.05)
yticks('');
xlabel('displacement (l_{opt})')


for k=1:height(x_sim)
    update_figure(sim_fig, Lce0+x_sim(k,end)/half_s_len_norm, ...
        (Lce0+Lse0)+x_sim(k,end-1)/half_s_len_norm, ...
        (Lce0+Lse0)+x_sim(k,1)/half_s_len_norm)
    xlim([0 Lce0+Lse0*1.5])
    axis equal
    pause(0.01)
end

%%

function [] = update_figure(fig, lce, lmtc, lmass)
set(fig.UserData.ce_handle, 'Xdata', lce*fig.UserData.ce_x);
set(fig.UserData.ce_handle, 'Ydata', ...
    fig.UserData.ce_y/sqrt((lce)));

set(fig.UserData.se_handle, 'Xdata', lce + (lmtc-lce)*fig.UserData.spring_x);
set(fig.UserData.mass_handle, 'Xdata', lmass);
end

function Xd = dAllStates(t, X, model, parms, Ca)

% X(1): mass position, X(2): mass velocity, X(3:end): muscle states

parms.vmtc = X(2);
[Xmusd,F_mus] = model(t, X(3:end), parms, Ca);
F_ext = 0.5; 

% mass dynamics % 
Xd = [X(2); (F_ext-F_mus)*100; Xmusd];

end