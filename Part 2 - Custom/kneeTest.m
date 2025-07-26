%% initialize simulation parameters 
clear; clc; 

addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

parms.forcible_detachment = 0;
% parms.kse = parms.kse*0.25;
parms.no_tendon = 1;

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

% parms.f = parms.f*0.5;
parms.JF = parms.JF*2;
parms.k11 = parms.k11 * 0.15;
parms.k12 = parms.k12 * 0.3;
parms.k21 = parms.k21 * 0.22;
parms.k22 = parms.k22 * 0.3;

nbins = 800;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0.0909;
parms.xss(end-1) = -40;
parms.xss(end) = -40;
parms.xi0 = linspace(10,60,nbins);
parms.xi = parms.xi0;
parms0 = parms;

%% run sumulation isometric and pre-move conditions
model = @fiber_dynamics;
dt = 0.005;
odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;

fig = figure;
pCa = 7.2; % <<-- change here to try different activation conditions

Ca = 10^(-pCa+6);
fig.UserData.parms = parms;
for k=1:2
    parms = parms0;
    x0 = parms.xss;

    % run prescribed motion (or hold) prior to knee drop
    if(k==1) % << isometric during first 1 sec
        prescribed_a = spline([-1,0],[0,0]); 
    else  % prescribed pre-movement.
        temp_t = -1:0.001:0;
        prescribed_a = spline(temp_t,-cos(temp_t/1*pi*2)*pi^2/2*pi*2); 
    end
    % run differential equation solver to get XB distributions
    [t_sim1,x_sim1] = ode15s(@(t,x)dAllStates_prescribed(t,x,model,parms,Ca,prescribed_a), ...
        (-1:dt:0),[pi/2, 0, x0],odeopt);    
    % calculate XB force
    F_sim1 = nan(height(x_sim1),1);
    for kk = 1:height(x_sim1)
        parms.vmtc = -x_sim1(kk,2)*10;
        [~,F_sim1(kk)] = model(t_sim1(kk), x_sim1(kk,3:end)', parms, Ca);
    end

    x0 = x_sim1(end,:); % << limb state at onset of drop

    % run knee drop simulation
    [t_sim2,x_sim2] = ode15s(@(t,x)dAllStates(t,x,model,parms,Ca), ...
        0:dt:10, x0,odeopt);
    % calculate XB force
    F_sim2 = nan(height(x_sim2),1);
    for kk = 1:height(x_sim2)
        parms.vmtc = -x_sim2(kk,2)*10;
        [~,F_sim2(kk)] = model(t_sim2(kk), x_sim2(kk,3:end)', parms, Ca);
    end

    subplot(3,1,1)
    hold on
    plot([t_sim1;t_sim2], [x_sim1(:,1);x_sim2(:,1)]*180/pi-90, 'LineWidth', 3-k)

    subplot(3,1,2)
    hold on
    plot([t_sim1;t_sim2], [F_sim1(:,1);F_sim2(:,1)], 'LineWidth', 3-k)

    if k==1
        fig.UserData.x_total_iso = [x_sim1;x_sim2];
        fig.UserData.t_total_iso = [t_sim1;t_sim2];
    else 
        fig.UserData.x_total_pre = [x_sim1;x_sim2];
        fig.UserData.t_total_pre = [t_sim1;t_sim2];
    end
end
xlim([-1 10])

ax1 = subplot(3,1,1);
ylabel('angle (deg)'), xlabel('time (s)')
fig.UserData.A_instance = plot([0,0], [-120 20]);
set(ax1,'ButtonDownFcn', ...
    @(s,e)update_time(s,e), 'HitTest','on')

ax2 = subplot(3,1,2);
ylabel('force'), xlabel('time (s)')
hold on
fig.UserData.F_instance = plot([0,0], [0 0.2]);
set(ax2,'ButtonDownFcn', ...
    @(s,e)update_time(s,e), 'HitTest','on')

linkaxes([ax1, ax2],'x');
xlim([-1 10])

fig.UserData.ax_dist = subplot(3,1,3);
fig.UserData.iso_dist = plot(nan,nan);
hold on
fig.UserData.pre_dist = plot(nan,nan);
xlabel('\Delta x')
ylabel('distribution of bound XB');
legend('isometric', 'pre-move')

%% differential equation for prescribed movement 
function Xd = dAllStates_prescribed(t, X, muscle_model, parms, Ca, prescribed_a)
% X(1): mass position, X(2): mass velocity, X(3:end): muscle states
parms.vmtc = -X(2)*10;
[Xmusd,~] = muscle_model(t, X(3:end), parms, Ca);

% mass dynamics % 
Xd = [X(2); ppval(prescribed_a,t); Xmusd];
end

%% differential equation for unprescribed movement 
function Xd = dAllStates(t, X, muscle_model, parms, Ca)
% X(1): mass position, X(2): mass velocity, X(3:end): muscle states

parms.vmtc = -X(2)*10;
[Xmusd,F_mus] = muscle_model(t, X(3:end), parms, Ca);
F_ext = -9.81*sin(X(1))*3.5 - 0.8*X(2); 

% mass dynamics % 
Xd = [X(2); (F_ext+F_mus*150); Xmusd];
end


function [] = update_time(src,eventData)
% load current simulation result that is shown in the figure
fig = src.Parent;
t_total = fig.UserData.t_total_iso;
x_total_iso = fig.UserData.x_total_iso;
x_total_pre = fig.UserData.x_total_pre;

% convert user's click to an instance of a simulation 
coords = eventData.IntersectionPoint;
t = min(max(-1,coords(1)),t_total(end));
idx = find(t_total>=t,1,'first');

% update XB distribution subplot to the selected instance 
set(fig.UserData.iso_dist, 'ydata', x_total_iso(idx,4:end-3), ...
    'xdata', fig.UserData.parms.xi0 + x_total_iso(idx,end));

set(fig.UserData.pre_dist, 'ydata', x_total_pre(idx,4:end-3), ...
    'xdata', fig.UserData.parms.xi0 + x_total_pre(idx,end));

set(fig.UserData.F_instance, 'xdata', [1 1]*t_total(idx))
set(fig.UserData.A_instance, 'xdata', [1 1]*t_total(idx))

set(fig.UserData.ax_dist, 'ylim', [0 max(max(x_total_iso(:,4:end-3)))], 'xlim', [-10 10]);

end