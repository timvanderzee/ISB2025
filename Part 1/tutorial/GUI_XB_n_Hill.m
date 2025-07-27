%% initialize - load parameters and set up 
clear;
addpath(genpath([pwd,'/../..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

half_s_len_norm = parms.s/2/parms.h;
nbins = 600;
pCa = 4.5;
Ca = 10^(-pCa+6);

% parms.forcible_detachment = 0;
% parms.kse = 0;
% parms.kpe = 0;
% parms.no_tendon = 1;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.xi0 = linspace(-10,10,nbins);
parms.xi = parms.xi0;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0.0909;
parms.pCa_ran = 4.5;

model = @fiber_dynamics;
parms0 = parms;

%% set up figures for GUI
fig = figure;

% Protocol data will be passed as fig.UserData between subplots.
% We start making figures using a protocol. 
fig.UserData.parms = parms;
fig.UserData.protocol_duration = [0.05,0.1,0.1,0.1,0.1,0.1]; % should be integer multiplication of dt=0.001
fig.UserData.protocol_v = [0, 0,-10, 10, 10, -10];
fig.UserData.protocol_pCa = [9 [1 1 1 1 1]*parms.pCa_ran];
load('hill_properties_pCa_45.mat');
fig.UserData.hill_properties = hill_properties;
sim_XB(fig, model,parms);

% set GUI range for attachment/detachment parameters  
parm_range = [100, 1, 30, 3, 30, 3, 4.5];

% run button - if clicked, GUI will run XB and Hill simulation 
ax = subplot(5,18,11:12);
fig.UserData.ax_run = ax;
set(fig.UserData.ax_run,'ButtonDownFcn', ...
    @(s,e)run_sim_callback(s,e,model),...
    'HitTest','on');
fig.UserData.statusHandle = text(0.5, 0.5, 'RUN',...
    'HorizontalAlignment','center','PickableParts','none');
xlim([0 1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.Box = 'off';

% generate Hill button - if clicked, GUI will run F-l, F-v, F-pCa 
% simulation and update Hill type model properties. 
ax = subplot(13,6,6);
fig.UserData.ax_hill = ax;
set(fig.UserData.ax_hill,'ButtonDownFcn', ...
    @(s,e)sim_flfvfpca(s,[],fig.UserData.parms),...
    'HitTest','on');
fig.UserData.hillHandle = text(0.5, 0.5, 'generate Hill',...
    'HorizontalAlignment','center','PickableParts','none');
xlim([0 1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.Box = 'off';

% GUI to chose stretch/shorten protocol
fig.UserData.ax_protocol = subplot(5,9,[4,5]);
plot(1,0,'.','markerSize',20);
title('protocol')
xlim([0 3]);
set(fig.UserData.ax_protocol,'ButtonDownFcn', ...
    @(s,e)len_protocol_callback(s,e),...
    'HitTest','on');
xticks(0:3)
xticklabels({'shorten','lengthen', ...
    sprintf('s cycle'),sprintf('l cycle')})
fig.UserData.ax_protocol.YAxis.Visible = 'off';

% GUI to adjust attachment/detachment parameters 
ax = subplot(5,3,1);
fig.UserData.ax_parm = ax;
plot(ax,[0:1, 3:6, 8],...
    [parms.f,parms.w,parms.k11,parms.k12,parms.k21,parms.k22,(pCa-4.5)]./parm_range,...
    '.', 'markerSize', 20, 'PickableParts', 'none');
xlim([-1 9])
xticks([0:1, 3:6, 8])
xticklabels({'f','w','k_{11}','k_{12}','k_{21}','k_{22}','pCa'});
text(ax, 0, 1.3, 'attachment parms', 'fontSize', 10, 'HorizontalAlignment','center');
text(ax, 4.5, 1.3, 'detachment parms', 'fontSize', 10, 'HorizontalAlignment','center');
text(ax, 8, 1.3, 'activation', 'fontSize', 10, 'HorizontalAlignment','center');

% subplot showing attachment function based on selected parameters
fig.UserData.ax_ffunc = subplot(5,3,4);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('attachment rate (f)')

% subplot showing detachment function based on selected parameters
fig.UserData.ax_gfunc = subplot(5,3,7);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('detachment rate (g)')

% subplot showing protocol-imposed pCa change 
fig.UserData.ax_pCa = subplot(5,3,8);
plot(fig.UserData.plot_t_protocol, fig.UserData.plot_pCa_protocol, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [4 9.5], 'PickableParts', 'none');
xlabel('time (s)')
ylabel('pCa')

% subplot showing protocol-imposed length change 
fig.UserData.ax_len = subplot(5,3,5);
plot(nan, nan, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [-10 10], 'PickableParts', 'none');
xlabel('time (s)')
ylabel('\Delta length (ps)')

% subplot showing simulated force outcome 
fig.UserData.ax_F = subplot(5,3,[11,14]);
hxb = plot(fig.UserData.t_total, fig.UserData.F_total, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [0 max(fig.UserData.F_total)*3], 'PickableParts', 'none');
hhill = plot(fig.UserData.t_total, fig.UserData.hillF_total, 'PickableParts', 'none');
legend([hxb hhill], 'XB', 'Hill', 'Location', 'southeast')
xlabel('time (s)')
ylabel('force (F_0)')

% subplot showing simulated cross-bridge distributions
fig.UserData.ax_dist = subplot(5,3,[10,13]);
plot_k = 100;
ax_dist = plot(fig.UserData.parms.xi0, ...
    fig.UserData.x_total(plot_k, 1:fig.UserData.parms.nbins));
xlabel('\Delta x (ps)')
ylabel('distribution of bound XB')

% allow user to click time-{imposed length, imposed pCa, force} subplots to
% see cross-bridge distribution at chosen time 
set(fig.UserData.ax_len,'ButtonDownFcn', ...
    @(s,e)update_time(s,e), 'HitTest','on')
set(fig.UserData.ax_F,'ButtonDownFcn', ...
    @(s,e)update_time(s,e), 'HitTest','on')
set(fig.UserData.ax_pCa,'ButtonDownFcn', ...
    @(s,e)update_time(s,e), 'HitTest','on')

% allow user to click ax_parm subplot to adjust attachment/dettachment
% parameters. After changing parameters, simulated force will disappear
% and user will need to "Run" again to get updated force. 
set(fig.UserData.ax_parm,'ButtonDownFcn', ...
    @(s,e)update_parms(s,e,parm_range),...
    'HitTest','on');

% initialize some plots using existing figure callbacks 
dummyE.IntersectionPoint = [1 -100];
len_protocol_callback(fig.UserData.ax_len, dummyE);

dummyE.IntersectionPoint = [0.1 0];
update_time(fig.UserData.ax_len, dummyE)

dummyE.IntersectionPoint = [-2,0];
update_parms(fig.UserData.ax_parm,dummyE,parm_range);

linkaxes([fig.UserData.ax_len, fig.UserData.ax_F, fig.UserData.ax_pCa],'x');
run_sim_callback(fig.UserData.ax_run, [], model);

%% Callback function when user clicks "Run" button. 
function [] = run_sim_callback(src,~,model)
fig = src.Parent;
set(fig.UserData.statusHandle, 'string', sprintf('running \n please wait'));
pause(0.01);
sim_XB(fig,model,fig.UserData.parms); % <<-- XB (and Hill) simulation 
dummyE.IntersectionPoint = [10 0];
update_time(fig.UserData.ax_len,dummyE)

set(fig.UserData.statusHandle, 'string', 'RUN');
end

%% Function to simulate cross-bridge model. 
% It will also call "sim_Hill" at the end. 
function [] = sim_XB(fig,model,parms)

odeopt = odeset('maxstep',1e-2);
x0 = parms.xss;
dt = 0.001;

% load selected protocol 
protocol_duration = fig.UserData.protocol_duration;
protocol_v = fig.UserData.protocol_v;
protocol_pCa = fig.UserData.protocol_pCa;

% initialize variables to save simulation. 
t_total = nan(round(sum(protocol_duration)/dt)+length(protocol_duration),1);
x_total = nan(length(t_total),length(x0));
v_total = nan(size(t_total));
F_total = nan(size(t_total));
pCa_total = nan(size(t_total));

% variables defined for protocol plotting purpose 
plot_t_protocol = nan(length(protocol_pCa)*2, 1);
plot_pCa_protocol = nan(length(protocol_pCa)*2, 1);

% iterate through sub-protocols 
curr_idx = 1;
curr_time = 0;
for p_itr = 1:length(protocol_v)
    plot_t_protocol((-1:0)+p_itr*2) = [curr_time, curr_time+protocol_duration(p_itr)];
    plot_pCa_protocol((-1:0)+p_itr*2) = protocol_pCa(p_itr);

    % set up sub-protocl using defined velocity and pCa. 
    parms.vmtc = protocol_v(p_itr);
    Ca = 10.^(-protocol_pCa(p_itr)+6);

    % solve differential equation 
    [t_sim,x_sim] = ode15s(model, ...
        curr_time:dt:curr_time+protocol_duration(p_itr), x0, odeopt, ...
        parms, Ca);

    % calculate force using ODE output 
    F_sim = nan(size(t_sim));
    for k_itr = 1:length(t_sim)
        [~,F_sim(k_itr)] = model(t_sim(k_itr), x_sim(k_itr,:)', parms, Ca);
    end

    % save sub-protocol results
    t_total(curr_idx:curr_idx+length(t_sim)-1) = t_sim;
    pCa_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_pCa(p_itr);
    v_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_v(p_itr);
    F_total(curr_idx:curr_idx+length(t_sim)-1) = F_sim;
    x_total(curr_idx:curr_idx+length(t_sim)-1,:) = x_sim;

    % prepare for the next iteration 
    x0 = x_sim(end,:);
    curr_idx = curr_idx + length(t_sim);
    curr_time = t_sim(end);
end

% pass simulation outcome through fig.UserData 
fig.UserData.t_total = t_total;
fig.UserData.x_total = x_total;
fig.UserData.F_total = F_total;
fig.UserData.pCa_total = pCa_total;
fig.UserData.v_total = v_total;
fig.UserData.plot_t_protocol = plot_t_protocol;
fig.UserData.plot_pCa_protocol = plot_pCa_protocol;

% call Hill type simulation 
sim_Hill(fig,fig.UserData.hill_properties,...
    parms);
end

function sim_Hill(fig,hill_p,parms)
% approximate force using Hill-type model. 
% This function will use Hill type model properties generated by clicking 
% "generate Hill" button, or use loaded properties if it hasn't been
% generated yet. Default Hill properties corresponds to default XB properties. 
% Derived for no tendon, no elastic elements for this example. 

half_s_len_norm = parms.s/2/parms.h;

% load Hill type model properties 
l_total = fig.UserData.x_total(:,end); 
pCa_total = fig.UserData.pCa_total;
v_total = fig.UserData.v_total;

% use spline interpolation to calculate force from F-l, F-v, F-pCa curves. 
fig.UserData.hillF_total = ppval(hill_p.FL_spline, l_total/half_s_len_norm).*...
    ppval(hill_p.FV_spline, -v_total/half_s_len_norm).*...
    ppval(hill_p.FPca_spline, pCa_total);
end

%% Callback that is called when "generate Hill" button is clicked. 
% It will run F-l, F-v, F-pCa protocols and save the outcome as splines. 
function [] = sim_flfvfpca(src,~,parms)
% load parameters and set up simulations 
fig = src.Parent;
set(fig.UserData.hillHandle, 'string', 'running - please wait');
odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
pause(0.01);
ref_pCa = 4.5; % << pCa during F-L, F-v simulation 

% variable to store hill properties 
hill_properties = struct();

for cond_itr = 1:3 % 1: F-L, 2: F-v, 3: F-pCa
   
    if(cond_itr==1) % F-L
        l0 = (-0.5:0.1:0.5) *half_s_len_norm;
    elseif(cond_itr==2) % F-v
        l0 = (-1:0.2:1) *half_s_len_norm / 20;
    else % F-pCa
        pCa_list = [4.5,5,5.5:0.1:6.5,7:0.5:9];
        l0 = zeros(size(pCa_list));
    end

    % subplot to monitor time-length 
    ax = subplot(6,6,cond_itr*12-1-6, 'colorOrder', winter(length(l0))); 
    cla(ax)
    hold on
    % subplot to monitor time-force
    ax = subplot(6,6,cond_itr*12-1, 'colorOrder', winter(length(l0))); 
    cla(ax)
    hold on
    F0 = nan(size(l0));

    % run the simulation 
    for k=1:length(l0)
        model = @fiber_dynamics;

        % activation protocol 
        if(cond_itr==3)
            pCa = pCa_list(k);
        else
            pCa = ref_pCa;
        end
        Ca = 10^(-pCa+6);

        % initial length 
        parms.xss(end-1) = l0(k);
        parms.xss(end) = l0(k);
        parms.lce0 = l0(k);

        % velocity protocol 
        if(cond_itr==1||cond_itr==3) % isometric for F-L & F-pCa
            us = [0,0];
            Ts = [2,1];
        else % isokinetic for F-v
            us = [0,-l0(k)*20,0];
            Ts = [3,0.05,3];
        end

        x0 = parms.xss;
        t = 0;
        x = x0;
       
        % solve ODE to get XB distributions over time 
        temp_idx = 1;
        for i = 1:length(us)
            parms.vmtc = us(i);
            T = Ts(i);
            [tnew,xnew] = ode15s(model, [0 T], x0, odeopt, parms, Ca);
            x0 = xnew(end,:);
            t = [t; tnew(2:end)+t(end)];
            x = [x; xnew(2:end,:)];
            if(i==2), temp_idx = height(x); end
        end
        
        % calculate force
        F = nan(1,height(x));
        for i = 1:height(x)
            [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
        end

        % plot time-length
        subplot(6,6,cond_itr*12-1-6)
        plot(t, x(:,end)/half_s_len_norm)
        hold on

        % plot time-force
        subplot(6,6,cond_itr*12-1)
        plot(t,F)
        hold on
        pause(0.1)
        F0(k) = F(temp_idx);
    end

    % plot F-L, F-v, F-pCa curves obtained from simulations and store 
    % the curves as splines 
    subplot(13,6,(-18:6:0)+cond_itr*24+6)
    if(cond_itr==1) % F-L
        plot(l0/half_s_len_norm, F0);
        xlabel('\Delta length (l_{opt})')
        ylabel('F (F_0)')
        title('F-l')
        hill_properties.FL_spline = spline(l0/half_s_len_norm, F0);
    elseif(cond_itr==2) % F-v
        plot(l0/half_s_len_norm/0.05, F0);
        xlabel('velocity (l_{opt}/s)')
        ylabel('F (F_0)')
        title('F-v')
        hill_properties.FV_spline = spline(l0/half_s_len_norm/0.05, F0/F0(abs(l0)<0.0001));
    else % F-pCa
        plot(pCa_list, F0);
        xlabel('pCa')
        ylabel('F (F_0)')
        title('F-pCa')
        xlim([4.5 9])
        hill_properties.FPca_spline = spline(pCa_list, F0/F0(1));
    end

    % label subplots 
    ax1 = subplot(6,6,cond_itr*12-1-6);
    ylabel('\Delta length (l_{opt})')
    ax2 = subplot(6,6,cond_itr*12-1);
    ylabel('F (F_0)')
    xlabel('time (s)')
    linkaxes([ax1, ax2],'x');
    if(cond_itr == 2)
        xlim([3, 3.1])
        plot(ones(size(F0))*3.05, F0, '.', 'markerSize', 10);
    end
end

% pass simulated Hill properties to other functions through fig.UserData
fig.UserData.hill_properties = hill_properties;

% get simulation protocol that is shown in the entire figure
l_total = fig.UserData.x_total(:,end); 
pCa_total = fig.UserData.pCa_total;
v_total = fig.UserData.v_total;
% apply new Hill type propetry we just generated, update Hill force 
fig.UserData.hillF_total = ppval(hill_properties.FL_spline, ...
    l_total/half_s_len_norm).*...
    ppval(hill_properties.FV_spline, -v_total/half_s_len_norm).*...
    ppval(hill_properties.FPca_spline, pCa_total);
set(fig.UserData.ax_F.Children(1), 'xdata', fig.UserData.t_total, ...
    'ydata', fig.UserData.hillF_total);

set(fig.UserData.hillHandle, 'string', 'generate Hill');
end

%% Callback to update attachment/detachment parameters 
function [] = update_parms(src,eventData, parm_range)
fig = src.Parent;

coords = eventData.IntersectionPoint;
x = coords(1);
y = max(0,min(1,coords(2)));

y_s = get(fig.UserData.ax_parm.Children(end),'ydata');

% depending on click location, update a paramter. 
if(x<-1)
else
    if(x<0.5) % f of attachment function
        fig.UserData.parms.f = y * parm_range(1);
        y_s(1) = y;
    elseif(x<1.5) % w of attachment function
        fig.UserData.parms.w = y * parm_range(2);
        y_s(2) = y;
    elseif(x<3.5) % k11 of detachment function
        fig.UserData.parms.k11 = y * parm_range(3);
        y_s(3) = y;
    elseif(x<4.5) % k12 of detachment function
        fig.UserData.parms.k12 = y * parm_range(4);
        y_s(4) = y;
    elseif(x<5.5) % k21 of detachment function
        fig.UserData.parms.k21 = y * parm_range(5);
        y_s(5) = y;
    elseif(x<6.5) % k22 of detachment function
        fig.UserData.parms.k22 = y * parm_range(6);
        y_s(6) = y;
    elseif(x<8.5) % pCa level during activation 
        fig.UserData.parms.pCa_ran = y * parm_range(7) + 4.5;
        y_s(7) = y;
    end

    % remove force and XB distribution figures as they are no longer valid
    % for the new sets of parameters 
    fig.UserData.x_total = fig.UserData.x_total*nan;
    fig.UserData.F_total = fig.UserData.F_total*nan;
    fig.UserData.hillF_total = fig.UserData.hillF_total*nan;
    % update using an existing callback 
    dummyE.IntersectionPoint = [0.1 0];
    update_time(fig.UserData.ax_len,dummyE)
end

% update GUI slider location based on click 
set(fig.UserData.ax_parm.Children(end),'ydata',y_s)
set(fig.UserData.ax_parm, 'ylim',[0 1]);

% get updated parameters 
parms = fig.UserData.parms;

% plot attachment rate function 
set(fig.UserData.ax_ffunc.Children(1), 'yData', ...
    parms.f_func(parms.xi, parms.f, parms.w));
set(fig.UserData.ax_ffunc, 'ylim', [0 300], 'xlim', [min(parms.xi0) max(parms.xi0)]);

% plot detachment rate function 
set(fig.UserData.ax_gfunc.Children(1), 'yData', ...
    parms.g_func(parms.xi, parms.k11, -parms.k12)+...
    parms.g_func(parms.xi, parms.k21, parms.k22));
set(fig.UserData.ax_gfunc, 'ylim', [0 6e5], 'xlim', [min(parms.xi0) max(parms.xi0)]);

% plot pCa protocol and pass it to other fucntions  
pCa_temp = get(fig.UserData.ax_pCa.Children(2), 'yData');
pCa_temp(3:end) = fig.UserData.parms.pCa_ran;
set(fig.UserData.ax_pCa.Children(2), 'yData', pCa_temp);
fig.UserData.protocol_pCa = [9 [1 1 1 1 1]*fig.UserData.parms.pCa_ran];
end

%% Callback to update XB distribution subplot based on where user clicks. 
% This function is also called when simulation is updated, so always plot
% simulated force and protocol. 
function [] = update_time(src,eventData)

% load current simulation result that is shown in the figure
fig = src.Parent;
t_total = fig.UserData.t_total;
x_total = fig.UserData.x_total;
F_total = fig.UserData.F_total;

% convert user's click to an instance of a simulation 
coords = eventData.IntersectionPoint;
t = min(max(0,coords(1)),t_total(end));
idx = find(t_total>=t,1,'first');

% update XB distribution subplot to the selected instance 
set(fig.UserData.ax_dist.Children(1), 'ydata', x_total(idx,2:end-3), ...
    'xdata', fig.UserData.parms.xi0 + x_total(idx,end));

if(~sum(isnan(x_total)))
    set(fig.UserData.ax_dist, 'ylim', [0 max(max(x_total(:,2:end-3)))], 'xlim', [-10 10]);
end

% update simulated force and protocol, show what time is selected. 
set(fig.UserData.ax_F.Children(1), 'ydata', fig.UserData.hillF_total);
set(fig.UserData.ax_F.Children(2), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_F.Children(3), 'ydata', F_total);
set(fig.UserData.ax_len.Children(1), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_len, 'ylim', [-1 1]*6);
set(fig.UserData.ax_pCa.Children(1), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_pCa, 'ylim', [4 9.5]);

if(~sum(isnan(fig.UserData.F_total)))
    set(fig.UserData.ax_F, 'ylim', ...
        [0 max([fig.UserData.F_total; fig.UserData.hillF_total])*1.1]);
end

end

%% callback to select simulation protocol. 
% If manually editing the code, make sure that protocol_t, protocol_v,
% protocol_pCa all have same number of elements. 
function [] = len_protocol_callback(src,eventData)
% load selections 
fig = src.Parent;
coords = eventData.IntersectionPoint;
x = coords(1);

% set protocol pCa. It is set to be minimally activated for the first
% sub-trial, and then activated to selected level during the rest of the
% sub-trials. User may manually change this if they want to simulate a
% different activation scenario. 
fig.UserData.protocol_pCa = [9 [1 1 1 1 1]*fig.UserData.parms.pCa_ran];

% User may select from pre-defined profiles. 
if(x<0.5) % hold-shorten-hold
    set(src.Children(1),'XData',0);
    fig.UserData.protocol_v = [0, 0, 0, -50, 0, 0];
elseif(x<1.5) % hold-stretch-hold
    set(src.Children(1),'XData',1);
    fig.UserData.protocol_v = [0, 0, 0, 50, 0, 0];
elseif(x<2.5) % hold-shorten-stretch-shorten-stretch-hold
    set(src.Children(1),'XData',2);
    fig.UserData.protocol_v = [0, 0, -50, 50, -50, 50];
else % hold-stretch-shorten-stretch-shorten-hold
    set(src.Children(1),'XData',3);
    fig.UserData.protocol_v = [0, 0, 50, -50, 50, -50];
end

% This part of the code will be called by other functions to remove
% simulated force outputs when the user changed parameters and has not yet
% ran simulations. 
if(coords(2)>-10)
    fig.UserData.x_total = fig.UserData.x_total*nan;
    fig.UserData.F_total = fig.UserData.F_total*nan;
    fig.UserData.hillF_total = fig.UserData.hillF_total*nan;

    % plotting will be handled by update_time callback with a dummy input
    dummyE.IntersectionPoint = [0.1 0];
    update_time(fig.UserData.ax_len,dummyE)
end

protocol_l = cumsum(fig.UserData.protocol_duration.*fig.UserData.protocol_v);
protocol_t = cumsum(fig.UserData.protocol_duration);
set(fig.UserData.ax_len.Children(2), 'xdata', [0 protocol_t], 'ydata', [0 protocol_l])
end
