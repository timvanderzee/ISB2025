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

parms.forcible_detachment = 0;
parms.kse = 0;
parms.kpe = 0;
parms.no_tendon = 1;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.xi0 = linspace(-10,10,nbins);
parms.xi = parms.xi0;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0;%0.0909;
parms.pCa_ran = 4.5;

model = @fiber_dynamics;
parms0 = parms;

%%
fig = figure;

fig.UserData.parms = parms;
fig.UserData.protocol_duration = [0.05,0.1,0.1,0.1,0.1]; % should be integer multiplication of dt
fig.UserData.protocol_v = [0, -10, 10, 10, -10];
fig.UserData.protocol_pCa = [9 [1 1 1 1]*parms.pCa_ran];

load('hill_properties_pCa_45.mat');
fig.UserData.hill_properties = hill_properties;

sim_XB(fig, model,parms);
parm_range = [100, 1, 30, 3, 30, 3, 4.5];

ax = subplot(5,10,10);
fig.UserData.ax_run = ax;
set(fig.UserData.ax_run,'ButtonDownFcn', ...
    @(s,e)run_sim_callback(s,e,model),...
    'HitTest','on');
fig.UserData.statusHandle = text(0.5, 0.5, 'click here to RUN',...
    'HorizontalAlignment','center','PickableParts','none');
xlim([0 1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.Box = 'off';


ax = subplot(5,10,9);
fig.UserData.ax_save = ax;
set(fig.UserData.ax_save,'ButtonDownFcn', ...
    @(s,e)load_callback(s,e),...
    'HitTest','on');
fig.UserData.loadTextHandle = text(0.5, 0.5, 'load Hill',...
    'HorizontalAlignment','center','PickableParts','none');
xlim([0 1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.Box = 'off';


fig.UserData.ax_protocol = subplot(5,10,6:8);
plot(1,0,'.','markerSize',20);
title('protocol')
xlim([0 3]);
set(fig.UserData.ax_protocol,'ButtonDownFcn', ...
    @(s,e)len_protocol_callback(s,e),...
    'HitTest','on');
xticks(0:3)
xticklabels({'shorten','lengthen', ...
    sprintf('shorten cycle'),sprintf('lengthen cycle')})
fig.UserData.ax_protocol.YAxis.Visible = 'off';

ax = subplot(5,2,1);
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

fig.UserData.ax_ffunc = subplot(5,2,3);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('attachment rate (f)')

fig.UserData.ax_gfunc = subplot(5,2,5);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('detachment rate (g)')

fig.UserData.ax_pCa = subplot(5,2,6);
plot(fig.UserData.t_protocol, fig.UserData.pCa_protocol, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [4 9.5], 'PickableParts', 'none');
xlabel('time (s)')
ylabel('pCa')

fig.UserData.ax_len = subplot(5,2,4);
plot(nan, nan, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [-10 10], 'PickableParts', 'none');
xlabel('time (s)')
ylabel('\Delta length (ps)')

fig.UserData.ax_F = subplot(5,2,[8,10]);
hxb = plot(fig.UserData.t_total, fig.UserData.F_total, 'PickableParts', 'none');
hold on
plot([0.1 0.1], [0 max(fig.UserData.F_total)*3], 'PickableParts', 'none');
hhill = plot(fig.UserData.t_total, fig.UserData.hillF_total, 'PickableParts', 'none');
legend([hxb hhill], 'XB', 'Hill', 'Location', 'northeast')
xlabel('time (s)')
ylabel('force (F_0)')

fig.UserData.ax_dist = subplot(5,2,[7,9]);
plot_k = 100;
ax_dist = plot(fig.UserData.parms.xi0, fig.UserData.x_total(plot_k, 1:fig.UserData.parms.nbins));
xlabel('\Delta x (ps)')
ylabel('bound XB fraction')

set(fig.UserData.ax_len,'ButtonDownFcn', ...
    @(s,e)update_time(s,e),...
    'HitTest','on')
set(fig.UserData.ax_F,'ButtonDownFcn', ...
    @(s,e)update_time(s,e),...
    'HitTest','on')
set(fig.UserData.ax_pCa,'ButtonDownFcn', ...
    @(s,e)update_time(s,e),...
    'HitTest','on')

set(fig.UserData.ax_parm,'ButtonDownFcn', ...
    @(s,e)update_parms(s,e,parm_range),...
    'HitTest','on');

dummyE.IntersectionPoint = [1 -100];
len_protocol_callback(fig.UserData.ax_len, dummyE);

dummyE.IntersectionPoint = [0.1 0];
update_time(fig.UserData.ax_len, dummyE)

dummyE.IntersectionPoint = [-2,0];
update_parms(fig.UserData.ax_parm,dummyE,parm_range);
linkaxes([fig.UserData.ax_len, fig.UserData.ax_F, fig.UserData.ax_pCa],'x');

run_sim_callback(fig.UserData.ax_run, [], model);

%% 
function [] = load_callback(src,~)
fig = src.Parent;
uiopen('load');
if(exist('hill_properties'))
    fig.UserData.hill_properties = hill_properties;
    set(fig.UserData.ax_F.Children(1), 'ydata', fig.UserData.hillF_total*nan);
end
end

function [] = update_parms(src,eventData, parm_range)
fig = src.Parent;

coords = eventData.IntersectionPoint;
x = coords(1);
y = max(0,min(1,coords(2)));

y_s = get(fig.UserData.ax_parm.Children(end),'ydata');

if(x<-1)
else
    if(x<0.5)
        fig.UserData.parms.f = y * parm_range(1);
        y_s(1) = y;
    elseif(x<1.5)
        fig.UserData.parms.w = y * parm_range(2);
        y_s(2) = y;
    elseif(x<3.5)
        fig.UserData.parms.k11 = y * parm_range(3);
        y_s(3) = y;
    elseif(x<4.5)
        fig.UserData.parms.k12 = y * parm_range(4);
        y_s(4) = y;
    elseif(x<5.5)
        fig.UserData.parms.k21 = y * parm_range(5);
        y_s(5) = y;
    elseif(x<6.5)
        fig.UserData.parms.k22 = y * parm_range(6);
        y_s(6) = y;
    elseif(x<8.5)
        fig.UserData.parms.pCa_ran = y * parm_range(7) + 4.5;
        y_s(7) = y;
    end
    fig.UserData.x_total = fig.UserData.x_total*nan;
    fig.UserData.F_total = fig.UserData.F_total*nan;
    fig.UserData.hillF_total = fig.UserData.hillF_total*nan;

    dummyE.IntersectionPoint = [0.1 0];
    update_time(fig.UserData.ax_len,dummyE)
end
set(fig.UserData.ax_parm.Children(end),'ydata',y_s)
set(fig.UserData.ax_parm, 'ylim',[0 1]);

parms = fig.UserData.parms;

set(fig.UserData.ax_ffunc.Children(1), 'yData', ...
    parms.f_func(parms.xi, parms.f, parms.w));
set(fig.UserData.ax_ffunc, 'ylim', [0 300], 'xlim', [min(parms.xi0) max(parms.xi0)]);

set(fig.UserData.ax_gfunc.Children(1), 'yData', ...
    parms.g_func(parms.xi, parms.k11, -parms.k12)+...
    parms.g_func(parms.xi, parms.k21, parms.k22));
set(fig.UserData.ax_gfunc, 'ylim', [0 6e5], 'xlim', [min(parms.xi0) max(parms.xi0)]);

pCa_temp = get(fig.UserData.ax_pCa.Children(2), 'yData');
pCa_temp(3:end) = fig.UserData.parms.pCa_ran;
set(fig.UserData.ax_pCa.Children(2), 'yData', pCa_temp);
fig.UserData.protocol_pCa = pCa_temp;
end

%%
function [] = update_time(src,eventData)
fig = src.Parent;

t_total = fig.UserData.t_total;
x_total = fig.UserData.x_total;
F_total = fig.UserData.F_total;

coords = eventData.IntersectionPoint;
t = min(max(0,coords(1)),t_total(end));
idx = find(t_total>=t,1,'first');

set(fig.UserData.ax_dist.Children(1), 'ydata', x_total(idx,2:end-3), ...
    'xdata', fig.UserData.parms.xi0 + x_total(idx,end));

if(~sum(isnan(x_total)))
    set(fig.Children(1), 'ylim', [0 max(max(x_total(:,2:end-3)))], 'xlim', [-10 10]);
end

set(fig.UserData.ax_F.Children(1), 'ydata', fig.UserData.hillF_total);
set(fig.UserData.ax_F.Children(2), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_F.Children(3), 'ydata', F_total);

if(~sum(isnan(fig.UserData.F_total)))
    set(fig.UserData.ax_F, 'ylim', ...
        [0 max([fig.UserData.F_total; fig.UserData.hillF_total])*1.1]);
end
set(fig.UserData.ax_len.Children(1), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_len, 'ylim', [-1 1]*6);

set(fig.UserData.ax_pCa.Children(1), 'xdata', [1 1]*t_total(idx));
set(fig.UserData.ax_pCa, 'ylim', [4 9.5]);

end

%%
function [] = run_sim_callback(src,~,model)
fig = src.Parent;
set(fig.UserData.statusHandle, 'string', 'running - please wait');
pause(0.01);
sim_XB(fig,model,fig.UserData.parms);
dummyE.IntersectionPoint = [10 0];
update_time(fig.UserData.ax_len,dummyE)

set(fig.UserData.statusHandle, 'string', 'RUN');
end

%%
function [] = len_protocol_callback(src,eventData)
fig = src.Parent;

coords = eventData.IntersectionPoint;
x = coords(1);

fig.UserData.protocol_duration = [0.05, 0.1, 0.1, 0.1,0.1]; % should be integer multiplication of dt
fig.UserData.protocol_v = [0, 0, 50, 0, 0];
fig.UserData.protocol_pCa = [9 [1 1 1 1]*fig.UserData.parms.pCa_ran];

if(x<0.5)
    set(src.Children(1),'XData',0);
    fig.UserData.protocol_v = [0, 0, -50, 0, 0];
elseif(x<1.5) 
    set(src.Children(1),'XData',1);
    fig.UserData.protocol_v = [0, 0, 50, 0, 0];
elseif(x<2.5)
    set(src.Children(1),'XData',2);
    fig.UserData.protocol_v = [0, -50, 50, -50, 50];
else
    set(src.Children(1),'XData',3);
    fig.UserData.protocol_v = [0, 50, -50, 50, -50];
end

if(coords(2)>-10)
    fig.UserData.x_total = fig.UserData.x_total*nan;
    fig.UserData.F_total = fig.UserData.F_total*nan;
    fig.UserData.hillF_total = fig.UserData.hillF_total*nan;

    dummyE.IntersectionPoint = [0.1 0];
    update_time(fig.UserData.ax_len,dummyE)
end

protocol_l = cumsum(fig.UserData.protocol_duration.*fig.UserData.protocol_v);
protocol_t = cumsum(fig.UserData.protocol_duration);
set(fig.UserData.ax_len.Children(2), 'xdata', [0 protocol_t], 'ydata', [0 protocol_l])
end

%%
function [] = sim_XB(fig,model,parms)

odeopt = odeset('maxstep',1e-2);

x0 = parms.xss;
dt = 0.001;

protocol_duration = fig.UserData.protocol_duration;
protocol_v = fig.UserData.protocol_v;
protocol_pCa = fig.UserData.protocol_pCa;

t_total = nan(round(sum(protocol_duration)/dt)+length(protocol_duration),1);
x_total = nan(length(t_total),length(x0));
v_total = nan(size(t_total));
F_total = nan(size(t_total));
pCa_total = nan(size(t_total));

t_protocol = nan(length(protocol_pCa)*2, 1);
pCa_protocol = nan(length(protocol_pCa)*2, 1);

curr_idx = 1;
curr_time = 0;
for p_itr = 1:length(protocol_v)
    t_protocol((-1:0)+p_itr*2) = [curr_time, curr_time+protocol_duration(p_itr)];
    pCa_protocol((-1:0)+p_itr*2) = protocol_pCa(p_itr);

    parms.vmtc = protocol_v(p_itr);
    Ca = 10.^(-protocol_pCa(p_itr)+6);

    [t_sim,x_sim] = ode15s(model, ...
        curr_time:dt:curr_time+protocol_duration(p_itr), x0, odeopt, ...
        parms, Ca);
    F_sim = nan(size(t_sim));
    for k_itr = 1:length(t_sim)
        [~,F_sim(k_itr)] = model(t_sim(k_itr), x_sim(k_itr,:)', parms, Ca);
    end
    t_total(curr_idx:curr_idx+length(t_sim)-1) = t_sim;
    pCa_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_pCa(p_itr);
    v_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_v(p_itr);

    F_total(curr_idx:curr_idx+length(t_sim)-1) = F_sim;
    x_total(curr_idx:curr_idx+length(t_sim)-1,:) = x_sim;
    x0 = x_sim(end,:);
    curr_idx = curr_idx + length(t_sim);
    curr_time = t_sim(end);
end

fig.UserData.t_total = t_total;
fig.UserData.x_total = x_total;
fig.UserData.F_total = F_total;
fig.UserData.pCa_total = pCa_total;
fig.UserData.v_total = v_total;

fig.UserData.t_protocol = t_protocol;
fig.UserData.pCa_protocol = pCa_protocol;

sim_Hill(fig,fig.UserData.hill_properties,...
    parms);
end

function sim_Hill(fig,hill_p,parms)
% Hill model, derived for no tendon, no elastic elements 
half_s_len_norm = parms.s/2/parms.h;

l_total = fig.UserData.x_total(:,end); 
t_total = fig.UserData.t_total;
pCa_total = fig.UserData.pCa_total;
v_total = fig.UserData.v_total;

fig.UserData.hillF_total = ppval(hill_p.FL_spline, l_total/half_s_len_norm).*...
    ppval(hill_p.FV_spline, -v_total/half_s_len_norm).*...
    ppval(hill_p.FPca_spline, pCa_total);
end