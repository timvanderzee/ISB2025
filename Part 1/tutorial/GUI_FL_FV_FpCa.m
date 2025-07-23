% clear; clc; close all

addpath(genpath([pwd,'/../..']))

warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

parms.forcible_detachment = 0;
parms.kse = 0;
parms.kpe = 0;
parms.no_tendon = 1;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;
parms.nbins = 500;
parms.xi0 = linspace(-15,15,parms.nbins);
parms.xi = parms.xi0;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0.0909;

fig = figure;

fig.UserData.parms = parms;
parm_range = [100, 1, 30, 3, 30, 3];

ax = subplot(4,3,4);

fig.UserData.ax_parm = ax;
plot(ax,[0:1, 3:6],...
    [parms.f,parms.w,parms.k11,parms.k12,parms.k21,parms.k22]./parm_range,...
    '.', 'markerSize', 20, 'PickableParts', 'none');
xlim([-1 7])
xticks([0:1, 3:6])
xticklabels({'f','w','k_{11}','k_{12}','k_{21}','k_{22}'});
text(ax, 0, 1.3, 'attachment', 'fontSize', 10, 'HorizontalAlignment','center');
text(ax, 4.5, 1.3, 'detachment', 'fontSize', 10, 'HorizontalAlignment','center');


fig.UserData.ax_ffunc = subplot(4,3,7);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('attachment rate (f)')

fig.UserData.ax_gfunc = subplot(4,3,10);
plot(parms.xi, nan*parms.xi);
xlabel('\Delta x')
ylabel('detachment rate (g)')

set(fig.UserData.ax_parm,'ButtonDownFcn', ...
    @(s,e)update_parms(s,e,parm_range),...
    'HitTest','on');
dummyE.IntersectionPoint = [-2,0];
update_parms(fig.UserData.ax_parm,dummyE,parm_range);


ax = subplot(4,3,1);
fig.UserData.ax_run = ax;
set(fig.UserData.ax_run,'ButtonDownFcn', ...
    @(s,e)run_sim_callback(s,e),...
    'HitTest','on');
fig.UserData.statusHandle = text(0.5, 0.5, 'click here to RUN',...
    'HorizontalAlignment','center','PickableParts','none');
xlim([0 1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.Box = 'off';



%%
function [] = run_sim_callback(src,~)
fig = src.Parent;
set(fig.UserData.statusHandle, 'string', 'running - please wait');
pause(0.01);
sim_flfvfpca(fig,fig.UserData.parms);
set(fig.UserData.statusHandle, 'string', 'click here to RUN');
end

function [] = sim_flfvfpca(fig,parms)

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
save_hill_properties = 0;
hill_properties = struct();

ref_pCa = 4.5; % << pCa during F-L, F-v simulation 

for cond_itr = 1:3 % << pick which protocol you want to run 
    % 1: F-L, 2: F-v, 3: F-pCa
   
    if(cond_itr==1) % F-L
        l0 = [-0.5:0.1:0.5] *half_s_len_norm;
    elseif(cond_itr==2) % F-v
        l0 = [-1:0.2:1] *half_s_len_norm / 20;
    else % F-pCa
        pCa_list = [4.5,5,5.5:0.1:6.5,7:0.5:9];
        l0 = zeros(size(pCa_list));
    end

    ax = subplot(6,3,cond_itr*6-1-3, 'colorOrder', winter(length(l0))); 
    cla(ax)
    hold on
    ax = subplot(6,3,cond_itr*6-1, 'colorOrder', winter(length(l0))); 
    cla(ax)
    hold on

    F0 = nan(size(l0));

    for k=1:length(l0)
        if(cond_itr==3)
            pCa = pCa_list(k);
        else
            pCa = ref_pCa;
        end
        Ca = 10^(-pCa+6);
        model = @fiber_dynamics;

        parms.xss(end-1) = l0(k);
        parms.xss(end) = l0(k);

        parms.lce0 = l0(k);

         if(cond_itr==1||cond_itr==3) % F-L & F-pCa
             us = [0,0];
             Ts = [2,1];
         else % F-v
             us = [0,-l0(k)*20,0];
             Ts = [3,0.05,3];
         end

        x0 = parms.xss;
        t = 0;
        x = x0;

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

        F = nan(1,height(x));
        for i = 1:height(x)
            [~,F(i)] = model(t(i), x(i,:)', parms, Ca);
        end

        subplot(6,3,cond_itr*6-1-3)
        plot(t, x(:,end)/half_s_len_norm)
        hold on

        subplot(6,3,cond_itr*6-1)
        plot(t,F)
        hold on
        pause(0.1)
        F0(k) = F(temp_idx);
    end

    subplot(3,3,cond_itr*3)
    if(cond_itr==1) % F-L
        plot(l0/half_s_len_norm, F0);
        xlabel('\Delta length (l_{opt})')
        ylabel('F (F_0)')
        title('F-l')
        hill_properties.FL_spline = spline(l0/half_s_len_norm, F0);
    elseif(cond_itr==2) % F-v
        plot(l0/half_s_len_norm*10, F0);
        xlabel('velocity (l_{opt}/s)')
        ylabel('F (F_0)')
        title('F-v')
        hill_properties.FV_spline = spline(l0/half_s_len_norm*10, F0);
    else % F-pCa
        plot(pCa_list, F0);
        xlabel('pCa')
        ylabel('F (F_0)')
        title('F-pCa')
        hill_properties.FPca_spline = spline(pCa_list, F0);
    end
    ax1 = subplot(6,3,cond_itr*6-1-3);
    ylabel('\Delta length (l_{opt})')
    ax2 = subplot(6,3,cond_itr*6-1);
    ylabel('F (F_0)')
    xlabel('time (s)')
    linkaxes([ax1, ax2],'x');
    if(cond_itr == 2)
        xlim([3, 3.1])
        plot(ones(size(F0))*3.05, F0, '.', 'markerSize', 10);
    end

end

if(save_hill_properties)
    save(sprintf('hill_properties_pCa_%d', ref_pCa*10), "hill_properties");
end



end

%% 
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
    end
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
end