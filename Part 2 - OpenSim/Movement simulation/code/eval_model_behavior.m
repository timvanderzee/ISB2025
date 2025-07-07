function[] = eval_model_behavior(vmax, RT, SRS_rel, V_rel, parms, out)

%% Isometric
parms.ti = [0 1];
parms.vts = [0 0];
x0 = 1e-3 * ones(5,1);
xp0 = zeros(5,1);
odeopt = odeset('maxstep', 1e-3);

sol0 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 .3], x0(:), xp0(:), odeopt);
X00 = sol0.y(:,end);
% F0 = sol0.y(1,:) + sol0.y(2,:);

%% Test force-velocity
% formulation and coefficients from De Groote et al. (2016)
d = [ -0.3183   -8.1492   -0.3741    0.8856];
vH = linspace(-1, 1, 20);
FH = d(1) * log((d(2)*vH + d(3)) + sqrt((d(2)*vH + d(3)).^2 + 1)) + d(4);

Fss = nan(1, length(vH));
for i = 1:length(vH)
    parms.ti = [0 .12 / max(abs(vH(i)*vmax), .6)];
    parms.vts = [vH(i) vH(i)] * vmax;
    
    sol1 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), parms.ti, X00, xp0, odeopt);
    F1 = sol1.y(1,:) + sol1.y(2,:);
    
    Fss(i) = F1(end) * 2;
   
end


%% Test SRS
% simulate a stretch at V_rel
parms.ti = [0 .01];
parms.vts = [V_rel V_rel] * vmax;
sol1 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 .01], X00, xp0, odeopt);
F1 = sol1.y(1,:) + sol1.y(2,:);
SRS1 = (F1(end)-F1(1)) / .01;

% now, simulate the entire stretch-shorten protocol
vt = [0 V_rel -V_rel 0] * vmax;
ts = [.3 .1 .1 10];

Ts = [0 cumsum(ts)];
N = 10000;
toc = linspace(0,sum(ts),N);
vts = zeros(1,N);

for i = 1:(length(Ts)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

parms.ti = toc;
parms.vts = vts;

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], x0, xp0, odeopt);
t = osol.x;
[~,xdot] = deval(osol, t);

RTs = logspace(-5,1, 20);
SRS2 = nan(length(RTs), 1);
for j = 1:length(RTs) 
    
    X0 = interp1(t(:), osol.y', Ts(4) + RTs(j));
    Xp0 = interp1(t(:), xdot', Ts(4) + RTs(j));
    
    parms.ti = [0 .01];
    parms.vts = [V_rel V_rel] * vmax;

    nsol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 .01], X0, Xp0(:), odeopt);

    nt = nsol.x;
    [~,nxdot] = deval(nsol, nt);
    Fdot2 = nxdot(1,:) + nxdot(2,:);
    
    F2 = nsol.y(1,:) + nsol.y(2,:);
%     hold on
%     plot(nsol.x + Ts(4) + RTs(j), F2);
    
%     id2 = nt < .01;
    SRS2(j) = (F2(end) - F2(1)) / .01;

end

% thixotropy: SRS reduction with pre-movement
thix = SRS2 ./ SRS1;

%% plotting
figure(2)
color = get(gca,'colororder');

subplot(121);
plot(out.v, out.Ft, 'o', 'color', color(3,:)); hold on
plot(vH * vmax, FH, '--', 'color', color(2,:)); hold on
plot(out.v, out.F,'o','color',color(1,:),'markerfacecolor', color(1,:))
plot(vH * vmax, Fss, '-','linewidth',1, 'color', color(1,:)); hold on

legend('Target' ,'Hill','Biophysical', 'location','best')
legend boxoff

xlabel('Velocity (L_0 / s)')
ylabel('Force')
box off
title('Force-velocity')


subplot(122);
plot(RT, SRS_rel, 'o', 'color', color(3,:)); hold on
yline(1,'--','color',color(2,:))

plot(RT, out.SRS, 'o', 'color', color(1,:),'markerfacecolor', color(1,:))
semilogx(RTs, thix,'-','color',color(1,:),'linewidth',1); hold on
xline(RT, '-.','color',color(3,:)); hold on

yline(SRS_rel, '-.','color',color(3,:)); hold on


set(gca,'XScale','log')
legend('Target', 'Hill','Biophysical', 'location','best')
legend boxoff

xlabel('Recovery time (s)')
ylabel('Relative short-range stiffness')
box off
title('History dependence')

ylim([0 1.2])

set(gcf,'units','normalized','position', [.1 .3 .8 .4])