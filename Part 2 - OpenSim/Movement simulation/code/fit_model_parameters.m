function[parms] = fit_model_parameters(opti, optparms, w, vmax, RT, SRS_rel, V_rel, parms)

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2'};

for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i}))])
end

% optparms = {'f', 'k11', 'k22', 'k21'};
lb = [1 1 0 1 1];
ub = [2e3 2e3 5 1e3 200];

for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1)'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),')']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),')']);
end

%% design velocity input vector
% this is for testing both force-velocity and history-dependent properties
N = 1000; % number of nodes 
[vts, Fts, toc, idF, idFd] = design_length_input_vector(vmax, RT, V_rel, N);
dt1 = mean(diff(toc));
% 
% plot(toc, vts); hold on
% plot(toc(idF), vts(idF),'o');

%% obtain initial guess
% intial guess is obtained through running a forward simulation with the
% initial parameter values

% first, simulate an isometric contraction
parms.vts = [0 0];
parms.ti = [0 1];

x0 = 1e-3 * ones(5,1);
xp0 = zeros(size(x0));
odeopt = odeset('maxstep', 1e-3);
[~, x0] = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 1], x0, xp0, odeopt);

% next, simulate response to specified velocity input vector
parms.vts = vts;
parms.ti = toc;

sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], x0(end,:), xp0, odeopt);
[~,xdot] = deval(sol, sol.x);

% interpolate solution to time nodes
Q0i     = interp1(sol.x, sol.y(1,:), toc); % zero-order moment
Q1i     = interp1(sol.x, sol.y(2,:), toc); % first-order moment
Q2i     = interp1(sol.x, sol.y(3,:), toc); % second-order moment
Noni    = interp1(sol.x, sol.y(4,:), toc); % thin filament activation
DRXi    = interp1(sol.x, sol.y(5,:), toc); % thick filament activation

dQ0dti  = interp1(sol.x, xdot(1,:), toc); % zero-order moment time derivative
dQ1dti  = interp1(sol.x, xdot(2,:), toc); % first-order moment time derivative
dQ2dti  = interp1(sol.x, xdot(3,:), toc); % second-order moment time derivative
dNondti = interp1(sol.x, xdot(4,:), toc); % thin filament activation time derivative
dDRXdti = interp1(sol.x, xdot(5,:), toc); % thick filament activation time derivative

%% Fit cross-bridge rates using direct collocation
% define opti states (defined as above)
Q0  = opti.variable(1,N);
Q1  = opti.variable(1,N); 
Q2  = opti.variable(1,N); 
Non = opti.variable(1,N);
DRX = opti.variable(1,N);

% define extra variables
p  = opti.variable(1,N); % mean strain of the distribution
q  = opti.variable(1,N); % standard deviation strain of the distribution
  
% (Slack) controls (defined as above)
dQ0dt  = opti.variable(1,N);
dQ1dt  = opti.variable(1,N); 
dQ2dt  = opti.variable(1,N); 
dNondt = opti.variable(1,N); 
dDRXdt = opti.variable(1,N); 

% Inequality constraints
opti.subject_to(Q0 >= 0);
opti.subject_to(Q1 >= -Q0);
opti.subject_to(q >= 0);
opti.subject_to(Non >= 0);
opti.subject_to(Non <= 1);
opti.subject_to(DRX >= 0);
opti.subject_to(DRX <= 1);

% Extra constraints
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% Set initial guess states based on forward simulation results
opti.set_initial(Q0, Q0i);
opti.set_initial(Q1, Q1i);
opti.set_initial(Q2, Q2i);
opti.set_initial(Non, Noni);
opti.set_initial(DRX, DRXi);
opti.set_initial(p, Q1i./Q0i);
opti.set_initial(q, Q2i./Q0i - (Q1i./Q0i).^2);
opti.set_initial(dQ0dt, dQ0dti);
opti.set_initial(dQ1dt, dQ1dti);
opti.set_initial(dQ2dt, dQ2dti);
opti.set_initial(dNondt, dNondti);
opti.set_initial(dDRXdt, dDRXdti);

%% dynamics constraints
F = Q0 + Q1; % cross-bridge force
error = [];
error_thin  = ThinEquilibrium(parms.Ca, Q0, Non, dNondt, parms.kon, parms.koff, koop, parms.Noverlap); % thin filament dynamics     
error_thick = ThickEquilibrium(F, DRX, dDRXdt, J1, J2, JF, parms.Noverlap); % thick filament dynamics
error1      = MuscleEquilibrium(Q0, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22,  Non, vts, DRX); % cross-bridge dynamics
error       = [error; error_thin(:); error_thick(:); error1(:)];
opti.subject_to(error == 0);

%% derivative constraints
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt1/2 + Non(1:N-1) == Non(2:N));
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt1/2 + DRX(1:N-1) == DRX(2:N));
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt1/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt1/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt1/2 + Q2(1:N-1) == Q2(2:N));

%% cost
Frel = F * 2;
Freldot = dQ0dt + dQ1dt;

% cost function
J = 0;
J = J + w(1) * sum((Frel(idF) - Fts(idF)).^2); % force-velocity fitting
J = J + w(2) * sum((SRS_rel - Freldot(idFd(2))/Freldot(idFd(1))).^2); % history-dependence fitting
J = J + w(3) * (sum(dQ0dt(1).^2) + sum(dQ1dt(1).^2) + sum(dQ2dt(1).^2)); % regularization term

opti.minimize(J); 

%% Solve problem
% options for IPOPT
options.ipopt.tol = 1*10^(-6);          
options.ipopt.linear_solver = 'mumps';
% opti.solver('ipopt',options);

% Solve the OCP
p_opts = struct('expand',true);
s_opts = struct('max_iter', 500);
opti.solver('ipopt',p_opts,s_opts);

figure(100); 
opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

try
    sol = opti.solve();  
catch
    sol = opti.debug();
end

close all

%% Plot the result
% obtain the solution
R.Q0    = sol.value(Q0); 
R.Q1    = sol.value(Q1); 
R.Q2    = sol.value(Q2); 
R.dQ0dt = sol.value(dQ0dt); 
R.dQ1dt = sol.value(dQ1dt); 
R.dQ2dt = sol.value(dQ2dt); 
R.F     = sol.value(Frel); 
R.Fdot  = R.dQ0dt + R.dQ1dt;
R.t = 0:dt1:(N-1)*dt1;

figure(1)
subplot(311)
plot(R.t, vts, 'linewidth',1.5); hold on
ylabel('Velocity')

subplot(312)
plot(R.t, R.F, 'linewidth',1.5); hold on
plot(R.t(idF), Fts(idF), '.', 'markersize',10); 
ylabel('Force')

subplot(313);
plot(R.t, R.Fdot, 'linewidth',1.5);
ylabel('Force-rate')

for i = 1:3
    subplot(3,1,i); 
    hold on
    box off
    xlabel('Time (s)')
    xlim([0 max(toc)])
end

%% Test the result
% extract the parameters
for i = 1:length(optparms)
    parms.(optparms{i}) = eval(['sol.value(',optparms{i},');']);
end

% run a forward simulation
osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 max(toc)], x0(end,:), xp0, odeopt);
t = osol.x;
x = osol.y;
F = x(1,:) + x(2,:);
[~,xdot] = deval(osol, t);
Fdot = xdot(1,:) + xdot(2,:);
Fi = interp1(t, F, toc);

figure(1)
color = get(gca,'colororder');

subplot(312)
plot(t, F*2,':','color',color(1,:))
legend('Biophysical','Hill','location','best')
legend boxoff

subplot(313)
plot(t, Fdot*2,':','color',color(1,:))

set(gcf,'units','normalized','position', [.1 .1 .4 .8])

%% Visualize force-velocity
[vs, id] = sort(vts(idF));
Fs = Fi(idF)*2;

Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilda = linspace(0,1.5);
vH = vmax/e2*(sinh((FMvtilda-e4)/e1)-e3); % can be inverted = simpler

figure(2)
subplot(211);
plot(vs, Fs(id), '-','linewidth',2); hold on
plot(vH, FMvtilda, '--'); hold on
xlabel('Velocity (L_0 / s)')
ylabel('Force')
box off
title('Force-velocity')
yline(1,'k--')
xline(0,'k--')

%% Test SRS
% first, simulate isometric
parms.ti = [0 .3];
parms.vts = [0 0];
sol0 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 .3], x0(end,:), xp0, odeopt);
X00 = sol0.y(:,end);

% next, simulate a stretch
parms.ti = [0 .01];
parms.vts = [V_rel V_rel] * vmax;
sol1 = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 .01], X00, xp0, odeopt);
F1 = sol1.y(1,:) + sol1.y(2,:);
SRS1 = (F1(end)-F1(1)) / .01;

% now, simulate the entire stretch-shorten protocol
vt = [0 V_rel -V_rel 0] * vmax;
ts = [.3 .1 .1 10];

Ts = [0 cumsum(ts)];
toc = linspace(0,sum(ts),10000);
vts = zeros(1,N);

for i = 1:(length(Ts)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

parms.ti = toc;
parms.vts = vts;

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], x0(end,:), xp0, odeopt);
t = osol.x;
[~,xdot] = deval(osol, t);
% F = osol.y(1,:) + osol.y(2,:);
% Fdot = xdot(1,:) + xdot(2,:);
% id = t > Ts(2) & t < (Ts(2)+0.01);
% SRS1 = mean(Fdot(id));
% close all
% figure(10)
% plot(t, F); hold on

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

figure(2);
color = get(gca,'colororder');

subplot(212);
semilogx(RT, R.Fdot(idFd(2))/R.Fdot(idFd(1)), '.', 'markersize', 10); hold on
xline(RT, '--','color',color(1,:))
yline(SRS_rel, '--','color',color(1,:))
semilogx(RTs, thix,':','color',color(1,:))
xlabel('Recovery time (s)')
ylabel('Relative short-range stiffness')
box off
title('History dependence')
yline(1,'--','color',color(2,:))

set(gcf,'units','normalized','position', [.5 .1 .4 .8])

