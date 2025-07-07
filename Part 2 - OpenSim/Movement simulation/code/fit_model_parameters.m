function[parms, out] = fit_model_parameters(opti, optparms, w, vmax, RT, SRS_rel, V_rel, parms)

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2'};

for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i})),';'])
end

lb = [1 1 0 1 1];
ub = [2e3 2e3 5 1e3 200];

for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1);'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),');']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),');']);
end

%% design velocity input vector
% this is for testing both force-velocity and history-dependent properties
N = 1000; % number of nodes 
[vts, Fts, toc, idF, idFd] = design_length_input_vector(vmax, RT, V_rel, N);
dt = mean(diff(toc));

%% obtain initial guess
% intial guess is obtained through running a forward simulation with the
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

k = 20;
Nonc = log(1+exp(Non*k))/k;
DRXc = log(1+exp(DRX*k))/k;
Q0c = log(1+exp(Q0*k))/k;
Fc = log(1+exp(F*k))/k;

error = [];
error_thin  = ThinEquilibrium(parms.Ca, Q0c, Nonc, dNondt, parms.kon, parms.koff, koop, parms.Noverlap); % thin filament dynamics     
error_thick = ThickEquilibrium(Fc, DRXc, dDRXdt, J1, J2, JF, parms.Noverlap); % thick filament dynamics
error1      = MuscleEquilibrium(Q0c, p, q, dQ0dt, dQ1dt, dQ2dt, f, parms.w, k11, k12, k21, k22,  Nonc, vts, DRXc); % cross-bridge dynamics
error       = [error; error_thin(:); error_thick(:); error1(:)];
opti.subject_to(error == 0);

%% derivative constraints
opti.subject_to((dNondt(1:N-1) + dNondt(2:N))*dt/2 + Non(1:N-1) == Non(2:N));
opti.subject_to((dDRXdt(1:N-1) + dDRXdt(2:N))*dt/2 + DRX(1:N-1) == DRX(2:N));
opti.subject_to((dQ0dt(1:N-1) + dQ0dt(2:N))*dt/2 + Q0(1:N-1) == Q0(2:N));
opti.subject_to((dQ1dt(1:N-1) + dQ1dt(2:N))*dt/2 + Q1(1:N-1) == Q1(2:N));
opti.subject_to((dQ2dt(1:N-1) + dQ2dt(2:N))*dt/2 + Q2(1:N-1) == Q2(2:N));

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
R.t = 0:dt:(N-1)*dt;

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

%% output
out.v = vts(idF);
out.F = R.F(idF);
out.SRS = R.Fdot(idFd(2))/R.Fdot(idFd(1));
out.Ft = Fts(idF);

end
