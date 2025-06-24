clear all; close all; clc
cd('C:\Users\u0167448\Documents\GitHub\ISB2025\Part 2 - OpenSim\Movement simulation\input\common')

% Horslen parameters
load('parms.mat','parms')

x0 = 1e-3 * ones(5,1);
xp0 = zeros(size(x0));

parms.a = 1;
parms.vMtilda = 0;   
parms.d = 0;
parms.act = 1;

Cas = 10.^(-1:.1:.5);

colors = parula(length(Cas));

parms.f = 1e3;
parms.k11 = 55.3664;
parms.k12 = 2;
parms.k21 =  451.0874;
parms.k22 =  0.2328;
parms.vMtilda = 0;

parms.vts = [0 0];
parms.ti = [0 1];

parms.Ca = 1;

[t0, x0] = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 1], x0, xp0);
figure(1)
plot(t0,x0)

X0 = x0(end,:);


%%
close all
vt = [0 2 -2 2 -5 2 0];

toc = linspace(0,1,1000);
N = length(toc);
dt1 = mean(diff(toc));

% model constraints
Ns = floor(linspace(0, N, length(vt)+1));
vts = zeros(1,N);

for i = 1:(length(Ns)-1)
    vts(Ns(i)+1:Ns(i+1)) = vt(i);
end

Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];
vmax = 10;

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilda = linspace(0,2);
vH = vmax/e2*(sinh((FMvtilda-e4)/e1)-e3); % 0-1
Fts = interp1(vH, FMvtilda, vts);

parms.vts = vts;
parms.ti = toc;
parms.Ca = 1;

sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 1], X0, xp0);
t = sol.x;
x = sol.y;
[~,xdot] = deval(sol, t);
% 
% clear dX
% for i = 1:length(t)
%     [~, dX(i,:)] = fiber_dynamics_implicit_no_tendon(t(i), x(:,i), zeros(size(x(:,i))), parms);
% end

% figure(1)
% for i = 1:5
%     nexttile
%     
%     plot(t, xdot(i,:)); hold on
%     plot(t, dX(:,i),'--')
% end



F = x(1,:) + x(2,:);

% 
Q0 = x(1,:);
Q1 = x(2,:);
Q2 = x(3,:);
Non = x(4,:);
DRX = x(5,:);

dQ0dt = xdot(1,:);
dQ1dt = xdot(2,:);
dQ2dt = xdot(3,:);

dNondt = xdot(4,:);
dDRXdt = xdot(5,:);

% interpolate back
Q0i = interp1(t, x(1,:), toc);
Q1i = interp1(t, x(2,:), toc);
Q2i = interp1(t, x(3,:), toc);
Noni = interp1(t, x(4,:), toc);
DRXi = interp1(t, x(5,:), toc);

dQ0dti = interp1(t, xdot(1,:), toc);
dQ1dti = interp1(t, xdot(2,:), toc);
dQ2dti = interp1(t, xdot(3,:), toc);
dNondti = interp1(t, xdot(4,:), toc);
dDRXdti = interp1(t, xdot(5,:), toc);

pi = Q1i./Q0i;
qi = Q2i./Q0i - (Q1i./Q0i).^2;

Fi = interp1(t, F, toc);

figure(1)
subplot(131)
plot(toc, vts)

subplot(132)
plot(t, F*2); hold on
plot(toc, Fts, '--')

%
vi = interp1(toc, vts, t);
F0 = F(end);

subplot(133)
plot(vH, FMvtilda); hold on
plot(vi, F*2, '.')

% Check for errors
% close all
% % error_thin = ThinEquilibrium(parms.Ca, Q0, Non, dNondt, parms.kon, parms.koff, parms.koop, parms.Noverlap);
% % error_thick = ThickEquilibrium(Q0, F, DRX, dDRXdt, parms.J1, parms.J2, parms.JF, parms.Noverlap);
% 
% error_thin = ThinEquilibrium(parms.Ca, Q0i, Noni, dNondti, parms.kon, parms.koff, parms.koop, parms.Noverlap);
% error_thick = ThickEquilibrium(Q0i, Fi, DRXi, dDRXdti, parms.J1, parms.J2, parms.JF, parms.Noverlap);
% error_muscle = MuscleEquilibrium_alt_v2(Q0i, pi, qi, dQ0dti, dQ1dti, dQ2dti, parms.f, parms.k11, parms.k12, parms.k21, parms.k22,  Noni, vts, DRXi);
% 
% 
% close all
% 
% figure(3)
% plot(error_thin); hold on
% plot(error_thick,'--')
% plot(error_muscle',':')
% 
% error_thin2 = (dNondti(Ns(1)+1:Ns(end)-1) + dNondti(Ns(1)+2:Ns(end)))*dt1/2 + Noni(Ns(1)+1:Ns(end)-1) - Noni(Ns(1)+2:Ns(end));
% plot(error_thin2,'--')
%% Fit cross-bridge rates using direct collocation
% addpath(genpath('C:\Users\timvd\Documents\casadi-windows-matlabR2016a-v3.5.5'))
addpath(genpath('C:\Users\u0167448\Documents\GitHub\casadi-windows-matlabR2016a-v3.5.5'))
import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

% States
Q0          = opti.variable(1,N);
Q1          = opti.variable(1,N); 
Q2          = opti.variable(1,N); 
Non          = opti.variable(1,N);
DRX          = opti.variable(1,N);

p          = opti.variable(1,N); 
q          = opti.variable(1,N); 
  
% (Slack) controls
dQ0dt        = opti.variable(1,N);
dQ1dt        = opti.variable(1,N); 
dQ2dt        = opti.variable(1,N); 
dNondt       = opti.variable(1,N); 
dDRXdt       = opti.variable(1,N); 

opti.subject_to(Q0 > 0);
opti.subject_to(q > 0);
opti.subject_to(0 < Non <= 1);
opti.subject_to(.1 < DRX <= 1);

% extra constraints
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% opti.subject_to(dQ0dt(1) == 0);
% opti.subject_to(dQ1dt(1) == 0);
% opti.subject_to(dQ2dt(1) == 0);
% opti.subject_to(dNondt(1) == 0);
% opti.subject_to(dDRXdt(1) == 0);
%  opti.subject_to(Q0(1) == Q0i(1));
% opti.subject_to(Q1(1) == Q1i(1));
% opti.subject_to(Q2(1) == Q2i(1));
% opti.subject_to(Non(1) == Nonii(1));
% opti.subject_to(DRX(1) == DRXi(1));

% initial guess states
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

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2'};

for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i}))])
end

optparms = {'f', 'k11', 'k22', 'k21'};
lb = [1 1 0 1];
ub = [2e3 2e3 5 1e3];

for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1)'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),')']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),')']);
end


% dynamic constraints
error = [];
parms.Ca = 1;
parms.Noverlap = 1;

id = (Ns(1)+1):Ns(end);
F = (Q0(id) + Q1(id));

error_thin = ThinEquilibrium(parms.Ca, Q0(id), Non(id), dNondt(id), parms.kon, parms.koff, koop, parms.Noverlap);      
error_thick = ThickEquilibrium(Q0(id), F, DRX(id), dDRXdt(id), J1, J2, JF, parms.Noverlap);
error1 = MuscleEquilibrium(Q0(id), p(id), q(id), dQ0dt(id), dQ1dt(id), dQ2dt(id), f, k11, k12, k21, k22,  Non(id), vts(id), DRX(id));

error = [error; error_thin(:); error_thick(:); error1(:)];

opti.subject_to((dNondt(Ns(1)+1:Ns(end)-1) + dNondt(Ns(1)+2:Ns(end)))*dt1/2 + Non(Ns(1)+1:Ns(end)-1) == Non(Ns(1)+2:Ns(end)));
opti.subject_to((dDRXdt(Ns(1)+1:Ns(end)-1) + dDRXdt(Ns(1)+2:Ns(end)))*dt1/2 + DRX(Ns(1)+1:Ns(end)-1) == DRX(Ns(1)+2:Ns(end)));

opti.subject_to((dQ0dt(Ns(1)+1:Ns(end)-1) + dQ0dt(Ns(1)+2:Ns(end)))*dt1/2 + Q0(Ns(1)+1:Ns(end)-1) == Q0(Ns(1)+2:Ns(end)));
opti.subject_to((dQ1dt(Ns(1)+1:Ns(end)-1) + dQ1dt(Ns(1)+2:Ns(end)))*dt1/2 + Q1(Ns(1)+1:Ns(end)-1) == Q1(Ns(1)+2:Ns(end)));
opti.subject_to((dQ2dt(Ns(1)+1:Ns(end)-1) + dQ2dt(Ns(1)+2:Ns(end)))*dt1/2 + Q2(Ns(1)+1:Ns(end)-1) == Q2(Ns(1)+2:Ns(end)));

opti.subject_to(error == 0);

Frel = F * 2;
Freldot = dQ0dt + dQ1dt;

% cost function
% J = 100 * sumsqr(Frel - Fts);

id1 = Ns(1:end-1)+1;
id2 = Ns(2:end);
id3 = round(id1 + (id2-id1).* .5);

J = 100 * sumsqr(Frel(id3) - Fts(id3));
J = J + 1 * (sumsqr(dQ0dt(1))+sumsqr(dQ1dt(1))+sumsqr(dQ2dt(1))); 

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
% plot(toc, Fts, '--')
opti.callback(@(i) plot(toc, [Fts; opti.debug.value(Frel)]))

try
    sol = opti.solve();  
catch
    sol = opti.debug();
end

%% Get values
R.Q0 = sol.value(Q0); 
R.Q1 = sol.value(Q1); 
R.Q2 = sol.value(Q2); 

R.dQ0dt = sol.value(dQ0dt); 
R.dQ1dt = sol.value(dQ1dt); 
R.dQ2dt = sol.value(dQ2dt); 

Fdot = R.dQ0dt + R.dQ1dt;
R.F = sol.value(Frel); 

% parameters
for i = 1:length(optparms)
    parms.(optparms{i}) = eval(['sol.value(',optparms{i},')']);
end

R.t = 0:dt1:(N-1)*dt1;

close all
figure(1)
subplot(311)
plot(R.t, R.Q0)

subplot(312)
plot(R.t, R.Q1)

subplot(313)
plot(R.t, R.Q2)

figure(2)
subplot(311)
plot(R.t, vts)
ylabel('Velocity')

subplot(312)
plot(R.t, R.F); hold on
plot(R.t, Fts, '--'); 
% plot(R.t(id1), R.F(id1),'o')
% plot(R.t(id2), R.F(id2),'x')
% plot(R.t(id3), R.F(id3),'kx', 'markersize',10,'linewidth',2)
ylabel('Force')
legend('Biophysical','Hill','location','best')
legend boxoff

subplot(313);
plot(R.t, Fdot);
ylabel('Force-rate')

for i = 1:length(Ns)-1
    xline(toc(Ns(i)+1),'k--')
end

cost = 100 * sumsqr(R.F(Ns(2:end)-1) - Fts(Ns(2:end)-1)) + 1 * (sumsqr(R.dQ0dt(1))+sumsqr(R.dQ1dt(1))+sumsqr(R.dQ2dt(1)));

for i = 1:3
    subplot(3,1,i)
    box off
    xlabel('Time (s)')
end

return
%% Test
parms.f = R.f;
parms.Ca = 1;
[ti,xi] = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 1], y0, yp0);

figure(2)
plot(ti, xi(:,1)+xi(:,2))



