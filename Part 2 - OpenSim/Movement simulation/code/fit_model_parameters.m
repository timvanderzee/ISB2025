function[parms] = fit_model_parameters(opti)

% cd('C:\Users\u0167448\Documents\GitHub\ISB2025\Part 2 - OpenSim\Movement simulation\input\common')

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


%% Design velocity vector
vmax = 5;

vt = [0 -2 2 -4 4 -6 0 5 -5 0 5]/10 * vmax;

% for velocity part, we just need to evaluate until evaluation time
idv = [2 3 4 5 6];

ts = .2 * ones(size(vt));
for i = 1:length(idv)
    ts(i+1) = .12 / abs(vt(idv(i)));
end

idS = 7:length(vt);
ts(idS) = [.3 .1 .1 .1 .1];
Ts = [0 cumsum(ts)];

toc = linspace(0,sum(ts),1000);
N = length(toc);
dt1 = mean(diff(toc));

% model constraints
Ns = floor(linspace(0, N, length(vt)+1));

vts = zeros(1,N);
for i = 1:(length(Ns)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

% close all
% figure(1)
% plot(toc, vts); hold on

%
idF = nan(1, length(idv)+1);

idF(1) = find(toc < Ts(2), 1, 'last');
for i = 1:length(idv)
    idF(i+1) = find(toc < Ts(i+1), 1, 'last');
end

% plot(toc(idF),vts(idF),'o')

idFd = [find(toc > Ts(end-4),1); find(toc > Ts(end-1),1)];
% plot(toc(idFd), vts(idFd),'x')


%% Get initial state
close all
Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];

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

sol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], X0, xp0);
t = sol.x;
x = sol.y;
[~,xdot] = deval(sol, t);

F = x(1,:) + x(2,:);

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

% figure(1)
% subplot(131)
% plot(toc, vts)
% 
% subplot(132)
% plot(t, F*2); hold on
% plot(toc, Fts, '--')
% 
% %
% vi = interp1(toc, vts, t);
% F0 = F(end);
% 
% subplot(133)
% plot(vH, FMvtilda); hold on
% plot(vi, F*2, '.')


%% Fit cross-bridge rates using direct collocation
% addpath(genpath('C:\Users\timvd\Documents\casadi-windows-matlabR2016a-v3.5.5'))

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

opti.subject_to(Q0 >= 0);
opti.subject_to(Q1 >= -Q0);
opti.subject_to(q >= 0);
opti.subject_to(Non >= 0);
opti.subject_to(Non <= 1);
opti.subject_to(DRX >= 0);
opti.subject_to(DRX <= 1);

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

% make sure things can't be negative
k = 20;
Nonc = log(1+exp(Non*k))/k;
DRXc = log(1+exp(DRX*k))/k;
Q0c = log(1+exp(Q0*k))/k;
Fc = log(1+exp(F*k))/k;

error_thin = ThinEquilibrium(parms.Ca, Q0c(id), Nonc(id), dNondt(id), parms.kon, parms.koff, koop, parms.Noverlap);      
error_thick = ThickEquilibrium(Q0c(id), Fc, DRXc(id), dDRXdt(id), J1, J2, JF, parms.Noverlap);
error1 = MuscleEquilibrium(Q0c(id), p(id), q(id), dQ0dt(id), dQ1dt(id), dQ2dt(id), f, k11, k12, k21, k22,  Nonc(id), vts(id), DRXc(id));

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

% id1 = Ns(1:end-1)+1;
% id2 = Ns(2:end);
% id3 = round(id1 + (id2-id1).* .5);


J = 100 * sum((Frel(idF) - Fts(idF)).^2);
J = J + 100 * sum((0.7 - Freldot(idFd(2))/Freldot(idFd(1))).^2);
J = J + 1 * (sum(dQ0dt(1).^2) + sum(dQ1dt(1).^2) + sum(dQ2dt(1).^2)); 

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

close all

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

% close all
% figure(1)
% subplot(311)
% plot(R.t, R.Q0); hold on
% 
% subplot(312)
% plot(R.t, R.Q1); hold on
% 
% subplot(313)
% plot(R.t, R.Q2); hold on

figure(1)
subplot(311)
plot(R.t, vts); hold on
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

% cost = 100 * sumsqr(R.F(Ns(2:end)-1) - Fts(Ns(2:end)-1)) + 1 * (sumsqr(R.dQ0dt(1))+sumsqr(R.dQ1dt(1))+sumsqr(R.dQ2dt(1)));

for i = 1:3
    subplot(3,1,i); hold on
    box off
    xlabel('Time (s)')
end



%% Test protocol
osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], X0, xp0);
t = osol.x;
x = osol.y;
F = x(1,:) + x(2,:);
[~,xdot] = deval(osol, t);
Fdot = xdot(1,:) + xdot(2,:);

figure(1)
subplot(312)
plot(t, F*2,':')

subplot(313)
plot(t, Fdot*2,':')

%% Test force-velocity
vmax = 5;

vt = [0 -2 2 -3 3 -4 4 -5 5 -6 6 -7 7 -8 8 -9 9]/10 * vmax;

% for velocity part, we just need to evaluate until evaluation time
idv = 2:length(vt);

ts = .2 * ones(size(vt));
for i = 1:length(idv)
    ts(i+1) = .12 / abs(vt(idv(i)));
end

Ts = [0 cumsum(ts)];
toc = linspace(0,sum(ts),1000);
N = length(toc);

% model constraints
Ns = floor(linspace(0, N, length(vt)+1));

vts = zeros(1,N);
for i = 1:(length(Ts)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

parms.ti = toc;
parms.vts = vts;

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], X0, xp0);
t = osol.x;
x = osol.y;
F = x(1,:) + x(2,:);
[~,xdot] = deval(osol, t);
Fdot = xdot(1,:) + xdot(2,:);

idF = nan(1, length(vt));
for i = 1:length(vt)
    idF(i) = find(t < (Ts(i+1)-.001), 1, 'last');
end

Fss = F(idF) * 2;

% close all
% figure(3)
% 
% subplot(121);
% plot(t, F*2); hold on
% plot(t(idF), F(idF)*2,'o')

v = interp1(toc, vts, t);

[~, id] = sort(vt);

figure(2)
subplot(121);
% subplot(122)
plot(vt(id), Fss(id), '-','linewidth',2); hold on
% plot(v, F, '.');

% Test SRS
vmax = 5;
% close all

ISIs = logspace(-5,0, 10);
thix = nan(size(ISIs));

for j = 1:length(ISIs)

vt = [0 .5 -.5 0 .5] * vmax;
ts = [.2 .2 .2 ISIs(j) .2];


Ts = [0 cumsum(ts)];
toc = linspace(0,sum(ts),1000);
N = length(toc);

% model constraints
Ns = floor(linspace(0, N, length(vt)+1));

vts = zeros(1,N);

for i = 1:(length(Ts)-1)
    id = (toc > Ts(i)) & (toc <= Ts(i+1));
    vts(id) = vt(i);
end

parms.ti = toc;
parms.vts = vts;

osol = ode15i(@(t,y,yp) fiber_dynamics_implicit_no_tendon(t,y,yp, parms), [0 Ts(end)], X0, xp0);
t = osol.x;
x = osol.y;
F = x(1,:) + x(2,:);
[~,xdot] = deval(osol, t);
Fdot = xdot(1,:) + xdot(2,:);


idFd1 = find(t > Ts(end-4) & t < (Ts(end-4) + .01));
idFd2 = find(t > Ts(end-1) & t < (Ts(end-1) + .01));
% 
% close all
% figure(1)
% subplot(311);
% plot(toc, vts)
% 
% subplot(312);
% plot(t, F*2); hold on
% plot(t(idFd2), F(idFd2)*2,'.')
% plot(t(idFd1), F(idFd1)*2,'.')
% 
% subplot(313);
% plot(t, Fdot*2); hold on
% plot(t(idFd2), Fdot(idFd2)*2,'.')
% plot(t(idFd1), Fdot(idFd1)*2,'.')

thix(j) = mean(Fdot(idFd2)) / mean(Fdot(idFd1));

end


figure(2);

subplot(122);
semilogx(ISIs, thix,'linewidth',2)

%% make nice
figure(2)

subplot(121)
xlabel('Velocity (L_0 / s)')
ylabel('Force')
box off
title('Force-velocity')
yline(1,'k--')
xline(0,'k--')

subplot(122)
xlabel('Recovery time (s)')
ylabel('Relative short-range stiffness')
box off
title('History dependence')
yline(1,'k--')

