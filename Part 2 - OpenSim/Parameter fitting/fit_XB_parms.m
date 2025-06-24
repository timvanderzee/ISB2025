close all; clc; clear all

% Horslen parameters
load('parms.mat')

y0 = 1e-3 * ones(5,1);
yp0 = zeros(size(y0));

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
%         parms.JF = 500;
%         parms.koop = 80;


for i = 1:length(Cas)

    parms.Ca = Cas(i);

    [t0,x0] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 1], y0, yp0);

    F0s(i) = x0(end,1) + x0(end,2);
    x0s(i,:) = x0(end,:);

    subplot(221)
    plot(t0,x0(:,1) + x0(:,2), 'color', colors(i,:)); hold on

    subplot(222)
    plot(Cas(i), F0s(i), 'o', 'color', colors(i,:),'markerfacecolor', colors(i,:)); hold on
end


subplot(222)
semilogx(Cas, F0s,'-'); hold on    

parms.Ca = 10^0;
xline(parms.Ca,'k--')

vs = linspace(10,-10,length(Cas));

[min, id] = min(abs(Cas - parms.Ca));
F0 = F0s(id);
x0 = x0s(id, :);


for i = 1:length(vs)

    parms.vMtilda = vs(i);
     
    [ti,xi] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 .002], x0(end,:), yp0);
    
    Fi = xi(:,1) + xi(:,2);
    
    subplot(223)
    plot(ti, Fi, 'color', colors(i,:)); hold on

    Fss(i,1) = Fi(end);
    
    subplot(224)
    plot(vs(i), Fss(i)/F0, '.', 'color', colors(i,:)); hold on
   
end

Fvparam = [ -0.3183   -8.1492   -0.3741    0.8856];
vmax = 10;

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilda = linspace(0,2);
vH = vmax/e2*(sinh((FMvtilda-e4)/e1)-e3); % 0-1


figure(1)
subplot(224)
plot(vs(:), (Fss(:,1) + parms.d*vs(:))/F0); hold on
% plot(vs(:), (Fss(:,2) + parms.d*vs(:))/F0,'--'); hold on
plot(vH, FMvtilda,':')
xline(0,'k--')
yline(0,'k--')

%% Get initial guess
parms.vMtilda = 0;
% parms.a = 1;

parms.Ca = .1;
[ti,xi] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 1], y0, yp0);

%% Fit cross-bridge rates using direct collocation
% addpath(genpath('C:\Users\timvd\Documents\casadi-windows-matlabR2016a-v3.5.5'))
addpath(genpath('C:\Users\u0167448\Documents\GitHub\casadi-windows-matlabR2016a-v3.5.5'))
import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

toc = linspace(0,.1,1000);
N = length(toc);

% States
Q0          = opti.variable(1,N);
Q1          = opti.variable(1,N); 
Q2          = opti.variable(1,N); 

p          = opti.variable(1,N); 
q          = opti.variable(1,N); 
  
% (Slack) controls
dQ0dt          = opti.variable(1,N);
dQ1dt          = opti.variable(1,N); 
dQ2dt          = opti.variable(1,N); 

% xi

% vt = [0 -9 -2 2 5];
vt = 0;
Ft = interp1(vH, FMvtilda, vt);

opti.subject_to(Q0 > 0);
opti.subject_to(q > 0);

% extra constraints
opti.subject_to(Q1 - Q0 .* p == 0);
opti.subject_to(Q2 - Q0 .* (p.^2 + q) == 0);

% initial guess states
opti.subject_to(Q0(1,1) == xi(end,1));
opti.subject_to(Q1(1,1) == xi(end,2));
opti.subject_to(Q2(1,1) == xi(end,3));

opti.set_initial(Q0, xi(end,1)*ones(1,N));
opti.set_initial(Q1, xi(end,2)*ones(1,N));
opti.set_initial(Q2, xi(end,3)*ones(1,N));

opti.set_initial(p, (xi(end,2)/xi(end,1)) * ones(1,N));
opti.set_initial(q, (xi(end,3)/xi(end,1) - xi(end,2).^2/xi(end,1)) * ones(1,N));

opti.set_initial(dQ0dt, zeros(1,N));
opti.set_initial(dQ1dt, zeros(1,N));
opti.set_initial(dQ2dt, zeros(1,N));

if length(y0) > 3
    Non          = opti.variable(1,N);
    dNondt          = opti.variable(1,N); 
    opti.set_initial(Non, xi(end,4)*ones(1,N));
    opti.set_initial(dNondt, zeros(1,N));
    opti.subject_to(1e-4 < Non <= 1);
end

if length(y0) > 4
    DRX          = opti.variable(1,N);
    dDRXdt       = opti.variable(1,N); 
    opti.set_initial(DRX, xi(end,5)*ones(1,N));
    opti.set_initial(dDRXdt, zeros(1,N));
    opti.subject_to(.1 < DRX <= 1);
end

% parameters
allparms = {'f','k11','k12','k21','k22','JF','koop','J1','J2'};

for i = 1:length(allparms)
    eval([allparms{i}, ' = ', num2str(parms.(allparms{i}))])
end

optparms = {'k11','k21','k22'};
lb = [1 1 0 1];
ub = [2e3 2e3 5 1e3];

for i = 1:length(optparms)
    eval([optparms{i}, '= opti.variable(1)'])
    eval(['opti.subject_to(',num2str(lb(i)), '<', optparms{i}, '<', num2str(ub(i)),')']);
    eval(['opti.set_initial(',optparms{i},',', num2str(parms.(optparms{i})),')']);
end

% model constraints
Ns = floor(linspace(1, N, length(vt)+1));

% dynamic constraints
error = [];
dt1 = mean(diff(toc));

parms.Ca = 1;
parms.Noverlap = 1;

for i = 1:length(vt)
    
    id = Ns(i):Ns(i+1);
    F = (Q0(id) + Q1(id));

    [error_thin] = ThinEquilibrium(parms.Ca, Q0(id), Non(id), dNondt(id), parms.kon, parms.koff, koop, parms.Noverlap);      
    [error_thick] = ThickEquilibrium(Q0(id), F, DRX(id), dDRXdt(id), J1, J2, JF, parms.Noverlap);
    [error1] = MuscleEquilibrium_alt_v2(Q0(id), p(id), q(id), dQ0dt(id), dQ1dt(id), dQ2dt(id), f, k11, k12, k21, k22,  Non(id), vt(i), DRX(id));

    error = [error; error_thin(:); error_thick(:); error1(:)];
end

opti.subject_to((dNondt(Ns(1)+1:Ns(end)-1) + dNondt(Ns(1)+2:Ns(end)))*dt1/2 + Non(Ns(1)+1:Ns(end)-1) == Non(Ns(1)+2:Ns(end)));
opti.subject_to((dDRXdt(Ns(1)+1:Ns(end)-1) + dDRXdt(Ns(1)+2:Ns(end)))*dt1/2 + DRX(Ns(1)+1:Ns(end)-1) == DRX(Ns(1)+2:Ns(end)));
opti.subject_to((dQ0dt(Ns(1)+1:Ns(end)-1) + dQ0dt(Ns(1)+2:Ns(end)))*dt1/2 + Q0(Ns(1)+1:Ns(end)-1) == Q0(Ns(1)+2:Ns(end)));
opti.subject_to((dQ1dt(Ns(1)+1:Ns(end)-1) + dQ1dt(Ns(1)+2:Ns(end)))*dt1/2 + Q1(Ns(1)+1:Ns(end)-1) == Q1(Ns(1)+2:Ns(end)));
opti.subject_to((dQ2dt(Ns(1)+1:Ns(end)-1) + dQ2dt(Ns(1)+2:Ns(end)))*dt1/2 + Q2(Ns(1)+1:Ns(end)-1) == Q2(Ns(1)+2:Ns(end)));

% opti.subject_to((dNondt(Ns(i)+1:Ns(i+1)-1) + dNondt(Ns(i)+2:Ns(i+1)))*dt1/2 + Non(Ns(i)+1:Ns(i+1)-1) == Non(Ns(i)+2:Ns(i+1)));
% opti.subject_to((dDRXdt(Ns(i)+1:Ns(i+1)-1) + dDRXdt(Ns(i)+2:Ns(i+1)))*dt1/2 + DRX(Ns(i)+1:Ns(i+1)-1) == DRX(Ns(i)+2:Ns(i+1)));
% opti.subject_to((dQ0dt(Ns(i)+1:Ns(i+1)-1) + dQ0dt(Ns(i)+2:Ns(i+1)))*dt1/2 + Q0(Ns(i)+1:Ns(i+1)-1) == Q0(Ns(i)+2:Ns(i+1)));
% opti.subject_to((dQ1dt(Ns(i)+1:Ns(i+1)-1) + dQ1dt(Ns(i)+2:Ns(i+1)))*dt1/2 + Q1(Ns(i)+1:Ns(i+1)-1) == Q1(Ns(i)+2:Ns(i+1)));
% opti.subject_to((dQ2dt(Ns(i)+1:Ns(i+1)-1) + dQ2dt(Ns(i)+2:Ns(i+1)))*dt1/2 + Q1(Ns(i)+1:Ns(i+1)-1) == Q1(Ns(i)+2:Ns(i+1)));


% error = [error1 error2];
opti.subject_to(error == 0);

% cost function
J = 0;

F = opti.variable(1, length(vt));
for i = 1:length(vt)
    F(i) = (Q1(Ns(i+1)) + Q0(Ns(i+1)));
end

Frel = F ./ F(1);

for i = 1:length(vt)
    J = J + (Frel(i) - Ft(i))^2 ;
end

if length(y0) > 3
   J = J + 0.001 * (sumsqr(dQ0dt)+sumsqr(dQ1dt)+sumsqr(dQ2dt)+sumsqr(dNondt)+sumsqr(dDRXdt));
else
   J = J + 0.001 * (sumsqr(dQ0dt)+sumsqr(dQ1dt)+sumsqr(dQ2dt)); 
end

opti.minimize(J); 
    
%% Solve problem
% options for IPOPT
options.ipopt.tol = 1*10^(-6);          
options.ipopt.linear_solver = 'mumps';
% opti.solver('ipopt',options);

% Solve the OCP
p_opts = struct('expand',true);
s_opts = struct('max_iter', 1000);
opti.solver('ipopt',p_opts,s_opts);

sol = opti.solve();  


%% Get values
R.Q0 = sol.value(Q0); 
R.Q1 = sol.value(Q1); 
R.Q2 = sol.value(Q2); 
R.dQ0dt = sol.value(dQ0dt); 
R.dQ1dt = sol.value(dQ1dt); 
R.dQ2dt = sol.value(dQ2dt); 

R.F = sol.value(F);

R.p = sol.value(p);
R.q = sol.value(q);

R.t = 0:dt1:(N-1)*dt1;

% parameters
for i = 1:length(optparms)
    parms.(optparms{i}) = eval(['sol.value(',optparms{i},')']);
end

if length(y0)>3
    R.Non = sol.value(Non);
    R.dNondt = sol.value(dNondt);
    
    R.DRX = sol.value(DRX);
end

figure(1)
plot(R.t, [R.Q0; R.Q1; R.Q2])

return
%% Check whether p and q are okay
% close all
% subplot(121)
% plot(R.p); hold on
% plot(R.Q1./R.Q0, '--')
% 
% subplot(122)
% plot(R.q); hold on
% plot(R.Q2./R.Q0 - R.Q1.^2./R.Q0.^2, '--')

% Evaluate force-pCa
% close all
% figure(10)

n = nan(length(vs),2);
pt = nan(length(vs),2);

parms.vMtilda = 0;   

Cas = 10.^(-2:.1:1);
colors = parula(length(Cas));

states = {'Q_0','Q_1','Q_2','N_{on}','DRX', 'SRX'};

for i = 1:length(Cas)

    parms.Ca = Cas(i);

    [t0,x0] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 .1], y0, yp0);

    SRX = 1 - x0(:,1) - x0(:,end);
    
    X = [x0 SRX];
    
    figure(1)
    for j = 1:size(X,2)
        subplot(2,3,j)
        plot(t0, X(:,j), 'color', colors(i,:)); hold on
        title(states{j})
    end
    
    F0s(i) = x0(end,1) + x0(end,2);
 
end

%%
% close all

figure(2)
subplot(121)
semilogx(Cas, F0s); hold on
xline(1,'--')
% F0 = max(F0s);

% Evaluate force-velocity
vs = linspace(-10,5,100);
Fss = nan(length(vs),1);

parms.vMtilda = 0;
parms.Ca = 1;
parms.Noverlap = 1;

[t0,x0] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 1], y0, yp0);

F0 = (x0(end,1) + x0(end,2));


% if ishandle(2), close(2); end

acts = 1;

for j = 1:length(acts)
    parms.Noverlap = acts(j);
    
    for i = 1:length(vs)

    %     disp(vs(i))
        parms.vMtilda = vs(i);

        [ti,xi] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 1], x0(end,:), yp0);

        F = xi(:,1) + xi(:,2);

        Fss(i,j) = F(end);

    %     figure(2)
    %     plot(ti, F); hold on


    end
end

Fss_rel = Fss / F0;


% FMv  = e1*log((e2*vs./vmax+e3)+sqrt((e2*vs./vmax+e3).^2+1))+e4; % extensor
%%

% if ishandle(2), close(2); end

% figure(2)

subplot(122)
plot(vH, FMvtilda,'linewidth',2); hold on
plot(vs(:), Fss_rel(:,1),'--','linewidth',2); hold on
% plot(vs(:), Fss_rel(:,2),'--','linewidth',2); hold on
% plot(vs, parms.d * vs(:),':')
% plot(vt, Fs_rel, 'x')
% plot(vs, Fss(:),':')
% plot(ones(size(F)) * vs(end), F/F0, '.')

plot(vt, R.F/R.F(1), 's')

for i = 1:length(vt)
xline(vt(i),'k--')
end

legend('Hill','Biophysical','location','best')
legend boxoff
xlabel('Velocity (L_0/s)')
ylabel('Force (-)')
box off

%% test for different activation

close all
delta = 1/F0;
parms.vMtilda = 0;

acts = linspace(0,1,10);
overlaps = ones(1,10);

parms.f = 1000;

% acts = ones(1,20);
% overlaps = linspace(0,1,20);

F00 = nan(1,length(acts));
for j = 1:2
    if j == 1
        acts = linspace(0,1,10);
        overlaps = ones(1,10);
    else
       acts = ones(1,10);
        overlaps = linspace(0,1,10);
    end
    
for i = 1:length(acts)
    disp(i)
    
    parms.act = acts(i);
    parms.Noverlap = overlaps(i);
    
    [t0,x0] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 1], y0, yp0);

    SRX = 1 - x0(:,end);
    
    X = [x0 SRX];
    
    F = (x0(:,1) + x0(:,2)) * delta;
    
    F00(i) = F(end);
%     F0
    
    figure(1)
    subplot(1,2,j)
    plot(t0, F, 'color', colors(i,:)); hold on
end

if j == 1
    figure(2)
    subplot(121)
    plot(acts, F00);
else
    
    figure(2)
    subplot(122)
    plot(overlaps, F00);
end
end

return

%% Plot
% Simulate

ti = 0;

if length(y0)>3
    xi = [R.Q0(1) R.Q1(1) R.Q2(1) R.Non(1)];
else
    xi = [R.Q0(1) R.Q1(1) R.Q2(1)];
end

for i = 1:length(vt)
    parms.vMtilda = vt(i);
    [t1,x1] = ode15i(@(t,y,yp) OdeFun_FV(t,y,yp, parms), [0 max(toc)/length(vt)], [R.Q0(Ns(i)+1) R.Q1(Ns(i)+1) R.Q2(Ns(i)+1)]', [R.dQ0dt(Ns(i)+1) R.dQ1dt(Ns(i)+1) R.dQ2dt(Ns(i)+1)]');
    
    ti = [ti; t1+ti(end)];
    xi = [xi; x1];
end


%%
if ishandle(2), close(2); end
figure(2)

Q = [R.Q0; R.Q1; R.Q2]; 

for i = 1:3
    subplot(1,3,i)
    plot(toc, Q(i,:))

    hold on
    plot(ti, xi(:,i),'--')
end
