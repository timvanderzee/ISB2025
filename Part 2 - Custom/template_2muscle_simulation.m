% clear; clc; close all

addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;
parms.forcible_detachment = 0;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

% approximation model %
parms.xss = zeros(1,7);
parms.n_func = @(xi,Q,eps)Q(1)...
    ./(sqrt(2*pi)*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps))))*...
    exp(-((xi-(Q(2)/max(Q(1),eps))).^2)...
    /(2*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps)))^2));
 
% discrete bin model -- can be used instead of "approximation model" % 
% nbins = 500;
% parms.nbins = nbins;
% parms.xss = zeros(1,parms.nbins + 4);
% parms.xi0 = linspace(-20,20,nbins);
% parms.xi = parms.xi0;

% have two identital muscles as agonist-antagonist pair 
parms1 = parms;
parms2 = parms;

%%
model = @fiber_dynamics;
pCa = 6.8;
Ca = 10^(-pCa+6);
num_Mstate = length(parms.xss);

mass_IC = [0.1, 0]; % initial condition of the mass

[t_sim,x_sim] = ode15s(@(t,x)dAllStates(t,x,model,parms1,parms2,...
    Ca,Ca,num_Mstate), ...
    [-1:0.001:0.6], ... simulation time
    [[mass_IC], parms1.xss, parms2.xss], ... initial condition
    odeopt);

F_sim = nan(height(x_sim),2);
for i = 1:height(x_sim)
    [~,F_sim(i,1)] = ...
        model(t_sim(i), x_sim(i,3:num_Mstate+2)', parms1, Ca);
    [~,F_sim(i,2)] = ...
        model(t_sim(i), x_sim(i,num_Mstate+3:end)', parms2, Ca);
end


subplot(3,1,1)
plot(t_sim, x_sim(:,1)*180/pi)
ylabel('angle (degree)')
subplot(3,1,2)
plot(t_sim, F_sim/(F_sim(1000-1,1)))
ylabel('muscle force')
subplot(3,1,3)
plot(t_sim, x_sim(:,end-num_Mstate)/half_s_len_norm, t_sim, x_sim(:,end)/half_s_len_norm)
legend('muscle 1', 'muscle 2')
ylabel('muscle length')
xlabel('time')
linkaxes(get(gcf,'children'), 'x')

%%

function Xd = dAllStates(t,X, muscle_model, parms1, parms2, ...
    Ca1, Ca2, num_Mstate)
% X(1): segment position, X(2): segment velocity, 
% X(3:num_Mstate+2): muscle 1 states, 
% X(num_Mstate+3:end): muscle 2 states 

parms1.vmtc = -X(2)*60; % <<- segment-muscle1 constraint 
[Xmusd1,F_mus1] = muscle_model(t, X(3:num_Mstate+2), parms1, Ca1);

parms2.vmtc = X(2)*60;  % <<- segment-muscle2 constraint 
[Xmusd2,F_mus2] = muscle_model(t, X(num_Mstate+3:end), parms2, Ca2);

Xmassd = dsegment(t,X);

Xd = [Xmassd+[0; (F_mus1-F_mus2)*200]; ...
    Xmusd1; Xmusd2]; 
end

function Md = dsegment(t,X)
% define passive dynamics of the mass
% example below is for an inverted pendulum. 
accleration = +9.81*sin(X(1));
Md = [X(2); accleration];
end