%% initialize 
clear; clc; 

addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
parms.forcible_detachment = 0;
parms.kse = parms.kse*0.25;
parms.no_tendon = 0;

parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.k11 = parms.k11 * 0.12;
parms.k12 = parms.k12 * 0.42;
parms.k21 = parms.k21 * 0.25;
parms.k22 = parms.k22 * 0.4;

% parms.xss = zeros(1,7);
% parms.n_func = @(xi, Q, eps) Q(1) ./ (sqrt(2*pi)*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))) * exp(-((xi-(Q(2)/max(Q(1), eps))).^2) / (2*(sqrt(max(Q(3)/Q(1) - (Q(2)/Q(1))^2, eps)))^2)); 

nbins = 1000;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0.0909;
parms.xss(end-1) = -40;
parms.xss(end) = -40;
parms.xi0 = linspace(-20,50,nbins);
parms.xi = parms.xi0;

parms0 = parms;

%% run sumulation for 2 pCa conditions x 2 pre-move conditions
figure
hold on

model = @fiber_dynamics;

for pCa_cond = [6.8] % <<-- change it to try different activation conditions
% try 4.5 (max) vs. 7.2 (medium) ... 9 minimum activation
    pCa = pCa_cond;
    Ca = 10^(-pCa+6);
    for k=1:2
        parms = parms0;
        x0 = parms.xss;

        % run prescribed motion / hold prior to knee drop 
        if(k==1)
            prescribed_a = spline([-1,0],[0,0]); % <<-- hold during first 1 sec
        else
            temp_t = -1:0.001:0;
            prescribed_a = spline(temp_t,-cos(temp_t/1*pi*2)*pi^2/2*pi*2);  % <<-- prescribed pre-movement.
        end
        [t_sim1,x_sim1] = ode15s(@(t,x)dAllStates_prescribed(t,x,model,parms,Ca, prescribed_a), [-1:0.005:0], ...
            [pi/2, 0, x0],odeopt);
        F_sim1 = nan(height(x_sim1),1);
        x0 = x_sim1(end,:); % << limb state at onset of drop 

        % run knee drop simulation 
        [t_sim,x_sim] = ode15s(@(t,x)dAllStates(t,x,model,parms,Ca), [0:0.005:10], ...
            x0,odeopt);
        
        plot([t_sim1;t_sim], [x_sim1(:,1);x_sim(:,1)]*180/pi-90, 'LineWidth', 3-k)

    end
end
xlim([-1 10])
legend('pCa = 8, no pre-move', 'pCa = 8, pre-move', 'pCa = 6.85, no pre-move', 'pCa = 6.85, pre-move')
xlabel('time (s)')

%% differential equation for prescribed movement 
function Xd = dAllStates_prescribed(t, X, muscle_model, parms, Ca, prescribed_a)
% X(1): mass position, X(2): mass velocity, X(3:end): muscle states
parms.vmtc = -X(2)*30;
[Xmusd,~] = muscle_model(t, X(3:end), parms, Ca);

% mass dynamics % 
Xd = [X(2); ppval(prescribed_a,t); Xmusd];
end

%% differential equation for unprescribed movement 
function Xd = dAllStates(t, X, muscle_model, parms, Ca)
% X(1): mass position, X(2): mass velocity, X(3:end): muscle states

parms.vmtc = -X(2)*30;
[Xmusd,F_mus] = muscle_model(t, X(3:end), parms, Ca);
F_ext = -9.81*sin(X(1))*3.5 - 0.8*X(2); 

% mass dynamics % 
Xd = [X(2); (F_ext+F_mus*380); Xmusd];
end