%% Initialize simulation 
addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;
% parms.forcible_detachment = 0;
% parms.kse = 0;
% parms.pse = 0;
% parms.no_tendon = 1;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.xss = zeros(1,7);
parms.n_func = @(xi,Q,eps)Q(1)...
    ./(sqrt(2*pi)*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps))))*...
    exp(-((xi-(Q(2)/max(Q(1),eps))).^2)...
    /(2*(sqrt(max(Q(3)/Q(1)-(Q(2)/Q(1))^2,eps)))^2));
% 
% nbins = 500;
% parms.nbins = nbins;
% parms.xss = zeros(1,parms.nbins + 4);
% parms.xi0 = linspace(-20,20,nbins);
% parms.xi = parms.xi0;

parms.xss(end-2) = 0.0909;
parms.xss(end-1) = -20;
parms.xss(end) = -20;

% have two identital muscles as agonist-antagonist pair 
parms1 = parms;
parms2 = parms;

parms1.vRatio = -60;
parms2.vRatio = 60;

model = @fiber_dynamics;

%% define cart perturbation protocol 
temp_t = -1:0.001:1;
temp_acc = zeros(size(temp_t));

perturb_num = 10;
temp_acc( (0:perturb_num)+1000) = -cos((0:perturb_num)*2*pi/perturb_num)+1;
temp_acc( (0:perturb_num)+1500) = cos((0:perturb_num)*2*pi/perturb_num)-1;
cart_acc_spline = spline(temp_t,temp_acc*50);

%% run simulation 
figure; hold on
for muscle_itr = 0:1

    if(muscle_itr==0) % free fall - no muscle 
        [t_sim,x_sim] = ode45(@(t,x)dMassStates(t,x,cart_acc_spline), ...
            [-1:0.001:0.6], ... simulation time
            [0, 0], ... initial condition
            odeopt);
    else % simulation with a muscle pair 
        pCa = 6.8; % <<- change here to try different baseline activation 
        Ca = 10^(-pCa+6);
        num_Mstate = length(parms.xss);

        [t_sim,x_sim] = ode15s(@(t,x)dAllStates(t,x,model,parms1,parms2,...
            Ca,Ca,num_Mstate,cart_acc_spline), ...
            [-1:0.001:0.6], ... simulation time
            [0, 0, parms1.xss, parms2.xss], ... initial condition
            odeopt);

        % calculate force 
        F_sim = nan(height(x_sim),2);
        for i = 1:height(x_sim)
            parms1.vmtc = x_sim(i,2)*parms1.vRatio;
            parms2.vmtc = x_sim(i,2)*parms2.vRatio;
            
            [~,F_sim(i,1)] = ...
                model(t_sim(i), x_sim(i,3:num_Mstate+2)', parms1, Ca);
            [~,F_sim(i,2)] = ...
                model(t_sim(i), x_sim(i,num_Mstate+3:end)', parms2, Ca);
        end

        subplot(3,1,2)
        plot(t_sim, F_sim/(F_sim(1000-1,1)), 'color', [muscle_itr 0 1-muscle_itr])
        hold on

        subplot(3,1,3)
        plot(t_sim, x_sim(:,end)/half_s_len_norm, 'color', [muscle_itr 0 1-muscle_itr])
        hold on
    end

    subplot(3,1,1)
    plot(t_sim, x_sim(:,1)*180/pi, 'color', [muscle_itr 0 1-muscle_itr])
    hold on
end

linkaxes(get(gcf,'children'), 'x')
xlabel('time (s)')
xlim([-0.1 0.6])

subplot(3,1,1)
ylabel('angle (deg)')

subplot(3,1,2)
ylabel('muscle force (A.U.)')

subplot(3,1,3)
ylabel('muscle length')
xlabel('time (s)')

%% differential equation with two muscles 
function Xd = dAllStates(t,X, muscle_model, parms1, parms2, ...
    Ca1, Ca2, num_Mstate, cart_acc_spline)
% X(1): mass position, X(2): mass velocity, 
% X(3:num_Mstate+2): muscle 1 states, X(num_Mstate+3:end): muscle 2 states 

parms1.vmtc = X(2)*parms1.vRatio;
[Xmusd1,F_mus1] = muscle_model(t, X(3:num_Mstate+2), parms1, Ca1);

parms2.vmtc = X(2)*parms2.vRatio;
[Xmusd2,F_mus2] = muscle_model(t, X(num_Mstate+3:end), parms2, Ca2);

Xmassd = dMassStates(t,X,cart_acc_spline);

Xd = [Xmassd+[0; (F_mus1-F_mus2)*200]; ...
    Xmusd1; Xmusd2]; 
end

%% differential equation for the dybnamics of the mass 
function Md = dMassStates(t,X, cart_acc_spline)
cart_acc = ppval(cart_acc_spline,t);                              
acc_gravity = +9.81*sin(X(1));
Md = [X(2); acc_gravity+cart_acc*cos(X(1))];
end