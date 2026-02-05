%% Initialize simulation 
addpath(genpath([pwd,'/..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

odeopt = odeset('maxstep',1e-2);
half_s_len_norm = parms.s/2/parms.h;
nbins = 500;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.no_tendon = 0;
% parms.kT = 0.000100;
parms.kse = parms.kse*1000;
parms.kpe = 0;

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

parms.xss(end-2) = 0.0909;
parms.xss(end-1) = -20; % initial point on F-l curve
parms.xss(end) = -20;

% have two identital muscles as agonist-antagonist pair 
parms1 = parms;
parms2 = parms;

parms1.vRatio = -60; % dmuscle_len / dtheta
parms2.vRatio = 60;

model = @fiber_dynamics;

%% define cart perturbation protocol 
temp_t = -1:0.001:1;
temp_acc = zeros(size(temp_t));

perturb_num = 50;
temp_acc( (0:perturb_num)+1000) = -cos((0:perturb_num)*2*pi/perturb_num)+1;
temp_acc( (0:perturb_num)+1500) = cos((0:perturb_num)*2*pi/perturb_num)-1;
cart_acc_spline = spline(temp_t,temp_acc*10);

%% run simulation 
figure; hold on
for muscle_itr = 0:1

    if(muscle_itr==0) % free fall - no muscle 
        [t_sim,x_sim] = ode45(@(t,x)dMassStates(t,x,cart_acc_spline,0), ...
            [-1:0.001:0.6], ... simulation time
            [0, 0], ... initial condition
            odeopt);
    else % simulation with a muscle pair 


        for kkkkkkkkkkk=1:10
        pCa0 = 6.8; % <<- change here to try different baseline activation 
        pCa1 = pCa0;
        pCa2 = pCa0;

        Ca1 = 10^(-pCa1+6);
        Ca2 = 10^(-pCa2+6);

        num_Mstate = length(parms.xss);

        time_step = 0.001;
        t_sim = -1:time_step:3;
        
        x_t = [0, 0, parms1.xss, parms2.xss];

        x_sim = nan(length(t_sim), length(x_t));
        F_sim = nan(length(t_sim), 2);
        muscle_sense = nan(length(t_sim), 3, 2);
        pCa_sim = nan(length(t_sim), 2);

        muscle_sense(1,:,:) = 0;
        pCa_sim(1,:) = pCa0;
        
        muscleF_acc = 0;
        delay_num = round(0.1/time_step);

        tic
        for ode_itr = 1:length(t_sim)
            xdot = dAllStates(t_sim(ode_itr),x_t,model,parms1,parms2,...
                Ca1,Ca2,num_Mstate,cart_acc_spline,muscleF_acc);

            x_sim(ode_itr,:) = x_t;
            [~,F_sim(ode_itr,1)] = ...
                model(t_sim(ode_itr), x_sim(ode_itr,3:num_Mstate+2)', parms1, Ca1);
            [~,F_sim(ode_itr,2)] = ...
                model(t_sim(ode_itr), x_sim(ode_itr,num_Mstate+3:end)', parms2, Ca2);

            muscleF_acc = (F_sim(ode_itr,1) - F_sim(ode_itr,2))*1000;
                        
            muscle_sense(ode_itr+1, :, 1) = ...
                [x_t(2+num_Mstate), xdot(2+num_Mstate), ...
                (xdot(2+num_Mstate)-muscle_sense(ode_itr,2,1))/time_step ];
            muscle_sense(ode_itr+1, :, 2) = ...
                [x_t(2+num_Mstate*2), xdot(2+num_Mstate*2),...
                (xdot(2+num_Mstate*2)-muscle_sense(ode_itr,2,2))/time_step];
            
            x_t = x_t + time_step*xdot';

            if(parms.no_tendon)
                K = [500 40 1]'/60000;
            else
                K = [2800 70 0.5]'/60000;
            end

            if(ode_itr>1000)
                num_diff_a1 = [1 -8 8 -1]*muscle_sense(ode_itr-delay_num+[-2 -1 1 2], 2, 1)/12; % improve numerical differentiation 
                num_diff_a2 = [1 -8 8 -1]*muscle_sense(ode_itr-delay_num+[-2 -1 1 2], 2, 2)/12;

                muscle_sense(ode_itr-delay_num, 3, 1) = num_diff_a1;
                muscle_sense(ode_itr-delay_num, 3, 2) = num_diff_a2;

                pCa1 = pCa0 - max(min((muscle_sense(ode_itr-delay_num, :, 1)-muscle_sense(1000, :, 1))*K, pCa0-4.5),0);
                pCa2 = pCa0 - max(min((muscle_sense(ode_itr-delay_num, :, 2)-muscle_sense(1000, :, 2))*K, pCa0-4.5),0);
            end
            Ca1 = 10^(-pCa1+6);
            Ca2 = 10^(-pCa2+6);

            pCa_sim(ode_itr,:) = [pCa1, pCa2];

            if(ode_itr==length(t_sim)), break; end
            if(abs(x_t(1))>pi*3), break; end
            % final iteration is only to get final XB state, so stop before
            % assigning ySim ySim(ode_itr+1,:) = [x_t' xdot(2)];
        end
        toc
        end

        subplot(5,1,2)
        plot(t_sim, F_sim/(F_sim(1000,1)))
        hold on

        subplot(5,1,3)
        plot(t_sim, x_sim(:,[-num_Mstate 0]+end)/half_s_len_norm)
        hold on

        subplot(5,1,4)
        plot(t_sim, (x_sim(:,[-num_Mstate 0]+end-1))/half_s_len_norm-0.016)
        hold on

        subplot(5,1,5)
        plot(t_sim, pCa_sim)
        hold on
    end

    subplot(5,1,1)
    plot(t_sim, x_sim(:,1)*180/pi, 'color', [muscle_itr 0 1-muscle_itr])
    hold on
end

linkaxes(get(gcf,'children'), 'x')
xlabel('time (s)')
xlim([-0.1 3])

subplot(5,1,1)
ylabel('angle (deg)')

subplot(5,1,2)
ylabel('muscle force (A.U.)')

subplot(5,1,3)
ylabel('muscle length')
xlabel('time (s)')

%% differential equation with two muscles 
function Xd = dAllStates(t,X, muscle_model, parms1, parms2, ...
    Ca1, Ca2, num_Mstate, cart_acc_spline, muscleF_acc)
% X(1): mass position, X(2): mass velocity, 
% X(3:num_Mstate+2): muscle 1 states, X(num_Mstate+3:end): muscle 2 states 

parms1.vmtc = X(2)*parms1.vRatio;
[Xmusd1,F_mus1] = muscle_model(t, X(3:num_Mstate+2), parms1, Ca1);

parms2.vmtc = X(2)*parms2.vRatio;
[Xmusd2,F_mus2] = muscle_model(t, X(num_Mstate+3:end), parms2, Ca2);

Xmassd = dMassStates(t,X,cart_acc_spline, muscleF_acc);

Xd = [Xmassd+[0; (F_mus1-F_mus2)*200]; ...
    Xmusd1; Xmusd2]; 
end

%% differential equation for the dynamics of the mass 
function Md = dMassStates(t,X,cart_acc_spline, muscleF_acc)
cart_acc = ppval(cart_acc_spline,t);                              
gravity_acc = +9.81*sin(X(1));
Md = [X(2); gravity_acc+cart_acc*cos(X(1))+muscleF_acc];
end