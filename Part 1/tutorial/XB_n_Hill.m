clear;
addpath(genpath([pwd,'/../..']))
warning('off')
load('parms.mat')
load('protocol.mat')
warning('on')

half_s_len_norm = parms.s/2/parms.h;
nbins = 600;
odeopt = odeset('maxstep',1e-2);

parms.forcible_detachment = 0;
parms.kse = 0;
parms.kpe = 0;
parms.no_tendon = 1;
parms.act = 1;
parms.cosa = 1;
parms.Noverlap = 1;

parms.xi0 = linspace(-10,10,nbins);
parms.xi = parms.xi0;
parms.nbins = nbins;
parms.xss = zeros(1,parms.nbins + 4);
parms.xss(end-2) = 0.0909;
activation_pCa = 4.5;

model = @fiber_dynamics;
parms0 = parms;

load('hill_properties_pCa_45.mat');
fig.UserData.hill_properties = hill_properties;

% define your protocol. 
% make sure that duration, velocity, pCa all has same number of elements. 
protocol_duration = [0.05, 0.1, 0.1, 0.1, 0.1,0.1]; % should be integer multiplication of dt
protocol_pCa = [9 [1 1 1 1 1]*activation_pCa];

% protocol_v = [0, 0,  0, 50, 0, 0]; % lengthening
% protocol_v = [0, 0, 0, -50, 0, 0]; % shortening
% protocol_v = [0, 0, -50, 50, -50, 50]; % shortening cycle
protocol_v = [0, 0, 50, -50, 50, -50]; % lengthening cycle


x0 = parms.xss;
dt = 0.001;


t_total = nan(round(sum(protocol_duration)/dt)+length(protocol_duration),1);
x_total = nan(length(t_total),length(x0));
v_total = nan(size(t_total));
F_total = nan(size(t_total));
pCa_total = nan(size(t_total));

plot_t_protocol = nan(length(protocol_pCa)*2, 1);
plot_pCa_protocol = nan(length(protocol_pCa)*2, 1);

curr_idx = 1;
curr_time = 0;
for p_itr = 1:length(protocol_v)
    plot_t_protocol((-1:0)+p_itr*2) = [curr_time, curr_time+protocol_duration(p_itr)];
    plot_pCa_protocol((-1:0)+p_itr*2) = protocol_pCa(p_itr);

    parms.vmtc = protocol_v(p_itr);
    Ca = 10.^(-protocol_pCa(p_itr)+6);

    [t_sim,x_sim] = ode15s(model, ...
        curr_time:dt:curr_time+protocol_duration(p_itr), x0, odeopt, ...
        parms, Ca);
    F_sim = nan(size(t_sim));
    for k_itr = 1:length(t_sim)
        [~,F_sim(k_itr)] = model(t_sim(k_itr), x_sim(k_itr,:)', parms, Ca);
    end
    t_total(curr_idx:curr_idx+length(t_sim)-1) = t_sim;
    pCa_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_pCa(p_itr);
    v_total(curr_idx:curr_idx+length(t_sim)-1) = protocol_v(p_itr);

    F_total(curr_idx:curr_idx+length(t_sim)-1) = F_sim;
    x_total(curr_idx:curr_idx+length(t_sim)-1,:) = x_sim;
    x0 = x_sim(end,:);
    curr_idx = curr_idx + length(t_sim);
    curr_time = t_sim(end);
end

%%
figure
subplot(3,1,1)
plot(plot_t_protocol, plot_pCa_protocol)
ylabel('pCa')
subplot(3,1,2)
plot(t_total, x_total(:,end))
ylabel('length')
subplot(3,1,3)
plot(t_total, F_total)
xlabel('time')
ylabel('Force')

%%
figure 
xi = repmat(parms.xi0,height(x_total),1)+...
    repmat(x_total(:,end), 1, width(parms.xi0));
plot(xi(1:10:end,:)', x_total(1:10:end,2:end-3)')
set(gca, 'colorOrder', winter(round(height(x_total)/10)+1))
xlabel('')